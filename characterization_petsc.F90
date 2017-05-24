!=============================================================================80
! PETSc take on reservoir fluid characterization
! Built from PETSc SNES Fortran example ex1f.F
!  Description: Uses the Newton method to solve a two-variable system.
!
!/*T
!  Concepts: SNES^basic uniprocessor example
!  Processors: 1
!T*/
!
! -----------------------------------------------------------------------

PROGRAM main
  USE roots
  IMPLICIT NONE


  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !                    Include files
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  !  The following include statements are generally used in SNES Fortran
  !  programs:
  !     petscsys.h       - base PETSc routines
  !     petscvec.h    - vectors
  !     petscmat.h    - matrices
  !     petscksp.h    - Krylov subspace methods
  !     petscpc.h     - preconditioners
  !     petscsnes.h   - SNES interface
  !  Other include statements may be needed if using additional PETSc
  !  routines in a Fortran program, e.g.,
  !     petscviewer.h - viewers
  !     petscis.h     - index sets
  !
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscsnes.h>
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !                   Variable declarations
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  !  Variables:
  !     snes        - nonlinear solver
  !     ksp        - linear solver
  !     pc          - preconditioner context
  !     ksp         - Krylov subspace method context
  !     x, r        - solution, residual vectors
  !     J           - Jacobian matrix
  !     its         - iterations for convergence

  !PetscFortranAddr ctx(6)
  ! Context Delivery Variables
  Vec      Vec_ctx
  PetscInt ctx_size
  
  ! Main Variable Collection
  SNES     snes
  PC       pc

  KSP      ksp
  Vec      Vec_x,Vec_r
  Mat      Mat_J
  SNESLineSearch linesearch
  PetscErrorCode  ierr
  PetscInt its,n,divtol
  PetscMPIInt size,rank
  PetscReal   tol
  PetscBool   setls
!  INTEGER :: i
  ! Program input variables
  LOGICAL :: l
  CHARACTER(80) :: input_file
  CHARACTER(80),ALLOCATABLE :: compname(:)
  INTEGER :: nc, nphase
  REAL(KIND=8) :: P,T,feed,Pc_c7pi,Tc_c7pi,MW_c7pi,rho_c7pi,               &
       m_c7pi(1),omega_c7pi(1)

  ! EOS (SRK) Related Variables
  REAL(KIND=8), ALLOCATABLE :: Tcr(:),Pcr(:),omega(:),MW(:),rho(:),z(:),   &
       x(:),y(:),m(:),coeff(:),vol_shift(:)

  ! Characterization Variables
  PetscInt cmax, cmin
  INTEGER :: lastTBProw
  REAL(kind=8) :: z_cnplus,rho_cnminus,MW_cnplus,x_init(4),rho_cnplus

  ! Solution Related Variables
  PetscInt  ::  init_ind(4)
  






  !  Note: Any user-defined Fortran routines (such as FormJacobian)
  !  MUST be declared as external.

  EXTERNAL FormFunction, FormJacobian, MyLineSearch


  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !                 Beginning of program
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CALL MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
  CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
  IF (size .NE. 1) THEN
     IF (rank .EQ. 0) THEN
        WRITE(6,*) 'This is a uniprocessor example only!'
     ENDIF
     SETERRQ(PETSC_COMM_SELF,1,' ',ierr)
  ENDIF
  !-----------------------------------------------------------------------------80
  ! Read Command Line Arguments
  !-----------------------------------------------------------------------------80
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=6)
  CALL PetscOptionsGetString(Petsc_Null_Character,'-f',input_file,l,ierr)
#else
  CALL PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-f',      &
       input_file,l,ierr)
#endif  
  IF (.NOT. l) THEN
     CALL PrintMsg("Usage: [mpiexec -n <np>] PTE511_HW* -f <input_filename>")
     GOTO 9
  END IF

  !=============================================================================80
  ! Read in sim parameters
  !-----------------------------------------------------------------------------80
  CALL PrintMsg("Reading input ...")
  OPEN(10,file=input_file,status='old')        
  ! Read # phases, # components
  READ(10,*) nc, nphase

  ALLOCATE(Tcr(nc),Pcr(nc),omega(nc),MW(nc),rho(nc),z(nc),x(nc),y(nc),m(nc),   &
       compname(nc),coeff(13),vol_shift(11))

  CALL ReadParameters
  CLOSE(10)

  ! Load literature referenced values (N2, CO2, H2S, C1 to C6)
  CALL SRK_Coefficients
  CALL PropstoC6 
  !-----------------------------------------------------------------------------80

  !=============================================================================80
  ! Assign stuff that should be read in
  !-----------------------------------------------------------------------------80
  lastTBProw = 15
  cmin = 11
  cmax = 200

  ! PETSc Stuff
  n = 4
  tol = 1.0E-4
  divtol = 20
  ctx_size = 6
  
  ! Initial Values for Solver
  x_init(1) = -3.0  ! A
  x_init(2) = -0.1  ! B
  x_init(3) =  0.5  ! C
  x_init(3) =  0.1  ! D
  ! Initial Values for Solver / index
  init_ind(1) = 0
  init_ind(2) = 1
  init_ind(3) = 2
  init_ind(4) = 3      

  !=============================================================================80
  ! Flesh out non PETSc property estimates
  !-----------------------------------------------------------------------------80 
  ! Initial C7+ Estimates
  MW_c7pi  = SUM(z(12:)*MW(12:)) / SUM(z(12:))
  rho_c7pi = SUM(z(12:)*MW(12:)) / SUM( (z(12:)*MW(12:))/rho(12:) )
  Tc_c7pi  = coeff(1)*rho_c7pi + coeff(2)*LOG(MW_c7pi) + coeff(3)*MW_c7pi +    &
       coeff(4)/MW_c7pi
  Pc_c7pi  = EXP(coeff(5) + coeff(6)*(rho_c7pi**coeff(9)) + coeff(7)/MW_c7pi + &
       coeff(8)/(MW_c7pi**2))

  ! Find omega via roots
  m_c7pi   = coeff(10) + coeff(11)*MW_c7pi + coeff(12)*rho_c7pi +              &
       coeff(13)*(MW_c7pi**2)
  CALL omegasolve(1,m_c7pi,omega_c7pi)
  !       
  ! Estimates to C20/whatever

  Tcr(12:) = coeff(1)*rho(12:) + coeff(2)*LOG(MW(12:)) + coeff(3)*MW(12:) +    &
       coeff(4)/(MW(12:)**2)
  Pcr(12:) = EXP(coeff(5) + coeff(6)*(rho(12:)**coeff(9)) + coeff(7)/MW(12:) + &
       coeff(8)/(MW(12:)**2))
  m(12:)   = coeff(10) + coeff(11)*MW(12:) + coeff(12)*rho(12:) +              &
       coeff(13)*(MW(12:)**2)

  ! Solve for omega values. SRK coefficients for quadratic equation
  CALL omegasolve(nc-12+1,m(12:),omega(12:))
  !---------------------------------------------------------------------------80

  !===========================================================================80
  ! Load context array
  !---------------------------------------------------------------------------80
  CALL VecCreateSeq(PETSC_COMM_SELF,ctx_size,Vec_ctx,ierr)
  
  ! 1) z from after last TBP cut on
  z_cnplus = SUM(z(lastTBProw:))
  CALL VecSetValue(Vec_ctx,0,z_cnplus,INSERT_VALUES,ierr)
!  ctx(1) = z_cnplus

  ! 2) rho from last TPB cut
  rho_cnminus = rho(lastTBProw)
  CALL VecSetValue(Vec_ctx,1,rho_cnminus,INSERT_VALUES,ierr)
 ! ctx(2) = rho_cnminus
 
  ! 3) Cmin and 4) Cmax
  CALL VecSetValue(Vec_ctx,2,cmin,INSERT_VALUES,ierr)
  CALL VecSetValue(Vec_ctx,3,cmax,INSERT_VALUES,ierr)  
!  ctx(3) = cmin; ctx(4)=cmax

  ! 5) MW from after last TBP cut on
  MW_cnplus = SUM(z(lastTBProw:)*MW(lastTBProw:)) / SUM(z(lastTBProw:))
  CALL VecSetValue(Vec_ctx,4,MW_cnplus,INSERT_VALUES,ierr)
!  ctx(5) = MW_cnplus

  ! 6) rho from after last TBP cut on
  rho_cnplus = SUM(z(lastTBProw:)*MW(lastTBProw:)) /                               &
           SUM( (z(lastTBProw:)*MW(lastTBProw:))/rho(lastTBProw:) )
  CALL VecSetValue(Vec_ctx,4,rho_cnplus,INSERT_VALUES,ierr)
  CALL VecAssemblyBegin(Vec_ctx,ierr)
  CALL VecAssemblyEnd(Vec_ctx,ierr)
!  ctx(6) = rho_cnplus  
  !===========================================================================80
  ! Detemine Parameters A, B, C, D
  !---------------------------------------------------------------------------80
  ! Process begins after last TBP HC, from C11+ in this case.

  ! Create Solution/NL Function vectors
  CALL SNESCreate(PETSC_COMM_WORLD,snes,ierr)
  CALL VecCreateSeq(PETSC_COMM_SELF,n,Vec_x,ierr)
  CALL VecDuplicate(Vec_x,Vec_r,ierr)

  ! Create Jacobian Matrix
  CALL MatCreate(PETSC_COMM_SELF,Mat_J,ierr)
  CALL MatSetSizes(Mat_J,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
  CALL MatSetFromOptions(Mat_J,ierr)
  CALL MatSetUp(Mat_J,ierr)

  ! Set function evaluation routine and vector
  CALL SNESSetFunction(snes,Vec_r,FormFunction,Vec_ctx,ierr)

  ! Set Jacobian matrix data structure and evaluation routine
  CALL SNESSetJacobian(snes,Mat_J,Mat_J,FormJacobian,Vec_ctx,ierr)
  !---------------------------------------------------------------------------80

  !===========================================================================80

  !===========================================================================80
  ! Customize and run nonlinear solver
  !---------------------------------------------------------------------------80
  CALL SNESGetKSP(snes,ksp,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCNONE,ierr)
  CALL KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,divtol,  &
       ierr)

  CALL SNESSetFromOptions(snes,ierr)
  
  call PetscOptionsHasName(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-setls',    &
  setls,ierr)

  if (setls) then
      call SNESGetLineSearch(snes, linesearch, ierr)
      call SNESLineSearchSetType(linesearch, 'shell', ierr)
      call SNESLineSearchShellSetUserFunc(linesearch, MyLineSearch,   &
     &PETSC_NULL_OBJECT, ierr)
  endif
  
  
  CALL PrintMsg("Initial Guess")
  ! Feed in initial guesses
  call VecSetValues(Vec_x,n,init_ind,x_init,INSERT_VALUES,ierr)
  CALL VecAssemblyBegin(Vec_x,ierr)
  CALL VecAssemblyEnd(Vec_x,ierr)
  CALL PrintMsg("Solving")  
  call SNESSolve(snes,PETSC_NULL_OBJECT,x,ierr)
  call SNESGetIterationNumber(snes,its,ierr);
  if (rank .eq. 0) then
     write(6,100) its
  endif
100 format('Number of SNES iterations = ',i5)

  CALL PrintMsg("Finished")
  call VecDestroy(Vec_x,ierr)
  call VecDestroy(Vec_r,ierr)
  call MatDestroy(Mat_J,ierr)
  call SNESDestroy(snes,ierr)
  call PetscFinalize(ierr)  
9 CALL PetscFinalize(ierr)  

  !=============================================================================80
CONTAINS
  !=============================================================================80
  SUBROUTINE ReadParameters
    IMPLICIT NONE
    INTEGER :: i
    ! Open input file
    ! Read separation P and T, plus inflow rate
    ! P, bar; T, K; feed, mol
    READ(10,*) P, T, feed
    ! Read main table components
    DO i=1,nc
       ! z, frac total; MW, g/mol; rho, g/cc, y, frac gas;
       READ(10,*) compname(i),z(i),MW(i),rho(i),y(i)
    END DO
    ! Convert values to correct units (P ==> atm)
    P = P*0.98692327
  END SUBROUTINE ReadParameters
  !-----------------------------------------------------------------------------80

  !=============================================================================80
  SUBROUTINE SRK_Coefficients
    IMPLICIT NONE
    coeff(1)    =   304.143
    coeff(2)    =   48.4052
    coeff(3)    =   0.710774
    coeff(4)    =   3800.73
    coeff(5)    =   3.05081
    coeff(6)    =   -0.903352
    coeff(7)    =   233.768
    coeff(8)    =   -12715.4
    coeff(9)    =   0.25
    coeff(10)   =   0.496902
    coeff(11)   =   0.00558442
    coeff(12)   =   0.010564
    coeff(13)   =   -5.243E-6
  END SUBROUTINE SRK_Coefficients
  !-----------------------------------------------------------------------------80      
  SUBROUTINE PropstoC6
    IMPLICIT NONE
    ! MW, g/mol      ; Pc, atm        ; Tc, K         
    MW(1)  = 28.016  ; Pcr(1)  = 33.6  ; Tcr(1)  = 126.2; omega(1)  = 0.040 ! N2
    MW(2)  = 44.01   ; Pcr(2)  = 72.9  ; Tcr(2)  = 304.2; omega(2)  = 0.228 ! CO2
    MW(3)  = 34.08   ; Pcr(3)  = 88.9  ; Tcr(3)  = 373.2; omega(3)  = 0.081 ! H2S
    MW(4)  = 16.043  ; Pcr(4)  = 45.4  ; Tcr(4)  = 190.4; omega(4)  = 0.008 ! C1        
    MW(5)  = 30.069  ; Pcr(5)  = 48.2  ; Tcr(5)  = 305.4; omega(5)  = 0.098 ! C2      
    MW(6)  = 44.096  ; Pcr(6)  = 41.9  ; Tcr(6)  = 369.8; omega(6)  = 0.152 ! C3
    MW(7)  = 58.123  ; Pcr(7)  = 36.0  ; Tcr(7)  = 408.2; omega(7)  = 0.176 ! iC4
    MW(8)  = 58.123  ; Pcr(8)  = 37.5  ; Tcr(8)  = 425.2; omega(8)  = 0.193 ! nC4
    MW(9)  = 72.15   ; Pcr(9)  = 33.4  ; Tcr(9)  = 460.4; omega(9)  = 0.227 ! iC5
    MW(10) = 72.15   ; Pcr(10) = 33.3  ; Tcr(10) = 469.7; omega(10) = 0.251 ! nC5
    MW(11) = 86.177  ; Pcr(11) = 29.3  ; Tcr(11) = 507.5; omega(11) = 0.296 ! C6
  END SUBROUTINE PropstoC6
  !-----------------------------------------------------------------------------80

  ! Print message ==============================================================80
  SUBROUTINE PrintMsg(msg)
    IMPLICIT NONE
    CHARACTER(*) :: msg
    IF (rank==0) PRINT *,msg
  END SUBROUTINE PrintMsg
  !-----------------------------------------------------------------------------80   

  ! Gateway to quadratic root finder ===========================================80
  SUBROUTINE omegasolve(neq,m,omega)
    IMPLICIT NONE
    INTEGER :: neq
    REAL(KIND=8) :: m(neq),omega(neq)
    REAL(KIND=8), ALLOCATABLE :: root_coeff(:,:)
    INTEGER, ALLOCATABLE :: root_type(:)

    ALLOCATE( root_coeff(neq,2),root_type(neq) )
    root_coeff(:,1) = 1.574 / (-0.176)
    root_coeff(:,2) = (0.480 / (-0.176)) 
    root_coeff(:,2) = root_coeff(:,2) - m(:) / (-0.176)
    root_type(:) = 0.0
    CALL poleq2(neq,root_coeff,root_type)
    omega(:) = MINVAL(root_coeff,2)
    DEALLOCATE( root_coeff, root_type )
  END SUBROUTINE omegasolve
END PROGRAM main
! The subroutines below are external (e.g. not part of program main)
!=============================================================================80
!     ------------------------------------------------------------------------
!     
!     FormFunction - Evaluates nonlinear function, F(x).
!     
!     Input Parameters:
!     snes - the SNES context
!     x - input vector
!     dummy - optional user-defined context (not used here)
!     
!     Output Parameter:
!     f - function vector
!     
      SUBROUTINE FormFunction(snes,Vec_x,Vec_f,Vec_ctx,ierr)
      IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscsnes.h>

      SNES     snes
      Vec      Vec_x,Vec_f,Vec_ctx
      PetscErrorCode ierr
      PetscScalar, pointer :: pntr(:)
  ! ctx
  ! 1) z from after last TBP cut on             
  ! 2) rho from last TPB cut
  ! 3) Cmin
  ! 4) Cmax
  ! 5) MW from after last TBP cut on  
  ! 6) rho from after last TBP cut on

!  Declarations for use with local arrays

! Local Variables     
      real(kind=8), allocatable :: z_rs(:),MW_rs(:),Ci_rs(:),rho_rs(:)
      REAL(kind=8) :: z_cnplus,rho_cnminus,rho_cnplus,mw_cnplus,A,B,C,D
      integer :: i,ncmps,c_min,c_max
!     Get pointers to vector data.
!     - For default PETSc vectors, VecGetArray() returns a pointer to
!     the data array.  Otherwise, the routine is implementation dependent.
!     - You MUST call VecRestoreArray() when you no longer need access to
!     the array.
!     - Note that the Fortran interface to VecGetArray() differs from the
!     C version.  See the Fortran chapter of the users manual for details.
    
      
      CALL VecGetArrayReadF90(Vec_ctx,pntr,ierr)
!   Allocate
      z_cnplus      = pntr(1)
      rho_cnminus   = pntr(2)
      c_min         = pntr(3)
      c_max         = pntr(4)
      ncmps         = c_max - c_min
      mw_cnplus     = pntr(5)
      rho_cnplus    = pntr(6)
      
      CALL VecRestoreArrayReadF90(Vec_ctx,pntr,ierr)
      ALLOCATE( z_rs(ncmps),MW_rs(ncmps),Ci_rs(ncmps),rho_rs(ncmps) )

! Read in last iteration results
      CALL VecGetArrayReadF90(Vec_x,pntr,ierr)    
        
      A = pntr(1)
      B = pntr(2)
      C = pntr(3)
      D = pntr(4)  
      
      CALL VecRestoreArrayReadF90(Vec_x,pntr,ierr)
!     Compute function

    ! Build arrays of correlations
      do i=1,ncmps 
        Ci_rs(i)  = (i-1) + c_min
        MW_rs(i)  = Ci_rs(i)*14 - 4
        z_rs(i)   = EXP(A + B*Ci_rs(i))
        rho_rs(i) = C + D*LOG(Ci_rs(i))
      end do
      
    ! Update function
      CALL VecGetArrayF90(Vec_f,pntr,ierr)
      pntr(1) = z_cnplus - SUM(z_rs)
      pntr(2) = z_cnplus*mw_cnplus - SUM( z_rs*MW_rs )
      pntr(3) = rho_cnplus - SUM(z_rs*MW_rs) / SUM( (z_rs*MW_rs)/rho_rs)
      pntr(4) = rho_cnminus - C + D*LOG(REAL(c_min - 1))
      CALL VecRestoreArrayF90(Vec_f,pntr,ierr)
!     Restore vectors

      RETURN
      END SUBROUTINE

! ---------------------------------------------------------------------
!
!  FormJacobian - Evaluates Jacobian matrix.
!
!  Input Parameters:
!  snes - the SNES context
!  x - input vector
!  dummy - optional user-defined context (not used here)
!
!  Output Parameters:
!  A - Jacobian matrix
!  B - optionally different preconditioning matrix
!  flag - flag indicating matrix structure
!
      SUBROUTINE FormJacobian(snes,Vec_x,Mat_J,Mat_B,Vec_ctx,ierr)
      IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscsnes.h>

      SNES         snes
      Vec          Vec_x,Vec_ctx
      Mat          Mat_J,Mat_B
      PetscScalar  JACO(4,4)
      PetscScalar, pointer :: pntr(:)
      PetscErrorCode ierr
      PetscInt idx(4),n
!      PetscFortranAddr ctx(*)
  ! ctx
  ! 1) z from after last TBP cut on             
  ! 2) rho from last TPB cut
  ! 3) Cmin
  ! 4) Cmax
  ! 5) MW from after last TBP cut on  
  ! 6) rho from after last TBP cut on
!  Declarations for use with local arrays

! Local Variables     
      real(kind=8), allocatable :: z_rs(:),MW_rs(:),Ci_rs(:),rho_rs(:)
      REAL(KIND=8) :: A,B,C,D
      integer :: i,ncmps,c_min,c_max
!  Get pointer to vector data

      n = 4
      CALL VecGetArrayReadF90(Vec_ctx,pntr,ierr)
      !   Allocate
      c_min = pntr(3)
      c_max = pntr(4)
      ncmps = c_max - c_min
      CALL VecRestoreArrayReadF90(Vec_ctx,pntr,ierr)
      
      CALL VecGetArrayReadF90(Vec_x,pntr,ierr)
      A = pntr(1)
      B = pntr(2)
      C = pntr(3)
      D = pntr(4)
      CALL VecRestoreArrayReadF90(Vec_x,pntr,ierr)
      
      ALLOCATE( z_rs(ncmps),MW_rs(ncmps),Ci_rs(ncmps),rho_rs(ncmps) )

!     Compute function

    ! Build arrays of correlations
      do i=1,ncmps 
        Ci_rs(i)  = (i-1) + c_min
        MW_rs(i)  = Ci_rs(i)*14 - 4
        z_rs(i)   = EXP(A + B*Ci_rs(i))
        rho_rs(i) = C + D*LOG(Ci_rs(i))
      end do
      
      
!  Compute Jacobian entries and insert into matrix.
!   - Since this is such a small problem, we set all entries for
!     the matrix at once.
!   - Note that MatSetValues() uses 0-based row and column numbers
!     in Fortran as well as in C (as set here in the array idx).

    ! Manually set matrix indicies
      idx(1) = 0
      idx(2) = 1
      idx(3) = 2
      idx(4) = 3

      ! Manually construct jacobian matrix...it is not complicated

      ! Begin with setting all nonspecified values to zero 
      JACO(:,:) = 0.0

      ! dr1/dA
      JACO(1,1) = (-1.0)*SUM(z_rs)
      JACO(1,2) = (-1.0)*SUM(Ci_rs*z_rs)
      JACO(2,1) = (-1.0)*SUM(z_rs*MW_rs)
      JACO(2,2) = (-1.0)*SUM(Ci_rs*z_rs*MW_rs)
      
      JACO(3,3) = (-1.0)*( (SUM(z_rs*MW_rs)*SUM( (z_rs*MW_rs) / (rho_rs**2) ) ) / &
                      ( SUM( (z_rs*MW_rs)/rho_rs )**2) )
                      
      JACO(3,4) = (-1.0)*( (SUM(z_rs*MW_rs)*SUM( (z_rs*MW_rs*LOG(Ci_rs)) /(rho_rs**2))) / &
                      ( SUM( (z_rs*MW_rs)/rho_rs )**2) )
      JACO(4,3) = -1.0
      JACO(4,4) = REAL(c_min-1)  
    
      CALL MatSetValues(Mat_J,n,idx,n,idx,A,INSERT_VALUES,ierr)

!  Assemble matrix

      CALL MatAssemblyBegin(Mat_J,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Mat_J,MAT_FINAL_ASSEMBLY,ierr)

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------------80
      SUBROUTINE MyLineSearch(linesearch, lctx, ierr)
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscsnes.h>

      SNES              snes
      INTEGER           lctx
      Vec               x, f,g, y, w
      PetscReal         ynorm,gnorm,xnorm
      PetscBool         flag
      PetscErrorCode    ierr

      PetscScalar       mone

      mone = -1.0
      CALL SNESLineSearchGetSNES(linesearch, snes, ierr)
      CALL SNESLineSearchGetVecs(linesearch, x, f, y, w, g, ierr)
      CALL VecNorm(y,NORM_2,ynorm,ierr)
      CALL VecAXPY(x,mone,y,ierr)
      CALL SNESComputeFunction(snes,x,f,ierr)
      CALL VecNorm(f,NORM_2,gnorm,ierr)
      CALL VecNorm(x,NORM_2,xnorm,ierr)
      CALL VecNorm(y,NORM_2,ynorm,ierr)
      CALL SNESLineSearchSetNorms(linesearch, xnorm, gnorm, ynorm,      &
     & ierr)
      flag = PETSC_FALSE
      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------------80      

