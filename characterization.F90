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
  KSP :: Krylov
  PC :: PreCon
  ! Main Variable Collection
  SNES     snes
  Vec      Vec_v,Vec_r,Vec_v0,Vec_dV,Vec_RHS,Vec_r0
  Vec      Vec_v_Jac,Vec_v0_Jac,Vec_r_Jac,Vec_dv_Jac,Vec_r0_Jac,Vec_TX
  Mat      Mat_J
  PetscErrorCode  ierr
  PetscInt n,divtol,Jn,J_AB_ind(2),J_CD_ind(2)
  PetscMPIInt size,rank
  PetscReal   tol
  IS :: J_AB, J_CD
  INTEGER :: i
  ! Program input variables
  LOGICAL :: l
  CHARACTER(80) :: input_file
  CHARACTER(80),ALLOCATABLE :: compname(:)
  INTEGER :: nc, nphase,c_min,c_max,ncmps
  REAL(KIND=8) :: P,T,feed,Pc_c7pi,Tc_c7pi,MW_c7pi,rho_c7pi,               &
       m_c7pi(1),omega_c7pi(1)

  ! EOS (SRK) Related Variables
  REAL(KIND=8), ALLOCATABLE :: Tcr(:),Pcr(:),omega(:),MW(:),rho(:),z(:),   &
       x(:),y(:),m(:),coeff(:),vol_shift(:)

  ! Characterization Variables
  INTEGER :: lastTBProw
  REAL(kind=8) :: z_cnplus,rho_cnminus,MW_cnplus,x_init(4),rho_cnplus,f1

  ! Solution Related Variables
  PetscInt  ::  init_ind(4)
  

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
  c_min = 11
  c_max = 200
  ncmps         = c_max - c_min
  ! PETSc Stuff
  n = 4
  tol = 1.0E-4
  divtol = 20
  ctx_size = 6
  Jn = 2
  
  J_AB_ind(0) = 0
  J_AB_ind(1) = 1  
  J_CD_ind(0) = 2
  J_CD_ind(1) = 3  
  ! Initial Values for Solver
  x_init(1) = -3.0  ! A
  x_init(2) = -0.1  ! B
  x_init(3) =  0.5  ! C
  x_init(4) =  0.1  ! D
  ! Initial Values for Solver / index
  init_ind(1) = 0
  init_ind(2) = 1
  init_ind(3) = 2
  init_ind(4) = 3      

  ! Other nonsense
  f1 = 1.0d0

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

  ! 1) z from after last TBP cut on
  z_cnplus = SUM(z(lastTBProw:))

  ! 2) rho from last TPB cut
  rho_cnminus = rho(lastTBProw)

 


  ! 5) MW from after last TBP cut on
  MW_cnplus = SUM(z(lastTBProw:)*MW(lastTBProw:)) / SUM(z(lastTBProw:))

  ! 6) rho from after last TBP cut on
  rho_cnplus = SUM(z(lastTBProw:)*MW(lastTBProw:)) /                               &
           SUM( (z(lastTBProw:)*MW(lastTBProw:))/rho(lastTBProw:) )

  !===========================================================================80
  ! Detemine Parameters A, B, C, D
  !---------------------------------------------------------------------------80
  ! Process begins after last TBP HC, from C11+ in this case.

  ! Create Solution/NL Function vectors
  CALL SNESCreate(PETSC_COMM_WORLD,snes,ierr)
  CALL VecCreateSeq(PETSC_COMM_SELF,n,Vec_v,ierr)
  CALL VecDuplicate(Vec_v,Vec_r,ierr)
  CALL VecDuplicate(Vec_v,Vec_r0,ierr)  
  CALL VecDuplicate(Vec_v,Vec_v0,ierr)


  ! Create Jacobian Matrix
  CALL MatCreate(PETSC_COMM_SELF,Mat_J,ierr)
  CALL MatSetSizes(Mat_J,PETSC_DECIDE,PETSC_DECIDE,Jn,Jn,ierr)
  CALL MatSetFromOptions(Mat_J,ierr)
  CALL MatSetUp(Mat_J,ierr)
  
  ! Create Index Set framework because I screwed up and tried to solve it 
  ! all up
  CALL ISCreateGeneral(PETSC_COMM_SELF,Jn,J_AB_ind,PETSC_COPY_VALUES,J_AB,ierr)
  CALL ISCreateGeneral(PETSC_COMM_SELF,Jn,J_CD_ind,PETSC_COPY_VALUES,J_CD,ierr)
  ! Create retarded solver vectors
  CALL VecCreateSeq(PETSC_COMM_SELF,Jn,Vec_v_Jac,ierr)  
  CALL VecDuplicate(Vec_v_Jac,Vec_RHS,ierr)  
!  CALL VecDuplicate(Vec_v_Jac,Vec_v0_Jac,ierr)    
!  CALL VecDuplicate(Vec_v_Jac,Vec_r_Jac,ierr)
!  CALL VecDuplicate(Vec_v_Jac,Vec_r0_Jac,ierr)
  CALL VecDuplicate(Vec_v_Jac,Vec_dV,ierr)  
  !---------------------------------------------------------------------------80

  !===========================================================================80

  !===========================================================================80
  ! Customize and run nonlinear "solver"
  !---------------------------------------------------------------------------80
!  Create linear solver context

  call KSPCreate(PETSC_COMM_WORLD,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(Krylov,Mat_J,Mat_J,Different_Nonzero_Pattern,ierr)
#else
     call KSPSetOperators(Krylov,Mat_J,Mat_J,ierr)
#endif  
!  Set runtime options (e.g., -ksp_type <type> -pc_type <type>)
  call KSPSetFromOptions(Krylov,ierr)  
  CALL SetupKSPSolver

  CALL PrintMsg("Initial Guess")

  ! Zero information placeholders
  CALL VecZeroEntries(Vec_r0,ierr)
  CALL VecZeroEntries(Vec_r,ierr)
  CALL VecZeroEntries(Vec_v,ierr)
  CALL VecZeroEntries(Vec_v0,ierr)
  CALL VecZeroEntries(Vec_dV,ierr)  
  CALL MatZeroEntries(Mat_J,ierr)  

  ! Feed in initial guesses
  call VecSetValues(Vec_v,n,init_ind,x_init,INSERT_VALUES,ierr)
  CALL VecAssemblyBegin(Vec_v,ierr)
  CALL VecAssemblyEnd(Vec_v,ierr)
      ! Find r vector for first iteration
  CALL PrintMsg("Solving")  

  do i=1,1
       ! Find r vector for iteration
      CALL FormFunction(Vec_v,Vec_r,ierr)  
      CALL VecGetSubvector(Vec_r,J_AB,Vec_r_Jac,ierr)
      CALL VecGetSubvector(Vec_r0,J_AB,Vec_r0_Jac,ierr)
      CALL VecGetSubvector(Vec_v,J_AB,Vec_v_Jac,ierr)      
      CALL VecGetSubvector(Vec_v0,J_AB,Vec_v0_Jac,ierr)

      ! Form Jacobian Matrix
      CALL FormJacobian(Jn,Vec_r_Jac,Vec_r0_Jac,Vec_v_Jac,Vec_v0_Jac,Mat_J,ierr)
      ! Move r vector to RHS
      CALL VecCopy(Vec_r_Jac,Vec_RHS,ierr)
      CALL VecView(Vec_RHS,PETSC_VIEWER_STDOUT_SELF,ierr)
      CALL VecScale(Vec_RHS,-f1,ierr)
      ! Solve linear system [J][V] = -[r]
      CALL KSPSolve(Krylov,Vec_RHS,Vec_dV,ierr)
      ! Update residual vector
      CALL VecAXPY(Vec_v_Jac,f1,Vec_dV,ierr)


      ! Restore Vectors
      CALL VecRestoreSubVector(Vec_r,J_AB,Vec_r_Jac,ierr)
      CALL VecRestoreSubVector(Vec_r0,J_AB,Vec_r0_Jac,ierr)      
      CALL VecRestoreSubVector(Vec_v,J_AB,Vec_v_Jac,ierr)
      CALL VecRestoreSubVector(Vec_v0,J_AB,Vec_v0_Jac,ierr)

      ! Debug      
      CALL VecView(Vec_v,PETSC_VIEWER_STDOUT_SELF,ierr)      
      CALL VecView(Vec_r,PETSC_VIEWER_STDOUT_SELF,ierr)
      CALL MatView(Mat_J,PETSC_VIEWER_STDOUT_SELF,ierr)          


      ! Backup prior iteration
      CALL VecCopy(Vec_v,Vec_v0,ierr)
       ! Backup prior iteration
      CALL VecCopy(Vec_r,Vec_r0,ierr)  
  end do
 
 
 
 
 

  CALL PrintMsg("Finished")
  call VecDestroy(Vec_v,ierr)
  call VecDestroy(Vec_v0,ierr)
  call VecDestroy(Vec_v_Jac,ierr)
  call VecDestroy(Vec_v0_Jac,ierr)
  call VecDestroy(Vec_r_Jac,ierr)
  call VecDestroy(Vec_r0_Jac,ierr)  
  call VecDestroy(Vec_r,ierr)
  call VecDestroy(Vec_RHS,ierr)  
  call MatDestroy(Mat_J,ierr)
  call SNESDestroy(snes,ierr)
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

!=============================================================================80      
  subroutine SetupKSPSolver
    implicit none
    call KSPSetType(Krylov,"gmres",ierr)
    call KSPGetPC(Krylov,PreCon,ierr)
    call PCSetType(PreCon,"asm",ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Double_Precision,        &
       Petsc_Default_Double_Precision,Petsc_Default_Integer,ierr)
#else
    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Real,Petsc_Default_Real, &
       Petsc_Default_Integer,ierr)
#endif
    call KSPSetFromOptions(Krylov,ierr)
  end subroutine SetupKSPSolver

!-----------------------------------------------------------------------------80      
END PROGRAM main
