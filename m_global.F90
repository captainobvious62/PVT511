! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module global
    use roots


    IMPLICIT NONE
   
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscvec.h90>    
   
    ! GLOBAL Variables
  INTEGER :: nc, nphase,c_min,c_max,ncmps
  REAL(KIND=8) :: P,T,feed,Pc_c7pi,Tc_c7pi,MW_c7pi,rho_c7pi,               &
       m_c7pi(1),omega_c7pi(1),eptol

  ! EOS (SRK) Related Variables
  REAL(KIND=8), ALLOCATABLE :: Tcr(:),Pcr(:),omega(:),MW(:),rho(:),z(:),   &
       x(:),y(:),m(:),coeff(:),vol_shift(:)

  ! Characterization Variables
  INTEGER :: lastTBProw
  REAL(kind=8) :: z_cnplus,rho_cnminus,MW_cnplus,x_init(4),rho_cnplus,f1
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

  CHARACTER(80),ALLOCATABLE :: compname(:)

  ! Solution Related Variables
  PetscInt  ::  init_ind(4)
  
  ! Non PETSc Characterization Variables
  REAL(KIND=8)  :: residual(2),solution(2),jacobian(2,2)
  
  ! Lazy FD Stencil
  REAL(KIND=8)  :: FD_E(2),FD_W(2)
!=============================================================================80
CONTAINS  
!=============================================================================80
    SUBROUTINE FindR_AB(A,B,zcp,mcp,cmin,cmax,rvec)
        IMPLICIT NONE
        REAL(KIND=8)    :: A,B,zcp,mcp,rvec(2)
        REAL(KIND=8),ALLOCATABLE :: r1sum(:),r2sum(:)
        INTEGER         :: i,j,k,cmin,cmax,ncomps,scn
        
        ncomps = cmax-cmin
        ALLOCATE( r1sum(ncomps),r2sum(ncomps) )
        
        do scn=cmin,cmax
            i = scn - (cmin-1)
            r1sum(i) = EXP(A + B*REAL(scn)) 
            r2sum(i) = EXP(A + B*REAL(scn))*(14.*REAL(scn) - 4.)
        end do
        rvec(1) = zcp - SUM(r1sum)
        rvec(2) = zcp*mcp - SUM(r2sum)
        DEALLOCATE( r1sum,r2sum )
    END SUBROUTINE FindR_AB
!-----------------------------------------------------------------------------80            
                    


    SUBROUTINE FormFunction(Vec_v,Vec_r,ierr)
      IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscvec.h90>

      Vec      Vec_r,Vec_v
      PetscErrorCode ierr
      PetscScalar, pointer :: pntr(:)

! Local Variables     
      real(kind=8), allocatable :: z_rs(:),MW_rs(:),Ci_rs(:),rho_rs(:)
      REAL(kind=8) :: A,B,C,D
      integer :: i
!     Get pointers to vector data.
!     - For default PETSc vectors, VecGetArray() returns a pointer to
!     the data array.  Otherwise, the routine is implementation dependent.
!     - You MUST call VecRestoreArray() when you no longer need access to
!     the array.
!     - Note that the Fortran interface to VecGetArray() differs from the
!     C version.  See the Fortran chapter of the users manual for details.
      
      ALLOCATE( z_rs(ncmps),MW_rs(ncmps),Ci_rs(ncmps),rho_rs(ncmps) )

! Read in last iteration results
      CALL VecGetArrayReadF90(Vec_v,pntr,ierr)    
        
      A = pntr(1)
      B = pntr(2)
      C = pntr(3)
      D = pntr(4)  
      
      CALL VecRestoreArrayReadF90(Vec_v,pntr,ierr)
!     Compute function
    ! Build arrays of correlations
      do i=1,ncmps 
        Ci_rs(i)  = REAL(i-1) + c_min
        MW_rs(i)  = Ci_rs(i)*14. - 4.
        z_rs(i)   = EXP(A + B*Ci_rs(i))
        rho_rs(i) = C + D*LOG(Ci_rs(i))
      end do
      
      ! Update Particulars
      
      
    ! Update function
      CALL VecGetArrayF90(Vec_r,pntr,ierr)
      pntr(1) = z_cnplus - SUM(z_rs)
      pntr(2) = z_cnplus*MW_cnplus - SUM( z_rs*MW_rs )
      pntr(3) = rho_cnplus - SUM(z_rs*MW_rs) / SUM( (z_rs*MW_rs)/rho_rs)
      pntr(4) = rho_cnminus - C + D*LOG(REAL(c_min - 1))
      CALL VecRestoreArrayF90(Vec_r,pntr,ierr)
      !CALL VecView(Vec_r,PETSC_VIEWER_STDOUT_SELF,ierr)
!     Restore vectors

      END SUBROUTINE FormFunction
!-----------------------------------------------------------------------------80      

      SUBROUTINE AnalyticalJacobian(Vec_v,Mat_J,ierr)
      IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscsnes.h>

      Vec          Vec_v
      Mat          Mat_J
      PetscScalar  JACO(4,4)
      PetscScalar, pointer :: pntr(:)
      PetscErrorCode ierr
! Local Variables     
      real(kind=8), allocatable :: z_rs(:),MW_rs(:),Ci_rs(:),rho_rs(:)
      REAL(KIND=8) :: A,B,C,D
      integer :: i
!  Get pointer to vector data
      
      
      CALL VecGetArrayReadF90(Vec_v,pntr,ierr)
      A = pntr(1)
      B = pntr(2)
      C = pntr(3)
      D = pntr(4)
      CALL VecRestoreArrayReadF90(Vec_v,pntr,ierr)
      
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
      CALL MatSetValues(Mat_J,n,init_ind,n,init_ind,JACO,INSERT_VALUES,ierr)

!  Assemble matrix

      CALL MatAssemblyBegin(Mat_J,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Mat_J,MAT_FINAL_ASSEMBLY,ierr)
      END SUBROUTINE      
!-----------------------------------------------------------------------------80
!=============================================================================80
      SUBROUTINE FormJacobian(n,Vec_r,Vec_r0,Vec_v,Vec_v0,Mat_J,ierr)
      IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscvec.h90>
      Vec          Vec_v,Vec_r,Vec_r0,Vec_v0
      Mat          Mat_J
      PetscScalar, pointer :: pntr(:)
      PetscErrorCode ierr
      PetscInt id1(2),id2(2),n
      PetscInt, allocatable :: idx(:)
      PetscScalar, Allocatable :: JAC(:,:)
      real(kind=8), allocatable :: V(:),V0(:),R(:),R0(:)
      integer :: i,j
      
      ALLOCATE( V(n),V0(n),R(n),R0(n),JAC(n,n),idx(n) )
      
!  Get pointer to vector data
      CALL VecGetArrayReadF90(Vec_v,pntr,ierr)
      do i=1,n
        V(i) = pntr(i)
      end do
      CALL VecRestoreArrayReadF90(Vec_v,pntr,ierr)

      CALL VecGetArrayReadF90(Vec_v0,pntr,ierr)
      do i=1,n
        V0(i) = pntr(i)
      end do
      CALL VecRestoreArrayReadF90(Vec_v0,pntr,ierr)

      CALL VecGetArrayReadF90(Vec_r,pntr,ierr)
      do i=1,n
        R(i) = pntr(i)
      end do
      CALL VecRestoreArrayReadF90(Vec_r,pntr,ierr)

      CALL VecGetArrayReadF90(Vec_r0,pntr,ierr)
      do i=1,n
        R0(i) = pntr(i)
      end do
      CALL VecRestoreArrayReadF90(Vec_r0,pntr,ierr)


    ! Manually set matrix indicies
      do i=1,n
        idx(i) = i-1
      end do        
      ! Manually construct jacobian matrix...it is not complicated

      ! Begin with setting all nonspecified values to zero 
      CALL MatZeroEntries(Mat_J,ierr)
      JAC = 0.0
      
      do j=1,n
        do i=1,n
            JAC(i,j) = (R(i)-R0(i)) / (V(j)-V0(j))        
        end do
      end do                    
      print *, R
      print *, V
      print *, JAC
      
!  Assemble matrix
      CALL MatSetValues(Mat_J,n,idx,n,idx,JAC,INSERT_VALUES,ierr)      
      CALL MatAssemblyBegin(Mat_J,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(Mat_J,MAT_FINAL_ASSEMBLY,ierr)
      END SUBROUTINE FormJacobian
!-----------------------------------------------------------------------------80    
!=============================================================================80

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
!-----------------------------------------------------------------------------80
       
!-----------------------------------------------------------------------------80
end module global
