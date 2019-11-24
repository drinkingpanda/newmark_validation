include "mkl_pardiso.f90"
    
module ODE_solver
    implicit none
    
!------------------------- Module variables-----------------------------------

    contains
  
!*********************************************************************************************
! sub-001           solve a FEM transient response ODE equations in BSR format
!
!                            MU(..)+CU(.)+KU=F
!using standard Newmark Method:
!
!   U(.)_n+1 = U(.)_n                 + [ (1.0-beta)*U(..)_n +  beta*U(..)_n+1 ]*det_t
!   U_n+1    = U_n     + U(.)_n*det_t + [(0.5-alpha)*U(..)_n + alpha*U(..)_n+1 ]*det_t*det_t  
!
!**********************************************************************************************
subroutine newMark_bsr( N,lb,nnz,delt_t ,ia,ja,M,C,K,b,U_old,U1_old,U2_old,U_new,U1_new,U2_new)
implicit none

include 'mkl_pardiso.fi'

integer :: i,idum(1)
real*8  :: ddum(1)
!input parameters

integer  ,intent(in) :: N,lb,nnz
integer  ,intent(in) :: ia(N+1),ja(nnz)

real*8               :: alpha  ,beta     !alpha,beta for newmark and pardiso
real*8   ,intent(in) :: delt_t
       
real*8   ,intent(in) :: M(nnz),C(nnz),K(nnz)                        ! value array of matrix
real*8   ,intent(in) :: U_old(N*lb) , U1_old(N*lb) , U2_old(N*lb)   ! U1 ,U2 are first ,second time derivative of U
real*8   ,intent(in) ::     b(N*lb)                                 ! F in fem equation

!tmp parameters
real*8               ::     A(nnz )
real*8               ::    v1(N*lb) , v2(N*lb) ,v3(N*lb) , v4(N*lb) , rhs(N*lb)
real*8               ::    a1,a2,a3,a4,a5,a6,a7,a8,a9

!output parameters
real*8   ,intent(out):: U_new(N*lb) , U1_new(N*lb) , U2_new(N*lb)

!pardiso variables
TYPE(MKL_PARDISO_HANDLE) ::  pt(64)
integer                  ::  maxfct, mnum, mtype, phase, nrhs, error, msglvl
integer                  ::  iparm(64)

!---------------------------Newmark method preprocessor------------------------------------------- 
! default newmark parameters
alpha = 0.25d0
beta  = 0.5d0

! assemble the Matrix
    A  = K + ( beta/alpha/delt_t )*C + ( 1.0d0/alpha/delt_t/delt_t )*M   ! coefficient matrix A

! assemble the rhs vector
    v1 = (1.0d0/alpha/delt_t/delt_t)*U_old + (1.0d0/alpha/delt_t)*U1_old + (1.0d0/2.0d0/alpha - 1.0d0)*U2_old
    v2 = (      beta/alpha/delt_t  )*U_old + (beta/alpha-1.0d0  )*U1_old + (beta/alpha/2.0d0  - 1.0d0)*U2_old*delt_t

    call mkl_dbsrgemv("N", N, lb, M, ia, ja, v1, v3)  ! v3 = M * v1
    call mkl_dbsrgemv("N", N, lb, C, ia, ja, v2, v4)  ! v4 = C * v2

    rhs = v3 +v4 +b

!---------------------------Pardiso parameter explanations---------------------------------------- 
!
! pt        ¡ª¡ª  integer(64). handle to internal data
! maxfct    ¡ª¡ª  1.
! mnum      ¡ª¡ª  which matrix to factoriza, 1.
! mtype     ¡ª¡ª  £¨1£¬structurally symmetrix)£¬£¨11£¬real non-symmetric£©
! phase     ¡ª¡ª  £¨1£¬analysis and symbolic factorize£©,(2¡ª¡ªnumerical factorize)£¬£¨3¡ª¡ªsolve£©
! N         ¡ª¡ª  number of equations
! a,ia,ja   ¡ª¡ª  BSR3
! perm      ¡ª¡ª  integer(N).
! nrhs      ¡ª¡ª  number of rhs.
!
!iparm      ¡ª¡ª pass various parameters
!msglvl     ¡ª¡ª whether output message
!b          ¡ª¡ª real(nrhs,N),rhs vectors. the array is replaced with solution if iparm(6)=1.
!x          ¡ª¡ª solution vector,if iparm(6)=0.
!error      ¡ª¡ª error indicator


!pardiso parameters

  DO i = 1, 64
     iparm(i)      = 0
     pt(i).DUMMY   = 0
  END DO

 iparm(1)  =  1 ! no solver default
 iparm(2)  =  2 ! fill-in reordering from METIS
 iparm(3)  =  2 ! numbers of processors
 iparm(4)  =  0 ! no iterative-direct algorithm
 iparm(8)  =  9 ! numbers of iterative refinement steps
 iparm(10) =  13 ! perturb the pivot elements with 1E-13
 iparm(11) =  1 ! use nonsymmetric permutation and scaling MPS
 iparm(13) =  1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
 iparm(18) = -1 ! Output: number of nonzeros in the factor LU
 iparm(19) = -1 ! Output: Mflops for LU factorization

 
 iparm(37) =  1 ! if iparm(37)=0,then CSR is used. if iparm(37)>0,then BSR is used and the number is block_size(6 of 6*6 block)
 
 error  = 0   ! initialize error flag
 msglvl = 0   ! print statistical information
 mtype  = 11  ! real unsymmetric
 nrhs   = 1   ! number of rhs
 maxfct = 1   ! max factorization
 mnum   = 1   ! the matrix to factorize


    ! phase = 13
    !
    !call pardiso (pt, maxfct, mnum, mtype, phase, N*lb, A, ia, ja, perm, nrhs, iparm, msglvl , rhs , U_new, error)
    !
    !
    
    !.. Reordering and Symbolic Factorization, This step also allocates
    ! all memory that is necessary for the factorization

    phase = 11 ! only reordering and symbolic factorization

    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
                  idum, nrhs, iparm, msglvl, ddum, ddum, error)
    
    !.. Factorization.
    phase = 22 ! only factorization
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
                  idum, nrhs, iparm, msglvl, ddum, ddum, error)
    
    !.. Back substitution and iterative refinement
    iparm(8) = 2 ! max numbers of iterative refinement steps
    phase = 33 ! only solving

    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
                  idum, nrhs, iparm, msglvl, rhs , U_new, error)
    
!----------------Post Pardiso part----------------------------------- 
    
    !compute U2_new
    
    a1 =  1.0d0/alpha/delt_t/delt_t
    a2 =  1.0d0/alpha/delt_t
    a3 =  1.0d0/2.0d0/alpha - 1.0d0
    
    U2_new = a1*U_new - ( a1*U_old + a2*U1_old + a3*U2_old )
    
    ! compute U1_new 
    
    a4 = beta/alpha/delt_t
    a5 = beta/alpha-1.0d0
    a6 = (beta/alpha/2.0d0-1.0d0)*delt_t
    
    U1_new = a4*U_new - ( a4*U_old + a5*U1_old +a6*U2_old    )
    
    !!.. Termination and release of memory
    !phase = -1 ! release internal memory
    !CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
    !           idum, nrhs, iparm, msglvl, rhs , U_new, error) 
    !
    
end subroutine newMark_BSR




end module ODE_solver
