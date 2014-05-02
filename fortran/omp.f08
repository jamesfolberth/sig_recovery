! omp.f08

module omp

   use compla ! also has BLAS interface for functions

   implicit none

   contains

   ! omp_algo
   subroutine omp_algo(s_hat,m,Phi,v,realloc)
   ! {{{
      real (kind=dblk) :: s_hat(:), Phi(:,:), v(:)
      integer (kind=intk) :: m
      logical :: realloc

      real (kind=dblk), allocatable, save :: Q(:,:),R(:,:),Phi_t(:,:),&
         resid(:),temp(:),rhs(:),x(:),a(:)
      integer (kind=intk), allocatable, save :: Lambda(:)
      integer (kind=intk) :: N,d,t

      N = size(Phi,1)
      d = size(Phi,2)

      if ( realloc ) then
         if ( allocated(Q) ) then
            deallocate(Q,R,Phi_t,resid,Lambda,temp,rhs,x,a)
         end if
         allocate(Q(N,N),R(N,m),Phi_t(N,m),resid(N),Lambda(m),temp(d),&
            rhs(N),x(m),a(N))
         realloc = .false.
      end if

      ! This didn't fix the problem
      !Q = 0.0_dblk
      !R = 0.0_dblk
      !resid = 0.0_dblk
      !Lambda = 0
      !temp = 0.0_dblk
      !rhs = 0.0_dblk
      !x = 0.0_dblk
      !a = 0.0_dblk
      !s_hat = 0.0_dblk

      ! 1 - initialization
      call dcopy(N,v,1,resid,1) ! resid <- v
      do t=1,m
         ! 2 - find index
         ! temp <- 1.0*Phi**T*resid + 0.0*temp
         call dgemv('T',N,d,1.0_dblk,Phi,N,resid,1,0.0_dblk,temp,1)
         Lambda(t) = idamax(d,temp,1)

         ! 3 - augment index set and matrix of atoms
         ! augment index set done in #2
         !Phi_t(:,t) = Phi(:,Lambda(t))
         call dcopy(N,Phi(1,Lambda(t)),1,Phi_t(1,t),1)

         ! 4 - solve least squares problem
         call qrappendcol(Q,R,t,Phi(:,Lambda(t)))
         call dgemv('T',N,N,1.0_dblk,Q,N,v,1,0.0_dblk,rhs,1)
         call dcopy(t,rhs,1,x,1)
         call back_solve_blk(R,x,t)

         ! 5 - update approximation and residual
         call dgemv('N',N,t,1.0_dblk,Phi_t,N,x,1,0.0_dblk,a,1)
         call dcopy(N,v,1,resid,1)
         call daxpy(N,-1.0_dblk,a,1,resid,1)

      end do
     
      ! 6 - construct approximate signal
      call dscal(d,0.0_dblk,s_hat,1)
      do t=1,m
         s_hat(Lambda(t)) = x(t)    
      end do

   end subroutine omp_algo
   ! }}}

   ! check_recovery
   function check_recovery(s,s_hat)
   ! {{{
      real (kind=dblk) :: s(:),s_hat(:)
      integer (kind=intk) :: check_recovery

      !integer (kind=intk) :: i
      !do i=1,size(s)
      !   print *, s(i),s_hat(i)
      !end do
      
      check_recovery = 0
      !print *,"norm = ", norm_p(s-s_hat,0)
      if ( norm_p(s-s_hat,0) < 10E-14_dblk ) then
         check_recovery = 1
      end if

   end function check_recovery
   ! }}}

   ! gen_sig
   ! generate sparse signal s with m non-zero entries
   ! if fillval_in is supplied, fill with fillval_in; otherwise uniform
   ! random from [-1,1]
   subroutine gen_sig(s,m,fillval_in)
   ! {{{
      real (kind=dblk) :: s(:)
      integer (kind=intk) :: m
      real (kind=dblk), optional :: fillval_in

      real (kind=dblk) :: temp
      integer (kind=dblk) :: i,sparse_inds(m),possible_ind

      if ( m > size(s) ) then
         stop('omp: gen_sig: cannot have an m-sparse signal with fewer than m elements')
      end if

      ! get the indexes we're going to fill
      i = 1
      sparse_inds = 0
      do while( i <= m) ! do until we fill sparse_inds
         ! sample from  ceiling([1/2,d-1/2]), which should give us
         ! a sampling from [1,2,3,...,d]  d = size(s)
         call random_number(temp)
         possible_ind = ceiling((dble(size(s))-1.0_dblk)*temp+0.5_dblk)

         ! possible_ind not in sparse_inds
         if ( .not. any(sparse_inds == possible_ind) ) then
            sparse_inds(i) = possible_ind
            i = i + 1
         end if
      end do

      ! zero out s
      call dscal(size(s),0.0_dblk,s,1)

      ! fill in the sparse vector 
      if ( present(fillval_in) ) then
         do i=1,size(sparse_inds)
            s(sparse_inds(i)) = fillval_in
         end do
      else
         do i=1,size(sparse_inds)
            call random_number(temp)
            s(sparse_inds(i)) = 2.0_dblk*temp-1.0_dblk ! uniform random on [-1,1]
         end do
      end if

   end subroutine gen_sig
   ! }}}

end module omp

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
