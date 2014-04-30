! compla.f08
! Computational Linear Algebra library
! James Folberth - Spring 2014

! TODO make `block_size` global (separate block sizes per algorithm)
! TODO make BLAS (ATLAS) calls where possible
! TODO overload print_array to print vectors as well
! TODO make vector 2-norm-squared function for QR decomp
! TODO take care of b in QR decomposition (loop and blas)

! wrap lib in module so it interfaces properly
module compla

   use HDF5 ! save stuff

   implicit none

   integer, parameter :: intk = kind(1)
   integer, parameter :: dblk = kind(1d0) ! TODO 4->intk, 8->dblk

   real (kind=8), parameter :: ZERO_TOL = 10**(-14)

   ! p-norm function overload
   public norm_p
   private norm_p_mat, norm_p_vec

   interface norm_p
      module procedure norm_p_mat, norm_p_vec
   end interface

   public write_dset
   private dwrite_dset_rank1, iwrite_dset_rank1,&
           dwrite_dset_rank2, dwrite_dset_rank0

   interface write_dset
      module procedure dwrite_dset_rank1, iwrite_dset_rank1,&
                       dwrite_dset_rank2, dwrite_dset_rank0
   end interface write_dset

   ! BLAS
   integer (kind=intk), external :: idamax


   contains

   !!!!!!!!!!!!!!!!!!!
   ! Cholesky Decomp !
   !!!!!!!!!!!!!!!!!!!
   ! {{{ 
   ! This uses the outer product formulation of the Cholesky decomposition
   ! Input matrix must by symmetric, positive definite.
   ! TODO optional input for transpose output or not?

   ! Outer product formulation
   subroutine chol(A)
   real (kind=8) :: A(:,:)

   integer (kind=4) :: i,j,k,Nc

   ! check that A is square
   if (size(A,1) /= size(A,2)) then
      print *, "error: compla.f08: chol: input matrix is not square"
      stop
   end if

   ! call check_sym(A)

   ! parameter to store size(A,1)
   Nc = size(A,1)

   ! main Cholesky loop
   ! Stores Cholesky factor in lower triangle of A, then transposes (column oriented)
   ! In my implementation, this is faster than the row oriented version (because Fortran stores column-wise)
   recurse: do i=1,Nc-1
      ! check a_{i,i} > 0
      if (A(i,i) <= ZERO_TOL) then
         print *, "error: compla.f08: chol: input matrix is not positive definite"
         stop
      end if

      A(i,i) = sqrt(A(i,i))
      A(i+1:Nc,i) = A(i+1:Nc,i) / A(i,i)
     
      ! Outer product A~ <- A^ - s*s'
      ! Only compute lower triangle
      row: do j=i+1,Nc
         col: do k=j,Nc
            A(k,j) = A(k,j)-A(j,i)*A(k,i)
         end do col
      end do row

   end do recurse

   if (A(Nc,Nc) <= ZERO_TOL) then
      print *, "error: compla.f08: chol: input matrix is not positive definite"
      stop
   end if
  
   A(Nc,Nc) = sqrt(A(Nc,Nc))

   ! transpose A, then zero out strictly lower triangle
   A = transpose(A)
   do j=1,Nc-1
      do i=j+1,Nc
        A(i,j) = 0d0
      end do
   end do
         
   end subroutine chol

   ! Outer product formulation with BLAS calls
   subroutine chol_blas(A)
   real (kind=8) :: A(:,:)

   integer (kind=4) :: i,j,Nc

   ! check that A is square
   if (size(A,1) /= size(A,2)) then
      print *, "error: compla.f08: chol_blas: input matrix is not square"
      stop
   end if

   ! call check_sym(A)

   ! parameter to store size(A,1)
   Nc = size(A,1)

   ! main Cholesky loop
   ! Stores Cholesky factor in lower triangle of A, then transposes
   recurse: do i=1,Nc-1
      ! check a_{i,i} > 0
      if (A(i,i) <= ZERO_TOL) then
         print *, "error: compla.f08: chol_blas: input matrix is not positive definite"
         stop
      end if

      A(i,i) = sqrt(A(i,i))
      ! upper triangular
      call dscal(Nc-i, 1d0/A(i,i), A(i,i+1), Nc) ! store in row
      ! lower triangular
      !call dscal(Nc-i, 1d0/A(i,i), A(i+1,i), 1) ! store in column
      ! I guess fortran passes the whole array or something

      ! Outer product A~ <- A^ - s*s**T
      ! upper triangular
      call dsyr('Upper', Nc-i, -1d0, A(i,i+1), Nc, A(i+1,i+1), Nc)
      ! lower triangular
      !call dsyr('Lower', Nc-i, -1d0, A(i+1,i), 1, A(i+1,i+1), Nc) 
      ! Note: don't pass the submatrix A(i+1:Nc,i+1:Nc), since that is WAY slower than giving dsyr its staring point, dimension, and
      ! leading dimension

   end do recurse

   if (A(Nc,Nc) <= ZERO_TOL) then
      print *, "error: compla.f08: chol_blas: input matrix is not positive definite"
      stop
   end if
  
   A(Nc,Nc) = sqrt(A(Nc,Nc))

   ! transpose A, then zero out strict lower triangle
   ! lower triangular only
   !A = transpose(A)
   do j=1,Nc-1
      !do i=j+1,Nc
      !  A(i,j) = 0d0
      !end do
      call dscal(Nc-j, 0d0, A(j+1,j), 1)
   end do
         
   end subroutine chol_blas


   ! Outer product formulation
   ! row oriented
   ! XXX don't use this; chol() is faster (at least for Nc>1000)
   subroutine chol_row(A)
   real (kind=8) :: A(:,:)

   integer (kind=4) :: i,j,k,Nc

   ! check that A is square
   if (size(A,1) /= size(A,2)) then
      print *, "error: compla.f08: chol_row: input matrix is not square"
      stop
   end if

   ! parameter to store size(A,1)
   Nc = size(A,1)

   ! main Cholesky loop
   ! Stores Cholesky factor in upper triangle of A; row oriented algorithm, even though FORTRAN stores column-wise
   recurse: do i=1,Nc-1
      ! check a_{i,i} > 0
      if (A(i,i) <= ZERO_TOL) then
         print *, "error: compla.f08: chol_row: input matrix is not positive definite"
         stop
      end if

      A(i,i) = sqrt(A(i,i))
      A(i,i+1:Nc) = A(i,i+1:Nc) / A(i,i)
     
      ! Outer product A~ <- A^ - s*s'
      ! Only compute upper triangle
      row: do j=i+1,Nc
         col: do k=i+1,j
            A(k,j) = A(k,j)-A(i,j)*A(i,k)
         end do col
      end do row
      !call print_array(A)

   end do recurse

   if (A(Nc,Nc) <= ZERO_TOL) then
      print *, "error: compla.f08: chol_row: input matrix is not positive definite"
      stop
   end if
  
   A(Nc,Nc) = sqrt(A(Nc,Nc))

   do j=1,Nc-1
      do i=j+1,Nc
         A(i,j) = 0d0
      end do
   end do
         
   end subroutine chol_row


   !subroutine chol_blk(A)
   ! TODO

   !subroutine check_sym()

   ! }}}


   !!!!!!!!!!!!!
   ! LU Decomp !
   !!!!!!!!!!!!!
   ! {{{
   subroutine lu(A,p)
      real (kind=8) :: A(:,:)
      integer (kind=4) :: p(:)
      
      integer (kind=4) :: i,j,k,Nc,Nr,m,temp_ind
      integer (kind=4) :: singular

      Nc = size(A,1)
      Nr = size(A,2)

      if (Nc /= Nr ) then
         print *, "error: compla.f08: lu: input matrix is not square"
         stop
      end if

      singular = 0
      
      ! this overwrites existing p
      p = (/ (k,k=1,Nc) /)

      row: do k=1,Nc-1

         m = col_find_absmax(A,p,k)
         !If pivot is near zero, matrix is singular
         if (abs(A(p(m),k)) < ZERO_TOL) then
            print *, "UNTESTED"
            print *, "warning: compla.f08: lu: input matrix is singular to within ZERO_TOL; U will be singular"
            singular = 1
            stop ! TODO move row to bottom of array and try to keep going

         else
            ! swap rows
            temp_ind = p(k)
            p(k) = p(m)
            p(m) = temp_ind

            multipliers: do i=k+1,Nc
               A(p(i),k) = A(p(i),k) / A(p(k),k)
               !print *,A(p(i),k)
            end do multipliers

            ! column oriented row operation
            row_op: do j=k+1,Nc
               do i=k+1,Nc
                  A(p(i),j) = A(p(i),j) - A(p(i),k)*A(p(k),j)
               end do
            end do row_op

         end if

         ! check A(end,end)
         if (abs(A(p(Nc),Nc)) < ZERO_TOL) then
            print *, "warning: compla.f08: lu: input matrix is singular to within ZERO_TOL; U will be singular"
            singular = 1
            stop
         end if
      end do row
         
   end subroutine lu


   ! LU without partial pivoting
   ! {{{
   subroutine lu_nopp(A)
      real (kind=8) :: A(:,:)
      integer (kind=4), allocatable :: p(:)      

      integer (kind=4) :: i,j,k,Nc,Nr,m
      integer (kind=4) :: singular

      Nc = size(A,1)
      Nr = size(A,2)

      if (Nc /= Nr ) then
         print *, "error: compla.f08: lu: input matrix is not square"
         stop
      end if

      singular = 0
      
      allocate(p(Nc))
      p = (/ (k,k=1,Nc) /)

      ! Note that I keep accessing (no swaps) the p 
      ! vector instead of removing it; I don't care
      ! about efficiency, since I'm never going to 
      ! use this in a real code
      row: do k=1,Nc-1

         m = col_find_absmax(A,p,k)
         !If pivot is near zero, matrix is singular
         if (abs(A(p(m),k)) < ZERO_TOL) then
            print *, "UNTESTED"
            print *, "warning: compla.f08: lu: input matrix is singular to within ZERO_TOL; U will be singular"
            singular = 1
            stop ! TODO move row to bottom of array and try to keep going

         else
            ! swap rows
            !temp_ind = p(k)
            !p(k) = p(m)
            !p(m) = temp_ind

            multipliers: do i=k+1,Nc
               A(p(i),k) = A(p(i),k) / A(p(k),k)
               !print *,A(p(i),k)
            end do multipliers

            ! column oriented row operation
            row_op: do j=k+1,Nc
               do i=k+1,Nc
                  A(p(i),j) = A(p(i),j) - A(p(i),k)*A(p(k),j)
               end do
            end do row_op

         end if

         ! check A(end,end)
         if (abs(A(p(Nc),Nc)) < ZERO_TOL) then
            print *, "warning: compla.f08: lu: input matrix is singular to within ZERO_TOL; U will be singular"
            singular = 1
            stop
         end if
      end do row

   end subroutine lu_nopp
   ! }}}


   function col_find_absmax(A,p,k)
      real (kind=8), intent(in) :: A(:,:)
      integer (kind=4), intent(in) :: p(:),k

      real (kind=8) :: aval,amax
      integer (kind=4) :: i,Nc,col_find_absmax

      Nc = size(A,1)

      amax = 0
      col_find_absmax = p(k)
      col: do i=k,Nc
         aval = abs(A(p(i),k))
         !print *,"find: ",p(i), aval
         if (aval > amax) then
            col_find_absmax=i
            amax = aval
         end if
      end do col
   end function col_find_absmax


   subroutine apply_perm_vector(A,p,trans)
      real (kind=8) :: A(:,:)
      integer (kind=4) :: p(:), trans

      ! ``compute'' P*
      if (trans == 0 ) then
         A(:,:) = A(p,:)
 
      ! `` compute P'*A
      else
         A(p,:) = A(:,:) 

      end if
   
   end subroutine apply_perm_vector


   subroutine form_LU(A,L,U)
      real (kind=8) :: A(:,:), L(:,:), U(:,:)

      integer (kind=4) :: i,j,Nc

      ! assume square
      Nc = size(A,1)

      L = 0d0
      U = 0d0

      row: do j=1,Nc
         L(j,j) = 1d0
         do i=j+1,Nc
            L(i,j) = A(i,j)
         end do 

         do i=j,1,-1
            U(i,j) = A(i,j)
         end do
      end do row

   end subroutine form_LU
   ! }}}


   !!!!!!!!!!!!!
   ! QR Decomp !
   !!!!!!!!!!!!!
   ! {{{
   ! QR decomp by reflectors
   ! Store u_k (to make Q) and R over A
   ! assume first entry of each u_k is identically 1
   subroutine qr(A,bin)
      real (kind=8) :: A(:,:)
      real (kind=8), optional :: bin(:)
     
      integer (kind=4) :: Nr,Nc,i,j,k
      real (kind=8) :: gamm, beta, tau
      real (kind=8), allocatable :: b(:),u(:),temp(:)

      Nr = size(A,1)
      Nc = size(A,2)

      allocate(u(Nr),temp(Nr),b(Nr))
     
      ! bin is the RHS of a linear system Ax=b
      if (present(bin)) then
         b = bin
      else 
         b = 0_dblk
      end if

      do k=1,Nc
         beta = maxval(abs(A(k:Nr,k)))
         if (beta <= ZERO_TOL) then
            ! gamma = 0, Q_k = I, no multiplication required
         else
            ! calculate the reflector (Watkins equation 3.2.35)
            u(k:Nr) = A(k:Nr,k) / beta
            tau = norm_p_vec(u(k:Nr),2)
            if (u(k)<0) tau = -tau
            u(k) = u(k) + tau
            gamm = u(k) / tau
            u(k+1:Nr) = u(k+1:Nr) / u(k)
            u(k) = 1d0
            tau = tau*beta

            if (abs(gamm) < ZERO_TOL) then
               print *, "error: compla.f08: qr: input matrix is singular to within ZERO_TOL"
               stop
            end if
           
            ! Compute Q_k*A_k:Nc,k+1:Nr
            temp = 0
            do j=k+1,Nc
               do i=k,Nr
                  temp(j) = temp(j)+u(i)*A(i,j)
               end do
            end do
            temp = gamm*temp

            do j=k+1,Nc
               do i=k,Nr
                  A(i,j) = A(i,j)-u(i)*temp(j)
               end do
            end do

            ! Compute Q_k*b_k:Nc
            temp(1) = 0
            do i=k,Nr
               temp(1) = temp(1) + u(i)*b(i)
            end do
            b(k:Nr) = b(k:Nr) - gamm*temp(1)*u(k:Nr)
               
            ! Store reflector vectors over A (u_k=1 is assumed)
            A(k,k) = -tau
            A(k+1:Nr,k) = u(k+1:Nr)
         end if
      end do

      if (present(bin)) bin = b

   end subroutine qr

   ! QR decomp by reflectors (BLAS calls)
   ! Store u_k's (to make Q) and R over A
   ! assume first entry of each u_k is identically 1
   subroutine qr_blas(A,bin)
      real (kind=8) :: A(:,:)
      real (kind=8), optional :: bin(:)
     
      integer (kind=4) :: Nr,Nc,i,k
      real (kind=8) :: gamm, beta, tau
      real (kind=8), allocatable :: b(:),u(:),v(:),temp(:)

      Nr = size(A,1)
      Nc = size(A,2)

      allocate(u(Nr),v(Nr),temp(Nr),b(Nr))
     
      ! bin is the RHS of a linear system Ax=b
      if (present(bin)) then
         b = bin
      else 
         b = 0_dblk
      end if

      do k=1,Nc
         !beta = maxval(abs(A(k:Nr,k)))
         beta = abs(A(k-1+idamax(Nr-k+1,A(k,k),1),k))
         if (beta <= ZERO_TOL) then
            ! gamma = 0, Q_k = I, no multiplication required
         else
            ! calculate the reflector (Watkins equation 3.2.35)
            !u(k:Nc) = A(k:Nc,k) / beta
            call dcopy(Nr-k+1,A(k,k),1,u(k),1)
            call dscal(Nr-k+1,1d0/beta,u(k),1)
            tau = norm_p_vec(u(k:Nr),2) 
            if (u(k)<0) tau = -tau
            u(k) = u(k) + tau
            gamm = u(k) / tau
            !u(k+1:Nc) = u(k+1:Nc) / u(k)
            call dscal(Nr-k,1d0/u(k),u(k+1),1)
            u(k) = 1d0
            tau = tau*beta

            if (abs(gamm) < ZERO_TOL) then
               print *, "error: compla.f08: qr_blas: input matrix is singular to within ZERO_TOL"
               stop
            end if
           
            ! Compute Q_k*A_k:Nc,k+1:Nr
            ! u**T*A = (A**T*u)**T
            !print *, Nr-k+1,Nc-k
            if (k < Nc) then
               call dgemv('T',Nr-k+1,Nc-k,1d0,A(k,k+1),Nr,u(k),1,0d0,temp(k+1),1)
               ! finish the rank-one update
               call dger(Nr-k+1,Nc-k, -gamm,u(k),1,temp(k+1),1,A(k,k+1),Nr)

               ! Compute Q_k*b_k:Nc
               temp(1) = 0
               do i=k,Nr
                  temp(1) = temp(1) + u(i)*b(i)
               end do
               b(k:Nr) = b(k:Nr) - gamm*temp(1)*u(k:Nr)
            end if
               
            ! Store reflector vectors over A (u_k=1 is assumed)
            A(k,k) = -tau
            !A(k+1:Nr,k) = u(k+1:Nr)
            call dcopy(Nr-k,u(k+1),1,A(k+1,k),1)
         end if
      end do

      if (present(bin)) bin = b

   end subroutine qr_blas


   ! Form Q,R from QR decomposition (this is expensive to call)
   subroutine form_qr(A,Q,R)
      real (kind=8) :: A(:,:), Q(:,:), R(:,:)

      integer (kind=4) :: Nr,Nc, j,k
      real (kind=8), allocatable :: u(:),temp(:)
      real (kind=8) :: gamm

      Nr = size(A,1)
      Nc = size(A,2)
      Q = eye(Nr,Nr)

      R = 0
      allocate(u(Nr),temp(Nr))

      ! R is just the upper triangle of A overwritten with \{u_k\} and R
      do j=1,Nc
         !do i=1,j
         !   R(i,j) = A(i,j)
         !end do
         call dcopy(j,A(1,j),1,R(1,j),1)
      end do
     
      ! Apply Q_k's (which form Q^T, k=Nc-1,..,1) to I to find Q*I
      ! Watkins exercise 3.2.44
      do k=Nc,1,-1
         temp = 0
         u(k) = 1d0
         !u(k+1:Nc) = A(k+1:Nc,k)
         call dcopy(Nr-k,A(k+1,k),1,u(k+1),1)
         gamm = 2/norm_p_vec(u(k:Nr),2)**2

         !do i=k,Nc
         !   do j=k,Nr
         !      temp(j) = temp(j)+u(i)*Q(i,j)
         !   end do
         !end do
         !temp = gamm*temp

         !do j=k,Nr
         !   do i=k,Nc
         !      Q(i,j) = Q(i,j)-u(i)*temp(j)
         !   end do
         !end do

         ! u**T*Q = (Q**T*u)**T
         call dgemv('T',Nr-k+1,Nr-k+1,1d0,Q(k,k),Nr,u(k),1,0d0,temp(k),1)
         ! finish the rank-one update
         call dger(Nr-k+1,Nr-k+1, -gamm,u(k),1,temp(k),1,Q(k,k),Nr)

      end do

   end subroutine form_qr

   ! }}}


   !!!!!!!!!!!!!!!!!!
   ! For/Back Solve !
   !!!!!!!!!!!!!!!!!!
   ! {{{
   subroutine back_solve(U,b)
      real (kind=8) :: U(:,:), b(:,:)
      integer (kind=4) :: i,j,Nc

      ! call check_square(U)
      ! call check_upper_tri(U)

      Nc=size(U,1)

      row: do j=Nc,1,-1
         ! zero on diagonal
         if (abs(U(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f08: back_solve: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/U(j,j)
         col: do i=j-1,1,-1
            b(i,1)=b(i,1) - U(i,j)*b(j,1)
         end do col
      end do row

   end subroutine back_solve

   subroutine back_solve_blk(U,b)
      real (kind=8) :: U(:,:), b(:,:)

      integer (kind=4), parameter :: blk_size=8
      integer (kind=4) :: i,j,s,Nc
      integer (kind=4) :: il, ih, jl, jh

      ! call check_square(U)
      ! call check_upper_tri(U)

      Nc = size(U,1)
      s = Nc / blk_size

      ! Column oriented backward substitution
      blk: do j=1,s
         if (j>1) then
            col_blk: do i=j-1,s-1
               il = Nc-blk_size*(i+1)+1
               ih = Nc-blk_size*i
               jl = Nc-blk_size*(j-1)+1
               jh = Nc-blk_size*(j-2)
               !print *, "subtract block: ", il,ih,jl,jh
               
               b(il:ih,1) = b(il:ih,1) - matmul(U(il:ih,jl:jh),b(jl:jh,1))

            end do col_blk

            ! top block (not necc. blk_size by blk_size)
            il = 1
            ih = Nc-blk_size*s
            jl = Nc-blk_size*(j-1)+1
            jh = Nc-blk_size*(j-2)
            !print *, "subtract top block: ", il,ih,jl,jh

            b(il:ih,1) = b(il:ih,1) - matmul(U(il:ih,jl:jh),b(jl:jh,1))

         end if

         ! call back_solve on the diagonal blocks
         jl = Nc-blk_size*j+1
         jh = Nc-blk_size*(j-1)
         !print *, "back_solve: ",jl,jh

         call back_solve(U(jl:jh,jl:jh),b(jl:jh,:))

      end do blk

      ! subtract final top block
      if (s>0) then
         il = 1
         ih = Nc-blk_size*s
         jl = Nc-blk_size*s+1
         jh = Nc-blk_size*(s-1)
         !print *, "subtract final top block: ", il,ih,jl,jh
         
         b(il:ih,1) = b(il:ih,1) - matmul(U(il:ih,jl:jh),b(jl:jh,1))
      end if

      ! Finish with regular back solve
      row_fin: do j=Nc-blk_size*s,1,-1
         ! zero on diagonal
         if (abs(U(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f08: back_solve: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/U(j,j)
         col_fin: do i=j-1,1,-1
            b(i,1)=b(i,1) - U(i,j)*b(j,1)
         end do col_fin
      end do row_fin

   end subroutine back_solve_blk

   
   subroutine for_solve(L,b)
      real (kind=8) :: L(:,:), b(:,:)
      integer (kind=4) :: i,j,Nc

      ! call check_square(L)
      ! call check_lower_tri(L)

      Nc = size(L,1)

      row: do j=1,Nc
         ! zero on diagonal means singular matrix
         if (abs(L(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f08: for_solve: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/L(j,j)
         col: do i=j+1,Nc
            b(i,1)=b(i,1) - L(i,j)*b(j,1) 
         end do col
      end do row

   end subroutine for_solve

   ! This is the same as the for_solve routine, but assumes
   ! that the main diagonal of L is filled with ones 
   ! expected input is A, overwritten with L,U from LU decomp.
   subroutine for_solve_lu(L,b)
      real (kind=8) :: L(:,:), b(:,:)
      integer (kind=4) :: i,j,Nc

      ! call check_square(L)
      ! call check_lower_tri(L)

      Nc = size(L,1)

      row: do j=1,Nc
         ! it is assumed that L(j,j)=1

         ! zero on diagonal means singular matrix
         !if (abs(L(j,j)) <= ZERO_TOL) then
         !   print *, "error: compla.f08: for_solve: input matrix is singular to tolerance"
         !   stop
         !end if

         !b(j,1) = b(j,1)/L(j,j) 
         col: do i=j+1,Nc
            b(i,1)=b(i,1) - L(i,j)*b(j,1) 
         end do col
      end do row

   end subroutine for_solve_lu

   subroutine for_solve_blk(L,b)
      real (kind=8) :: L(:,:), b(:,:)

      integer (kind=4), parameter :: blk_size=8
      integer (kind=4) :: i,j,s,Nc
      integer (kind=4) :: il,ih,jl,jh

      ! call check_square(L)
      ! call check_lower_tri(L)

      Nc = size(L,1)
      s = Nc / blk_size

      ! Do forward subs in square blocks as far as possible
      ! column oriented
      blk: do j=1,s
         !print *, "j=",j
         if (j>1) then
            col_blk: do i=j,s
               il = blk_size*(i-1)+1
               ih = blk_size*i
               jl = blk_size*(j-2)+1
               jh = blk_size*(j-1)
               ! print *, il,ih,jl,jh
               
               b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))

            end do col_blk

            ! subtract bottom block (not necc. of size blk_size by blk_size)
            il = blk_size*s+1
            ih = Nc
            jl = blk_size*(j-2)+1
            jh = blk_size*(j-1)
            ! print *, il,ih,jl,ih

            b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))

         end if
 
         !print *, "calling for_solve"
         jl = blk_size*(j-1)+1
         jh = blk_size*j

         call for_solve(L(jl:jh,jl:jh), b(jl:jh,:))

      end do blk

      ! subtract final bottom block
      if (s>0) then
         il = blk_size*s+1
         ih = Nc
         jl = blk_size*(s-1)+1
         jh = blk_size*s

         b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))
      end if

      ! Finish up with regular forward subs
      !print *, blk_size*s
      row_fin: do j=blk_size*s+1,Nc
         ! zero on diagonal means singular matrix
         if (abs(L(j,j)) <= ZERO_TOL) then
            print *, "error: compla.f08: for_solve_blk: input matrix is singular to tolerance"
            stop
         end if

         b(j,1) = b(j,1)/L(j,j)
         col_fin: do i=j+1,Nc
            b(i,1)=b(i,1) - L(i,j)*b(j,1) 
         end do col_fin
      end do row_fin

   end subroutine for_solve_blk

   ! This is the same as for_solve_blk, except that is assumes ones along the main diagonal of L
   ! Use this routine for L,U overwritten on A
   subroutine for_solve_lu_blk(L,b)
      real (kind=8) :: L(:,:), b(:,:)

      integer (kind=4), parameter :: blk_size=8
      integer (kind=4) :: i,j,s,Nc
      integer (kind=4) :: il,ih,jl,jh

      ! call check_square(L)
      ! call check_lower_tri(L)

      Nc = size(L,1)
      s = Nc / blk_size

      ! Do forward subs in square blocks as far as possible
      ! column oriented
      blk: do j=1,s
         !print *, "j=",j
         if (j>1) then
            col_blk: do i=j,s
               il = blk_size*(i-1)+1
               ih = blk_size*i
               jl = blk_size*(j-2)+1
               jh = blk_size*(j-1)
               ! print *, il,ih,jl,jh
               
               b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))

            end do col_blk

            ! subtract bottom block (not necc. of size blk_size by blk_size)
            il = blk_size*s+1
            ih = Nc
            jl = blk_size*(j-2)+1
            jh = blk_size*(j-1)
            ! print *, il,ih,jl,ih

            b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))

         end if
 
         !print *, "calling for_solve"
         jl = blk_size*(j-1)+1
         jh = blk_size*j

         ! XXX this should be the only difference between
         !     for_solve_blk and for_solve_lu_blk
         call for_solve_lu(L(jl:jh,jl:jh), b(jl:jh,:))

      end do blk

      ! subtract final bottom block
      if (s>0) then
         il = blk_size*s+1
         ih = Nc
         jl = blk_size*(s-1)+1
         jh = blk_size*s

         b(il:ih,1) = b(il:ih,1) - matmul(L(il:ih,jl:jh),b(jl:jh,1))
      end if

      ! Finish up with regular forward subs
      !print *, blk_size*s
      row_fin: do j=blk_size*s+1,Nc
         ! zero on diagonal means singular matrix
         ! L(j,j) = 1, by assumption
         !if (abs(L(j,j)) <= ZERO_TOL) then
         !   print *, "error: compla.f08: for_solve_blk: input matrix is singular to tolerance"
         !   stop
         !end if

         b(j,1) = b(j,1) ! L(j,j)=1
         col_fin: do i=j+1,Nc
            b(i,1)=b(i,1) - L(i,j)*b(j,1) 
         end do col_fin
      end do row_fin

   end subroutine for_solve_lu_blk

   subroutine fb_solve_chol(A,b,x)
      real (kind=8) :: A(:,:), b(:,:), x(:,:)
      real (kind=8), allocatable :: wrk(:,:)
      integer (kind=4) :: Nr,Nc

      Nr = size(A,1)
      Nc = size(A,2)
   
      if (Nr /= Nc) then
         print *, "error: compla.f08: fb_solve_chol: input matrix is not square"
         stop
      end if

      allocate(wrk(Nc,Nc))
      x = b
      wrk = A
      call chol(wrk) ! stores R in wrk
      wrk = transpose(wrk)
      call for_solve(wrk, x) !A*x = R'*R*x=b -> R'*y=b
      wrk = transpose(wrk)
      call back_solve(wrk,x) ! Rx=y

   end subroutine fb_solve_chol

   subroutine fb_solve_blk_chol(A,b,x)
      real (kind=8) :: A(:,:), b(:,:), x(:,:)
      real (kind=8), allocatable :: wrk(:,:)
      integer (kind=4) :: Nr,Nc

      Nr = size(A,1)
      Nc = size(A,2)
   
      if (Nr /= Nc) then
         print *, "error: compla.f08: fb_solve_blk_chol: input matrix is not square"
         stop
      end if

      allocate(wrk(Nc,Nc))
      x = b
      wrk = A
      call chol(wrk) ! stores R in wrk
      wrk = transpose(wrk)
      call for_solve_blk(wrk, x) !A*x = R'*R*x=b -> R'*y=b
      wrk = transpose(wrk)
      call back_solve_blk(wrk,x) ! Rx=y
   
   end subroutine fb_solve_blk_chol


   ! This routine is basically for testing purposes
   ! I'd probably want to use the block routine
   subroutine fb_solve_lu(A,b,x,pvec)
      ! A is assumed to be a ``work array''
      real (kind=8) :: A(:,:), b(:,:), x(:,:)
      integer (kind=4), allocatable, intent(out), optional :: pvec(:)

      integer (kind=4) :: Nr,Nc,i
      integer (kind=4), allocatable :: p(:)

      Nr = size(A,1)
      Nc = size(A,2)
   
      if (Nr /= Nc) then
         print *, "error: compla.f08: fb_solve_lu: input matrix is not square"
         stop
      end if

      allocate(p(Nc))
      p = (/ (i,i=1,Nc) /)
      call lu(A,p) ! stores LU in A (note that L,U wont be triangular until we apply the permutation vector)
      call apply_perm_vector(A,p,0) ! now make L,U triangular in memory

      x = b
      call apply_perm_vector(x,p,0) ! permute b ( PAx = LUx = Pb )
 
      call for_solve_lu(A,x)
      call back_solve(A,x)

      if (present(pvec)) pvec = p

   end subroutine fb_solve_lu

   subroutine fb_solve_blk_lu(A,b,x,pvec)
      ! A is assumed to be a ``work array''
      ! pvec is optional argument to output permutation vector
      real (kind=8) :: A(:,:), b(:,:), x(:,:)
      integer (kind=4), allocatable, intent(out), optional :: pvec(:)

      integer (kind=4) :: Nr,Nc,i
      integer (kind=4), allocatable :: p(:)

      Nr = size(A,1)
      Nc = size(A,2)
   
      if (Nr /= Nc) then
         print *, "error: compla.f08: fb_solve_blk_lu: input matrix is not square"
         stop
      end if

      allocate(p(Nc))
      p = (/ (i,i=1,Nc) /)
      call lu(A,p) ! stores LU in A
      call apply_perm_vector(A,p,0)

      x = b
      call apply_perm_vector(x,p,0) ! permute b

      call for_solve_lu_blk(A,x)
      call back_solve_blk(A,x)

      if (present(pvec)) pvec = p
   
   end subroutine fb_solve_blk_lu

   ! }}}


   !!!!!!!!!!!!!!!!!
   ! Vec/Mat Norms !
   !!!!!!!!!!!!!!!!!
   ! {{{
   ! Frobenius matrix norm/vector 2-norm
   function norm_f(A)
      real (kind=8), intent(in) :: A(:,:)
      real (kind=8) :: norm_f

      integer (kind=4) :: i,j,Nr,Nc

      Nr=size(A,1)
      Nc=size(A,2)

      norm_f = 0d0

      row: do j=1,Nc
         col: do i=1,Nr
            norm_f = norm_f + A(i,j)*A(i,j)
         end do col
      end do row

      norm_f = sqrt(norm_f)
   
   end function norm_f

   ! vector p-norm
   function norm_p_vec(A,p)
      real (kind=8), intent(in) :: A(:)
      integer (kind=4), intent(in) :: p
      real (kind=8) :: norm_p_vec

      ! 1 norm (max col sum)
      if (p == 1) then
         norm_p_vec = sum(abs(A(:)))

      ! Infinity norm (max row sum)
      else if (p == 0) then
         norm_p_vec = maxval(abs(A(:)))

      ! 2 norm 
      else if (p == 2) then
         norm_p_vec = norm2(A)

      else
         print *, "error: compla.f08: norm_p -> norm_p_vec: unsupported p-norm", p
         stop

      end if

   end function norm_p_vec

   ! Induced matrix p-norm
   function norm_p_mat(A,p)
      real (kind=8), intent(in) :: A(:,:)
      integer (kind=4), intent(in) :: p
      real (kind=8) :: norm_p_mat

      real (kind=8) :: temp
      integer (kind=4) :: i

      norm_p_mat = 0d0

      ! 1 norm (max col sum)
      if (p == 1) then
         row: do i=1,size(A,2)
            temp = sum(abs(A(:,i)))
            if (temp > norm_p_mat) then
               norm_p_mat = temp
            end if
         end do row

      ! Infinity norm (max row sum)
      else if (p == 0) then
         col: do i=1,size(A,1)
            temp = sum(abs(A(i,:)))
            if (temp > norm_p_mat) then
               norm_p_mat = temp
            end if
         end do col

      ! 2 norm (lol, do you really need the induced 2 norm)
      else
         if (size(A,2) == 1) then ! vector 2-norm
            norm_p_mat = norm2(A)
         else
            stop "2-norm not implemented"
         end if

      end if

   end function norm_p_mat

   function condest_lu(A,wrk,pvec,Tin)
      ! 1-norm condition number estimation
      ! using LU overwritten on A with permutation vector pvec
      !
      ! Example:
      !     A = rand_mat(100,100)
      !     wrk = A
      !     pvec = (/ (i,i=1,100) /)
      !     call lu(wrk,pvec)
      !     call apply_perm_vector(wrk,pvec,0)
      !     print *, condest(A,wrk,pvec)
       
      real (kind=8), intent(in) :: A(:,:), wrk(:,:)
      integer (kind=4), intent(in) :: pvec(:)
      integer (kind=4), optional :: Tin

      real (kind=8) :: condest_lu, temp
      real (kind=8), allocatable :: w(:,:), w_temp(:,:)
      integer (kind=4) :: T,i

      T = 2
      if (present(Tin)) T=Tin

      condest_lu = 0d0

      ! Try to estimate \|A^{-1}w\|_1 / \|w\|_1
      do i=1,T
         w = rand_mat(size(A,2),1)
         w_temp = w

         call apply_perm_vector(w_temp,pvec,0)
         call for_solve_lu_blk(wrk,w_temp)
         call back_solve_blk(wrk,w_temp)
         
         temp = norm_p(w_temp,1) / norm_p(w,1)
         if (temp > condest_lu) condest_lu = temp

      end do

      ! Finally, multiply by \|A\|_1
      condest_lu = condest_lu * norm_p(A,1)

   end function condest_lu



   ! }}}


   !!!!!!!!!!!!!!!!!
   ! Misc Routines !
   !!!!!!!!!!!!!!!!!
   ! {{{
   ! behaves like Octave's tic()
   subroutine tic(t)
      integer, intent(out) :: t
      call system_clock(t)
   end subroutine tic
   
   ! returns time in seconds from now to time described by t
   real function toc_return(t)
      integer, intent(in) :: t
      integer :: now, clock_rate
   
      call system_clock(now, clock_rate)
      toc_return = real(now - t)/real(clock_rate)
   end function toc_return
   
   ! prints time in seconds from now to time described by t
   subroutine toc(t)
      integer, intent(in) :: t
      real (kind=8) :: time
      integer :: now, clock_rate
   
      call system_clock(now, clock_rate)
      time = real(now - t)/real(clock_rate)
      print *,"Elapsed time is ",time," seconds."
   end subroutine toc
   
   ! print a rank 2  array
   subroutine print_array(A)
      real (kind=8), intent(in) :: A(:,:)

      integer (kind=4) :: i,j
      character (len=30) :: rowfmt
   
      write(rowfmt, "(A,I4,A)") "(",size(A,2),"(1X,SS,10Es13.4))"
      print *,
      row_print: do i=1,size(A,1)
         write(*, fmt=rowfmt) (A(i,j), j=1,size(A,2))
      end do row_print
   end subroutine print_array

   ! print a rand 1 array
   subroutine print_vector(A)
      real (kind=8), intent(in) :: A(:)

      integer (kind=4) :: i

      print *,
      row_print: do i=1,size(A,1)
         write(*, fmt="(1X,SS,10Es13.4)") A(i)
      end do row_print
   end subroutine print_vector

   function eye(Nr,Nc)
      integer (kind=4), intent(in) :: Nr,Nc
      real (kind=8), allocatable :: eye(:,:)

      integer (kind=4) :: i

      allocate(eye(Nr,Nc))
      eye = 0
      do i=1,min(Nr,Nc)
         eye(i,i) = 1d0
      end do

   end function eye

   ! k-th subdiagonal of A
   function diag(A,kin)
      real (kind=8), intent(in) :: A(:,:)
      integer (kind=4), intent(in), optional :: kin
      real (kind=8), allocatable :: diag(:)

      integer (kind=4) :: m,n,i,k
      m = size(A,1); n = size(A,2)

      if ( present(kin) ) then
         k = kin
      else
         k = 0
      end if

      allocate(diag(min(m-abs(k),n-abs(k))))

      if ( k >= 0 ) then
         do i=1,min(m-abs(k),n-abs(k))
            diag(i) = A(i+k,i)
         end do
      else
         do i=1,min(m-abs(k),n-abs(k))
            diag(i) = A(i,i-k)
         end do
      end if

   end function diag

   function rand_mat(Nr,Nc)
      integer (kind=4), intent(in) :: Nr, Nc
      real (kind=8), allocatable :: rand_mat(:,:)

      real (kind=8) :: rand
      integer (kind=4) :: i,j

      allocate(rand_mat(Nr,Nc))

      ! Populate with random numbers [0,1]
      row: do j=1,Nc
         col: do i=1,Nr
            call random_number(rand)
            rand_mat(i,j) = rand
         end do col
      end do row
   end function rand_mat

   function rand_spd_mat(N)
      integer (kind=4), intent(in) :: N
      real (kind=8), allocatable :: rand_spd_mat(:,:)
      
      real (kind=8) :: rand, m
      integer (kind=4) :: i,j

      allocate(rand_spd_mat(N,N))
      m=0d0

      ! Populate with random numbers [0,10]
      row: do j=1,N
         col: do i=1,N
            call random_number(rand)
            rand = 10d0*rand
            rand_spd_mat(i,j) = rand

            if (rand > m) m=rand

         end do col
      end do row

      ! Make symmetric
      rand_spd_mat = rand_spd_mat + transpose(rand_spd_mat)

      ! Make diagonally dominant (implies positive definite)
      diag: do j=1,N
         rand_spd_mat(j,j) = rand_spd_mat(j,j) + 10d0*N
      end do diag

   end function rand_spd_mat

  
   subroutine init_random_seed()
     ! Stolen from http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
     ! {{{
     implicit none
     integer, allocatable :: seed(:)
     integer :: i, n, un, istat, dt(8), pid, t(2), s
     integer(8) :: count, tms
   
     call random_seed(size = n)
     allocate(seed(n))
     ! First try if the OS provides a random number generator
     open(newunit=un, file="/dev/urandom", access="stream", &
          form="unformatted", action="read", status="old", iostat=istat)
     if (istat == 0) then
        read(un) seed
        close(un)
     else
        ! Fallback to XOR:ing the current time and pid. The PID is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        call system_clock(count)
        if (count /= 0) then
           t = transfer(count, t)
        else
           call date_and_time(values=dt)
           tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
           t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n >= 3) then
           seed(1) = t(1) + 36269
           seed(2) = t(2) + 72551
           seed(3) = pid
           if (n > 3) then
              seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
           end if
        else
           seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
     end if
     call random_seed(put=seed)
     ! }}}

   end subroutine init_random_seed
   ! }}}


   !!!!!!!!
   ! HDF5 !
   !!!!!!!!
   ! {{{

   ! These routines make the dspace, dset, and write data
   ! they all interface to the 'write_dset' generic routine
   subroutine dwrite_dset_rank0(file_id, scalar, dsetname)
      integer (kind=intk) :: file_id
      real (kind=dblk) :: scalar
      character (len=*) :: dsetname

      integer (kind=intk) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(1)

      dims(1) = 0
      call h5screate_simple_f(0,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,scalar,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)

   end subroutine dwrite_dset_rank0

   subroutine dwrite_dset_rank1(file_id, array, dsetname)
      integer (kind=intk) :: file_id
      real (kind=8) :: array(:)
      character (len=*) :: dsetname

      integer (kind=4) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(1)

      dims(1) = size(array)
      call h5screate_simple_f(1,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)

   end subroutine dwrite_dset_rank1

   subroutine iwrite_dset_rank1(file_id, array, dsetname)
      integer (kind=intk) :: file_id
      integer (kind=intk) :: array(:)
      character (len=*) :: dsetname

      integer (kind=intk) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(1)

      dims(1) = size(array)
      call h5screate_simple_f(1,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_INTEGER,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,array,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)

   end subroutine iwrite_dset_rank1


   subroutine dwrite_dset_rank2(file_id, array, dsetname)
      integer (kind=intk) :: file_id
      real (kind=dblk) :: array(:,:)
      character (len=*) :: dsetname

      integer (kind=intk) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(2)

      dims(1) = size(array,1); dims(2) = size(array,2)
      call h5screate_simple_f(2,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)
   
      end subroutine dwrite_dset_rank2
 

   ! }}}

end module compla

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
