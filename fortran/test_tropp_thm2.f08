program test_tropp_thm2
   
   use compla
   use omp
   
   implicit none

   !character (len=255), parameter :: savefile = "../data/tropp_fig1_data_fort.h5"
   integer (kind=intk), parameter :: d=256_intk
   real (kind=dblk), parameter    :: delta = 0.35_dblk
   real (kind=dblk), parameter    :: K=5.0_dblk
   integer (kind=intk), parameter :: num_sigs = 1000
   integer (kind=intk), parameter :: m = 20
   integer (kind=intk) :: N

   ! loop vars
   real (kind=dblk), allocatable :: Phi(:,:),s(:),v(:),s_hat(:)
   real (kind=dblk) :: mu_Phi,Sigma_Phi
   integer (kind=intk) :: i_m,i_N,sig_ind,num_recovered
   logical :: realloc

   ! N from Thm 2 of Tropp 2007
   N = ceiling(K*dble(m)*log(dble(d)/delta))
   N = 256_intk

   call init_random_seed()

   allocate(s(d),s_hat(d))

   ! generate signals and try to recover them
   allocate(Phi(N,d),v(N))
  
   num_recovered = 0.0_dblk
   mu_Phi = 0.0_dblk
   Sigma_Phi = 1.0_dblk/dble(N)
   call mvnrml_rand_mat(Phi,mu_Phi,Sigma_Phi)

   realloc = .true. ! flag to reallocate mem. for vars in omp_algo
   do sig_ind=1,num_sigs
      !call gen_sig(s,m_vec(i_m))
      call gen_sig(s,m,1.0_dblk)

      ! v = 1*Phi*s + 0*v
      call dgemv('N',N,d,1.0_dblk,Phi,N,&
         s,1,0.0_dblk,v,1)

      ! Do OMP recovery
      ! omp_algo sets realloc to .false. if it allocates
      call omp_algo(s_hat,m,Phi,v,realloc)

      num_recovered = num_recovered + check_recovery(s,s_hat)
      !print *, dble(sig_ind)/dble(num_sigs)*100.0_dblk,num_recovered,sig_ind
      print *, dble(sig_ind)/dble(num_sigs)*100.0_dblk,sig_ind - num_recovered

   end do

   print *, "N = ",N
   write (*,"(A,I5,A,I5)"), "Number Recovered: ",num_recovered," of ",num_sigs

   deallocate(Phi,v)
   deallocate(s,s_hat)

   contains
   
end program test_tropp_thm2

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
