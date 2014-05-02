program tropp_fig1_comp
   
   use compla
   use omp
   
   implicit none

   character (len=255), parameter :: savefile = "../data/tropp_fig1_data_fort.h5"
   integer (kind=intk), parameter :: d=256_intk
   real (kind=dblk), parameter    :: delta = 0.1_dblk
   integer (kind=intk), parameter :: K=5_intk
   integer (kind=intk), parameter :: num_sigs = 100
   integer (kind=intk), parameter, dimension(5) :: m_vec = (/ 4,12,20,28,36 /)

   ! pre-loop vars
   integer (kind=intk), allocatable :: N_vec(:)
   real (kind=dblk), allocatable :: percent_recovered(:,:)

   ! loop vars
   real (kind=dblk), allocatable :: Phi(:,:),s(:),v(:),s_hat(:)
   real (kind=dblk) :: mu_Phi,Sigma_Phi
   integer (kind=intk) :: i_m,i_N,sig_ind
   logical :: realloc


   call init_random_seed()

   ! qrappendcol cannot handle short and wide matrices
   ! so we must start at the maximum value of m (corresonds to square Phi_t)
   N_vec = ceiling(linspace(36.0_dblk,dble(d),ceiling(real((real(d)-1_dblk)/5_dblk))))

   allocate(percent_recovered(size(N_vec),size(m_vec)))
   allocate(s(d),s_hat(d))

   percent_recovered = 0_dblk
  
   ! generate signals and try to recover them
   do i_m=1,size(m_vec)
      do i_N=1,size(N_vec)

         ! print some stuff 
         write (*,"(A,I3,A,I3,A,I3)"),"m = ",m_vec(i_m),"  N = ",N_vec(i_N)," of ",N_vec(size(N_vec))   

         allocate(Phi(N_vec(i_N),d),v(N_vec(i_N)))
         
         mu_Phi = 0.0_dblk
         Sigma_Phi = 1_dblk/dble(N_vec(i_N))
         call mvnrml_rand_mat(Phi,mu_Phi,Sigma_Phi)

         realloc = .true. ! flag to reallocate mem. for vars in omp_algo
         do sig_ind=1,num_sigs
            !call gen_sig(s,m_vec(i_m))
            call gen_sig(s,m_vec(i_m),1.0_dblk)

            ! v = 1*Phi*s + 0*v
            call dgemv('N',N_vec(i_N),d,1.0_dblk,Phi,N_vec(i_N),&
               s,1,0.0_dblk,v,1)

            ! Do OMP recovery
            ! omp_algo sets realloc to .false. if it allocates
            call omp_algo(s_hat,m_vec(i_m),Phi,v,realloc)

            percent_recovered(i_N,i_m) = percent_recovered(i_N,i_m) + &
               dble(check_recovery(s,s_hat))

         end do

         print *, percent_recovered(i_N,i_m) / dble(num_sigs) * 100.0_dblk
         deallocate(Phi,v)
         stop('stuff')
      end do

      ! save data
      percent_recovered = percent_recovered / dble(num_sigs) * 100.0_dblk
      call save_data_fig1(savefile,N_vec,m_vec,percent_recovered,d)
   end do

   percent_recovered = percent_recovered / dble(num_sigs) * 100.0_dblk
   call save_data_fig1(savefile,N_vec,m_vec,percent_recovered,d)

   deallocate(percent_recovered,s,s_hat)

   contains
   
   ! save the stuff to the thing
   subroutine save_data_fig1(savefile,N_vec,m_vec,percent_recovered,d)
   ! {{{
      character (len=*) :: savefile
      integer (kind=intk) :: N_vec(:),m_vec(:),d
      real (kind=dblk) :: percent_recovered(:,:)
     
      ! HDF5 stuff
      integer (kind=intk) :: h5error, file_id

      call h5open_f(h5error)
      call h5fcreate_f(savefile, H5F_ACC_TRUNC_F, file_id, h5error)
 
      call write_dset(file_id, dble(N_vec), "N_vec") ! lol, hack hack hack
      call write_dset(file_id, dble(m_vec), "m_vec")
      call write_dset(file_id, percent_recovered, "percent_recovered")
      call write_dset(file_id, dble(d), "d")
 
      call h5fclose_f(file_id,h5error)
      call h5close_f(h5error)

   end subroutine save_data_fig1
   ! }}}


end program tropp_fig1_comp

! vim: set ts=3 sw=3 sts=3 et :
! vim: foldmarker={{{,}}}
