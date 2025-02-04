subroutine sub_stopping_criterion(n1,n2,n3,delta,epsilon,x_ref,x_dns,g,convergence_test,gap)

  implicit none

  !VARIABLE DECLARATIONS
  integer :: n1,n2,n3,i1,i2,i3    
  real,dimension(n1,n2,n3) :: x_ref,x_dns
  real,dimension(n1,n2,n3) :: g
  real :: delta,epsilon,Dijkx1,Dijkx2,Dijkx3,norm_Dijkx,norm_g,gap
  logical :: convergence_test


  gap=0.

  do i3=1,n3-1
     do i2=1,n2-1
        do i1=1,n1-1       
           Dijkx1=x_dns(i1+1,i2,i3)-x_dns(i1,i2,i3)
           Dijkx2=x_dns(i1,i2+1,i3)-x_dns(i1,i2,i3)
           Dijkx3=x_dns(i1,i2,i3+1)-x_dns(i1,i2,i3)

           norm_Dijkx=sqrt(Dijkx1**2+Dijkx2**2+Dijkx3**2)
           gap=gap+norm_Dijkx 
        enddo
     enddo
  enddo

  call normL2(n1,n2,n3,g,norm_g)
  gap=gap+delta*norm_g

  do i3=1,n3
     do i2=1,n2
        do i1=1,n1
           gap=gap-x_ref(i1,i2,i3)*g(i1,i2,i3)
        enddo
     enddo
  enddo


  !write(*,*) 'gap vaut :', gap
  convergence_test=(gap.le.epsilon)

end subroutine sub_stopping_criterion

subroutine normL2(n1,n2,n3,x,norm_x)

  implicit none

  !IN
  integer :: n1,n2,n3
  real,dimension(n1,n2,n3) :: x
  !IN/OUT
  real :: norm_x

  !Local variable
  integer :: i1,i2,i3

  norm_x=0.
  do i3=1,n3
     do i2=1,n2
        do i1=1,n1
           norm_x=norm_x+x(i1,i2,i3)**2
        enddo
     enddo
  enddo
  norm_x=sqrt(norm_x)

end subroutine normL2

subroutine norm_inf(n1,n2,n3,x,norm_x)
  
  implicit none
 
  !IN
  integer :: n1,n2,n3
  real,dimension(n1,n2,n3) :: x
  !IN/OUT
  real :: norm_x
   
  norm_x=maxval(abs(x(:,:,:)))

end subroutine norm_inf


subroutine sub_normalize(n1,n2,n3,x_ns,min_xns,norm_inf_xns)
  
  implicit none
  
  !IN
  integer :: n1,n2,n3
  !IN/OUT
  real,dimension(n1,n2,n3) :: x_ns 
  real :: min_xns,norm_inf_xns
  
  min_xns=minval(x_ns(:,:,:))
  x_ns(:,:,:)=x_ns(:,:,:)-min_xns
  call norm_inf(n1,n2,n3,x_ns,norm_inf_xns)
  x_ns(:,:,:)=x_ns(:,:,:)/(norm_inf_xns)
  
  write(*,*) 'MIN VAL OF INITIAL  NOISY IMAGE IS ', min_xns
  write(*,*) 'NORM INF INITIAL NOISY IMAGE IS ', norm_inf_xns
  

end subroutine sub_normalize

!===============================================================
! final scaling
!===============================================================
subroutine sub_post(n1,n2,n3,x_ns,x_dns,min_xns,norm_inf_xns)
  
  implicit none
  
  !IN
  integer :: n1,n2,n3
  !IN/OUT
  real,dimension(n1,n2,n3) :: x_ns,x_dns
  real :: min_xns,norm_inf_xns
  
  !RE-NORMALIZE AND RE-TRANSLATE NOISY AND DENOISED IMAGE
  x_ns(:,:,:)=x_ns(:,:,:)*norm_inf_xns  
  x_ns(:,:,:)=x_ns(:,:,:)+min_xns
  x_dns(:,:,:)=x_dns(:,:,:)*norm_inf_xns  
  x_dns(:,:,:)=x_dns(:,:,:)+min_xns
  
  
end subroutine sub_post
!*********************************************************************
!                          DENOISING ROUTINE                         !
!--------------------------------------------------------------------!
! Implementation of the TV denooising algorithm proposed by Hansen   !
! The algorithm is based on the first-order Nesterov algorihthm      !
!--------------------------------------------------------------------!

subroutine sub_denoise_TV(n1,n2,n3,epsilon,delta,epsilon_noise,x_ns,x_dns)

  implicit none

  !VARIABLE DECLARATIONS
  !Image dimensions (n1=first, n2=second, n3=third) IN
  integer :: n1,n2,n3    
  !Noisy image IN
  real,dimension(n1,n2,n3) :: x_ns
  !Denoised image IN/OUT
  real,dimension(n1,n2,n3) :: x_dns
  !Precision and noise level IN
  real :: epsilon,delta

  !LOCAL VARIABLES
  integer :: i1,i2,i3,cpt_iter,cpt_print
  real :: mu,Lmu,Dijkx1,Dijkx2,Dijkx3,norm_Dijkx,norm_temp,norm_w,den,Ak,alpha,tk
  real :: uijk1,uijk2,uijk3,epsilon_true,R,epsilon_noise
  logical :: convergence_test
  real,dimension(:,:,:),allocatable :: x_ref,g,temp,w,y,z
!  real,dimension(:,:,:,:),allocatable :: u
  real ::n1n2n3,gap,gap_prev

  R = maxval(x_ns(:,:,:))
  epsilon_true=R*n1*n2*epsilon
  n1n2n3=float(n1*n2*n3)
  !=================================== jean db noise 25. 15dB
  delta=delta*sqrt(n1n2n3)*epsilon_noise
  mu=epsilon_true/n1n2n3  
  Lmu=8./mu
  
  write(*,*) 'R is :',R
  write(*,*) 'epsilon_true is :',epsilon_true
  write(*,*) ' mu is :',mu
  write(*,*) ' Lmu is :',Lmu
  
  convergence_test=.false.

  allocate(g(n1,n2,n3))
!  allocate(u(n1,n2,n3,2))
  allocate(temp(n1,n2,n3),w(n1,n2,n3),y(n1,n2,n3),z(n1,n2,n3))
  allocate(x_ref(n1,n2,n3))
  g(:,:,:)=0.
!  u(:,:,:,:)=0.
  x_ref(:,:,:)=x_ns(:,:,:)
  x_dns(:,:,:)=x_ref(:,:,:)
  w(:,:,:)=0.
  y(:,:,:)=0.
  z(:,:,:)=0.
  Ak=0.5
  cpt_iter=0
  cpt_print=0

  do while (.not.convergence_test)

     !STEP 1 : COMPUTE g
     g(:,:,:)=0.
     !open(10,file='test')
     do i3=1,n3-1
        do i2=1,n2-1
           do i1=1,n1-1
              Dijkx1=x_dns(i1+1,i2,i3)-x_dns(i1,i2,i3)
              Dijkx2=x_dns(i1,i2+1,i3)-x_dns(i1,i2,i3)
              Dijkx3=x_dns(i1,i2,i3+1)-x_dns(i1,i2,i3)
              !write(10,*) 'Dijx1 ',Dijx1,x_dns(i1+1,i2),x_dns(i1,i2)
              !write(10,*) 'Dijx2 ',Dijx2,x_dns(i1,i2+1),x_dns(i1,i2)
              norm_Dijkx=sqrt(Dijkx1**2+Dijkx2**2+Dijkx3**2)           
              den=max(mu,norm_Dijkx)           
              uijk1=Dijkx1/den
              uijk2=Dijkx2/den
              uijk3=Dijkx3/den
              g(i1,i2,i3)=g(i1,i2,i3)-uijk1-uijk2-uijk3           
              g(i1+1,i2,i3)=g(i1+1,i2,i3)+uijk1
              g(i1,i2+1,i3)=g(i1,i2+1,i3)+uijk2
              g(i1,i2,i3+1)=g(i1,i2,i3+1)+uijk3
           enddo
        enddo
     enddo

     !close(10)
!     open(60,file='test_g',access='direct',recl=4*n1*n2*n3)
!     write(60,rec=1) g(:,:,:)
!     close(60)

     !STEP 2 : COMPUTE y     
     temp(:,:,:)=Lmu*(x_dns(:,:,:)-x_ref(:,:,:))-g(:,:,:)
     call normL2(n1,n2,n3,temp,norm_temp)         
     den=max(Lmu,norm_temp/delta)    
     y(:,:,:)=temp(:,:,:)/den+x_ref(:,:,:)

     !STEP 3 : COMPUTE z
     w(:,:,:)=w(:,:,:)+0.5*((cpt_iter+1)**2)*g(:,:,:)    
     call normL2(n1,n2,n3,w,norm_w)
     den=max(Lmu,norm_w/delta)
     z(:,:,:)=-1.0*w(:,:,:)/den+x_ref(:,:,:)

     !STEP 4 : update x_dns     
     alpha=0.5*((cpt_iter+2)**2)
     Ak=Ak+alpha
     tk=alpha/Ak
     x_dns(:,:,:)=tk*z(:,:,:)+(1.-tk)*y(:,:,:)

     !Compute previous gap
     if(cpt_iter.eq.0) then
        call sub_stopping_criterion(n1,n2,n3,delta,epsilon_true,x_ref,x_dns,g,convergence_test,gap_prev)
     else
        gap_prev=gap
     endif
     call sub_stopping_criterion(n1,n2,n3,delta,epsilon_true,x_ref,x_dns,g,convergence_test,gap)

     !Update cpt_iter
     cpt_iter=cpt_iter+1
     cpt_print=cpt_print+1
     if(cpt_print.eq.1) then
        write(*,*) 'Iteration :', cpt_iter
        write(*,*) 'Gap : ',gap
        write(*,*) 'den_z is : ',den
        cpt_print=0
     endif
     !convergence_test=(cpt_iter.ge.10)
     convergence_test=(cpt_iter.ge.50)
     !convergence_test=(cpt_iter.ge.5)
     !convergence_test=(gap>gap_prev)
  enddo

  deallocate(g)
  deallocate(temp,w,y,z)
  deallocate(x_ref)

end subroutine sub_denoise_TV


!*********************************************************************
!                          DENOISING CODE                            !
!--------------------------------------------------------------------!
! Implementation of the TV denooising algorithm proposed by Hansen   !
! The algorithm is based on the first-order Nesterov algorihthm      !
!--------------------------------------------------------------------!
!
!  epsilon ... convergence
!  delta (tau in the paper) ... around 1
!  epsilon_noise ... related to noise level (logarithmique scale)
!
!  cplx_wfld not introduced ... to be checked
!=====================================================================
subroutine nesterov(x_ns,x_dns,n1,n2,n3,epsilon,delta,epsilon_noise)
  implicit none

  !VARIABLE DECLARATIONS
  !Image dimensions (n1=first, n2=second, n3=third)
  integer :: n1,n2,n3
  !Noisy image
  real,dimension(n1,n2,n3) :: x_ns
  !Denoised image
  real,dimension(n1,n2,n3) :: x_dns
  !Precision
  real :: epsilon,delta,epsilon_noise,norm_inf_xns,min_xns
!  !Option
!  integer :: cplx_wfld
  
  !IF COMPLEX WAVEFIELD CONVERT IT IN A REAL NOISY IMAGE
  ! frozen   subroutine should be checked as not tested obviously

!  if(cplx_wfld) then     
!    write(*,*) ' complex field ... conversion '
!     call sub_wfld_to_image(filename,n1,n2,x_ns)
!  else
!     open(60,file=filename,access='direct',recl=4*n1*n2)
!     read(60,rec=1) x_ns(:,:)
!     close(60)
!  endif
    
  !NORMALIZE AND TRANSLATE IMAGE
  call sub_normalize(n1,n2,n3,x_ns,min_xns,norm_inf_xns)
    
  !DENOISE IMAGE
  call sub_denoise_TV(n1,n2,n3,epsilon,delta,epsilon_noise,x_ns,x_dns)
  
  !POST_TREATMENT  
  call sub_post(n1,n2,n3,x_ns,x_dns,min_xns,norm_inf_xns)
  
!  !IF COMPLEX WAVEFIELD CONVERT BACK TO WAVEFIELD 
!  if(cplx_wfld) then     
!     call sub_image_to_wfld(n1,n2,x_dns)
!  endif

end subroutine nesterov



  
