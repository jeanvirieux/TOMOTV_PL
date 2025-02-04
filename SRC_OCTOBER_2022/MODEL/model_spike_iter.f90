!=======================================================================
!
!***********************************************************************
!
!***********************************************************************
program model_spike_iter
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modspike

  !----------------------------- variables
  real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  integer(kind=4) n1inv,n2inv,n3inv,chx   ! dimensions for inversion grid 
  integer(kind=4) finput             ! flux entree et sortie

  integer(kind=4)                    :: i1,i2,i3
  character(len=1)                   :: carac
  character(len=4)                   :: etiq

  integer(kind=4) :: iopt
  real(kind=4) :: x1,x2,x3
  real(kind=4) :: xcent1,xcent2,xcent3,l1,l2,l3

  real(kind=8) :: arg
  real(kind=4) :: vpert

  real(kind=4) :: xstart1,xend1,xstep1
  real(kind=4) :: xstart2,xend2,xstep2
  real(kind=4) :: xstart3,xend3,xstep3

  integer(kind=4) :: nrep1,nrep2,nrep3
  integer(kind=4) :: in1,in2,in3,iter

  finput=51
  ! #################################### choix des ondes
  open(finput,file='model.head',status='old')
  read(finput,'(a)') carac
  read(finput,*) chx
  ! ####################################  input
  read(finput,'(a)') carac
  read(finput,*) x1inv,x2inv,x3inv
  read(finput,'(a)') carac
  read(finput,*) n1inv,n2inv,n3inv
  read(finput,'(a)') carac
  read(finput,*) h1inv,h2inv,h3inv
  close(finput)


  ALLOCATE(modspike(n1inv,n2inv,n3inv))
  ! #####################################  array P
  ALLOCATE(modinvp(n1inv,n2inv,n3inv))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then                  !###  array S
     ALLOCATE(modinvs(n1inv,n2inv,n3inv))
  endif

  ! #####################################  io
  open(7,file='modelP',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  !------------ if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     open(9,file='modelS',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  endif

  ! #####################################  reading the inverse model
  read(7,rec=1) modinvp             ! true grid in P
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     read(9,rec=1) modinvs             ! true grid in S
  endif
  close(7)
  close(9)

  write(*,*) ' iopt=0 relative or iopt=1 absolute values '
  read(*,*) carac
  read(*,*) iopt

  write(*,*) ' enter spreading of gaussian spike l1,l2,l3 '
  read(*,*) carac
  read(*,*) l1,l2,l3

  write(*,*) ' starting P models design'
!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!! P models

  read(*,*) carac
  read(*,*) vpert
  read(*,*) xstart1,xend1,nrep1
  read(*,*) xstart2,xend2,nrep2
  read(*,*) xstart3,xend3,nrep3
  if(nrep1 > 1) then
     xstep1=(xend1-xstart1)/(nrep1-1)
  else
     xstep1=0.
  endif
  if(nrep2 > 1) then
     xstep2=(xend2-xstart2)/(nrep2-1)
  else
     xstep2=0.
  endif
  if(nrep3 > 1) then
     xstep3=(xend3-xstart3)/(nrep3-1)
  else
     xstep3=0.
  endif
!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!! building spike models
  iter =0 
  do in3=1,nrep3
     xcent3=xstart3+xstep3*(in3-1)
     do in2=1,nrep2
        xcent2=xstart2+xstep2*(in2-1)
        do in1=1,nrep1
           xcent1=xstart1+xstep1*(in1-1)

           iter=iter+1
           write(etiq,'(i4.4)') iter
           write(*,*) ' P iter ',iter

           do i3=1,n3inv
              x3=x3inv+(i3-1)*h3inv
              do i2=1,n2inv
                 x2=x2inv+(i2-1)*h2inv
                 do i1=1,n1inv
                    x1=x1inv+(i1-1)*h1inv
                    arg=((x1-xcent1)/l1)**2+((x2-xcent2)/l2)**2+((x3-xcent3)/l3)**2
                    if(iopt == 1) then
                       modspike(i1,i2,i3)=modinvp(i1,i2,i3)+vpert*exp(-arg)   ! valeur vraie
                    else
                       modspike(i1,i2,i3)=exp(-arg)*vpert*20.  ! for color scaling
                    endif
                 enddo
              enddo
           enddo

           write(*,*) 'P valeur min ',minval(modspike),' P valeur max ',MAXVAL(modspike)
           open(7,file='modelP.spike.'//etiq,access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
           write(7,rec=1) modspike
           close(7)


        enddo
     enddo
  enddo

  write(*,*) ' end for P spike model '

  if(chx == 2) then

!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!! S  models

     read(*,*) carac
     read(*,*) vpert
     read(*,*) xstart1,xend1,nrep1
     read(*,*) xstart2,xend2,nrep2
     read(*,*) xstart3,xend3,nrep3
     if(nrep1 > 1) then
        xstep1=(xend1-xstart1)/(nrep1-1)
     else
        xstep1=0.
     endif
     if(nrep2 > 1) then
        xstep2=(xend2-xstart2)/(nrep2-1)
     else
        xstep2=0.
     endif
     if(nrep3 > 1) then
        xstep3=(xend3-xstart3)/(nrep3-1)
     else
        xstep3=0.
     endif

!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!! building spike models
     iter =0 

     do in3=1,nrep3
        xcent3=xstart3+xstep3*(in3-1)
        do in2=1,nrep2
           xcent2=xstart2+xstep2*(in2-1)
           do in1=1,nrep1
              xcent1=xstart1+xstep1*(in1-1)

              iter=iter+1
              write(etiq,'(i4.4)') iter
              write(*,*) ' S iter ',iter

              do i3=1,n3inv
                 x3=x3inv+(i3-1)*h3inv
                 do i2=1,n2inv
                    x2=x2inv+(i2-1)*h2inv
                    do i1=1,n1inv
                       x1=x1inv+(i1-1)*h1inv


                       arg=((x1-xcent1)/l1)**2+((x2-xcent2)/l2)**2+((x3-xcent3)/l3)**2
                       if(iopt == 1) then
                          modspike(i1,i2,i3)=modinvp(i1,i2,i3)+vpert*exp(-arg)  ! la valeur vraie
                       else
                          modspike(i1,i2,i3)=exp(-arg)*vpert*20.  ! for color scaling
                       endif
                    enddo
                 enddo
              enddo
              write(*,*) 'S valeur min ',minval(modspike),' S valeur max ',MAXVAL(modspike)
              open(7,file='modelS.spike.'//etiq,access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
              write(7,rec=1) modspike
              close(7)


           enddo
        enddo
     enddo

     write(*,*) ' end for S spike model '

  endif

end program model_spike_iter






