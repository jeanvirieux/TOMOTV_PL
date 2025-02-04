!=======================================================================
!
!***********************************************************************
!
!***********************************************************************
program model_spike
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modspike

  !----------------------------- variables
  real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  integer(kind=4) n1inv,n2inv,n3inv,chx   ! dimensions for inversion grid 
  integer(kind=4) flog,finput             ! flux entree et sortie

  integer(kind=4)                    :: i1,i2,i3
  character(len=1)                   :: carac

  integer(kind=4) :: iopt
  real(kind=4) :: x1,x2,x3
  real(kind=4) :: xcent1,xcent2,xcent3,l1,l2,l3

  real(kind=8) :: arg
  real(kind=4) :: vpert


  flog=51
  open(flog,file='flog.model_spike')
  finput=flog+1
  ! #################################### choix des ondes
  open(finput,file='model.head',status='old')
  read(finput,'(a)') carac
  read(finput,*) chx
  write(flog,*) ' flag 1 pour ondes P, 2 pour ondes P et S',chx
  ! ####################################  input
  read(finput,'(a)') carac
  read(finput,*) x1inv,x2inv,x3inv
  write(flog,*) ' (x1,x2,x3) origin of the inversion grid ',x1inv,x2inv,x3inv
  read(finput,'(a)') carac
  read(finput,*) n1inv,n2inv,n3inv
  write(flog,*) ' (n1,n2,n3) inversion grid ',n1inv,n2inv,n3inv
  read(finput,'(a)') carac
  read(finput,*) h1inv,h2inv,h3inv
  write(flog,*) ' (h1,h2,h3) grid sampling ', h1inv,h2inv,h3inv
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
  
! #################################### choix des ondes
  open(finput,file='spike.head',status='old')
  write(*,*) ' iopt=0 relative or iopt=1 absolute values '
  read(finput,*) carac
  read(finput,*) iopt

  write(*,*) ' enter spreading of gaussian spike l1,l2,l3 '
  read(finput,*) carac
  read(finput,*) l1,l2,l3

  write(*,*) ' starting P models design'
!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!! P models
  read(finput,*) carac
  read(finput,*) xcent1,xcent2,xcent3,vpert

  do i3=1,n3inv
     x3=x3inv+(i3-1)*h3inv
     do i2=1,n2inv
        x2=x2inv+(i2-1)*h2inv
        do i1=1,n1inv
           x1=x1inv+(i1-1)*h1inv
           arg=((x1-xcent1)/l1)**2+((x2-xcent2)/l2)**2+((x3-xcent3)/l3)**2
           if(iopt == 1) then
              modspike(i1,i2,i3)=modinvp(i1,i2,i3)+vpert*exp(-arg)
           else
              modspike(i1,i2,i3)=exp(-arg)*vpert*15.  ! for color scaling
           endif
        enddo
     enddo
  enddo
  write(*,*) ' P valeur max ',MAXVAL(modspike)
  open(7,file='modelP.spike',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  write(7,rec=1) modspike
  close(7)

  write(*,*) ' end for P spike model '

  if(chx == 2) then

!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!! S  models
     read(finput,*) carac
     read(finput,*) xcent1,xcent2,xcent3,vpert
     close(finput)

     do i3=1,n3inv
        x3=x3inv+(i3-1)*h3inv
        do i2=1,n2inv
           x2=x2inv+(i2-1)*h2inv
           do i1=1,n1inv
              x1=x1inv+(i1-1)*h1inv
              arg=((x1-xcent1)/l1)**2+((x2-xcent2)/l2)**2+((x3-xcent3)/l3)**2
              if(iopt == 1) then
                 modspike(i1,i2,i3)=modinvs(i1,i2,i3)+vpert*exp(-arg)
              else
                 modspike(i1,i2,i3)=exp(-arg)*vpert*20.  ! for color scaling
              endif
           enddo
        enddo
     enddo
     write(*,*) ' S valeur max ',MAXVAL(modspike)
     open(7,file='modelS.spike',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(7,rec=1) modspike
     close(7)

     write(*,*) ' end for S spike model '

  endif

  close(flog)

end program model_spike






