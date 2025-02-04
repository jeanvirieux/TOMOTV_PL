!=======================================================================
! difference entre 2 modeles
!***********************************************************************
!
!***********************************************************************
program model_diff
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:) :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:) :: modspikep,modspikes

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
  open(flog,file='flog.model_diff')
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


  ! #####################################  array P
  ALLOCATE(modspikep(n1inv,n2inv,n3inv))
  ALLOCATE(modinvp(n1inv,n2inv,n3inv))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then                  !###  array S
     ALLOCATE(modspikes(n1inv,n2inv,n3inv))
     ALLOCATE(modinvs(n1inv,n2inv,n3inv))
  endif

  ! #####################################  io
  open(7,file='modelP',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  open(77,file='modelP.ref',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  !------------ if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     open(9,file='modelS',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
     open(99,file='modelS.ref',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  endif

  ! #####################################  reading the inverse model
  read(7,rec=1) modspikep             ! true grid in P
  read(77,rec=1) modinvp
  close(7)
  close(77)
  modinvp(:,:,:)=max(0.,modspikep(:,:,:)-modinvp(:,:,:))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     read(9,rec=1) modspikes             ! true grid in S
     read(99,rec=1) modinvs
     close(9)
     close(99)
     modinvs(:,:,:)=max(0.,modspikes(:,:,:)-modinvs(:,:,:))
  endif
  open(7,file='modelP.dif',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  write(7,rec=1) modinvp
  close(7)
  write(*,*) maxval(modinvp),minval(modinvp)
  if(chx.eq.2) then
     open(9,file='modelS.dif',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(9,rec=1) modinvs
     close(9)
  endif

  close(flog)

end program model_diff






