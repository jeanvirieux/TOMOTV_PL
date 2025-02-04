!=======================================================================
! difference entre 2 modeles a ajouter a un modele courant
!***********************************************************************
!
!***********************************************************************
program model_diff_sum
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:) :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:) :: modspikep,modspikes
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:) :: moddiffp,moddiffs

  !----------------------------- variables
  real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  integer(kind=4) n1inv,n2inv,n3inv,chx   ! dimensions for inversion grid 
  integer(kind=4) flog,finput             ! flux entree et sortie

  character(len=1)                   :: carac

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
  ALLOCATE(moddiffp(n1inv,n2inv,n3inv))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then                  !###  array S
     ALLOCATE(modspikes(n1inv,n2inv,n3inv))
     ALLOCATE(modinvs(n1inv,n2inv,n3inv))
     ALLOCATE(moddiffs(n1inv,n2inv,n3inv))
  endif
  
  open(7,file='modelP.dif',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  read(7,rec=1) moddiffp
  close(7)
  write(*,*) 'P max and min',maxval(moddiffp),minval(moddiffp)
  if(chx.eq.2) then
     open(9,file='modelS.dif',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     read(9,rec=1) moddiffs
     close(9)
     write(*,*) 'S max and min',maxval(moddiffs),minval(moddiffs)
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
  moddiffp(:,:,:)=moddiffp(:,:,:)+max(0.,modspikep(:,:,:)-modinvp(:,:,:))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     read(9,rec=1) modspikes             ! true grid in S
     read(99,rec=1) modinvs
     close(9)
     close(99)
     moddiffs(:,:,:)=moddiffs(:,:,:)+max(0.,modspikes(:,:,:)-modinvs(:,:,:))
  endif
  open(7,file='modelP.dif',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  write(7,rec=1) moddiffp
  close(7)
  write(*,*) 'P max and min 2',maxval(moddiffp),minval(moddiffp)
  if(chx.eq.2) then
     open(9,file='modelS.dif',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(9,rec=1) moddiffs
     close(9)
     write(*,*) 'S max and min 2',maxval(moddiffs),minval(moddiffs)
  endif

  close(flog)

end program model_diff_sum






