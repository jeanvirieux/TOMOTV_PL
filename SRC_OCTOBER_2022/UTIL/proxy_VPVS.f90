!=======================================================================
!
!  compute Vp/VS and VP*VS et le rapport
!   
!***********************************************************************
program proxy_VPVS
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: ratio,pro,div
  !     integer pmodforwp,pmodforws
  !----------------------------- variables
  integer(kind=4) :: n1inv,n2inv,n3inv
  real(kind=4) :: h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  integer(kind=4) :: flog,finput                      ! flux entree et sortie

  !  character(len=132) name_acquiP,name_acquiS       ! acquisition files
  !  character(len=132) src_inv,src_tab               !
  !  character(len=132) name_modelp,name_forwardp     ! filename for model and FMM for P wave
  !  character(len=132) name_models,name_forwards     ! filename for model and FMM for S wave

  character(len=1)                   :: carac
  integer(kind=4) :: chx
  real(kind=4) :: rat_min,rat_max,pro_min,pro_max

  flog=51
  open(flog,file='flog.model')
  finput=flog+1
  ! #################################### choix des ondes
  open(finput,file='../model.head',status='old')
  read(finput,'(a)') carac
  read(finput,*) chx
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
  

  ! #####################################  arrays
  ALLOCATE(modinvp(n1inv,n2inv,n3inv))
  ALLOCATE(modinvs(n1inv,n2inv,n3inv))
  ALLOCATE(ratio(n1inv,n2inv,n3inv))
  ALLOCATE(pro(n1inv,n2inv,n3inv))
  ALLOCATE(div(n1inv,n2inv,n3inv))


  ! #####################################  io
  open(7,file='modelP.old',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  open(8,file='modelS.old',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  open(9,file='ratioPS',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  open(10,file='productPS',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  open(11,file='divisionPS',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  ! #####################################  reading the inverse model
  read(7,rec=1) modinvp             ! true grid in P
  ! ----------- if we have selected P & S velocities simultaneously
  read(8,rec=1) modinvs             ! true grid in S

  ratio(:,:,:)=modinvp(:,:,:)/modinvs(:,:,:)*1000.*10.
  pro(:,:,:)=log(modinvp(:,:,:)*modinvs(:,:,:))*1000.

  div(:,:,:)=ratio(:,:,:)/pro(:,:,:)*1000.

  rat_min=minval(ratio)
  rat_max=maxval(ratio)
  pro_min=minval(pro)
  pro_max=maxval(pro)

  write(*,*) ' ratio min max ',rat_min,rat_max
  write(*,*) ' product min max ',pro_min,pro_max
  
  write(9,rec=1) ratio
  write(10,rec=1) pro
  write(11,rec=1) div
  
  close(7)
  close(8)
  close(9)
  close(10)
  close(11)
  

end program proxy_VPVS






