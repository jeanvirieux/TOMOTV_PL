!***********************************************************************
! module : model.interpol
! ==============
! purpose : interpolation the inversion model into a new inversion model
!           with different spatial discretization
! This program builds as well files "timefd.par" and "inversion.par"
!
! One may select only P waves or P & S waves. Grids for P and S velocities
! are identical.
!
! model description : the model is described on a rectangular grid for
!           the moment ..... 
!                      algorithms are selected for possible extension to
!                      unstructured finite element grids
!
! the log file name is "flog.model.tomo.interpol" 
!
!-----------------------------------------------------------------------
! rebuilt by Jean Virieux in fortran 90 (summer 2009)
!***********************************************************************
program model_interpol
  use s_interpol_mod
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modforwp,modforws
  !     integer pmodforwp,pmodforws
  !----------------------------- variables
  real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  real(kind=4) h1,h2,h3,x1_fwd,x2_fwd,x3_fwd       ! spatial step for FMM and origin position
  integer(kind=4) n1inv,n2inv,n3inv,n1,n2,n3,chx   ! dimensions for inversion grid and for FMM grid
  integer(kind=4) flog                             ! output flux
  character(len=132) name_modelp,name_forwardp     ! filename for model and FMM for P wave
  character(len=132) name_models,name_forwards     ! filename for model and FMM for S wave

  real(kind=4) :: vmin,vmax

  integer(kind=4) :: iopt=1                        ! harmonic
  flog=51
  open(flog,file='flog.model.tomo.interpol')

  ! #################################### choix des ondes
  write(flog,*) ' 1. pour ondes P, 2. pour ondes P et S'
  write(*,*) ' 1. pour ondes P, 2. pour ondes P et S'
  read(*,*) chx
  write(flog,*) chx
  ! ####################################  input
  write(flog,*) ' enter the P-model name '
  read(*,'(a)') name_modelp
  write(flog,'(a)') name_modelp
  write(flog,*) ' enter the forward P-model name '
  read(*,'(a)') name_forwardp
  write(flog,'(a)') name_forwardp
  ! ----------- si on choisit P et S
  ! if(chx.eq.2) then
  write(flog,*) ' enter the S-model name even if not necessary '
  write(*,*) ' enter the S-model name even if not necessary '
  read(*,'(a)') name_models
  write(flog,'(a)') name_models
  write(flog,*) ' enter the forward S-model name '
  write(*,*) ' enter the forward S-model name '
  read(*,'(a)') name_forwards
  write(flog,'(a)') name_forwards
  ! endif
  ! ------------   
  write(flog,*) ' enter (x1,x2,x3) origin of the grid ' 
  write(*,*) ' enter (x1,x2,x3) origin of the grid ' 
  read(*,*) x1inv,x2inv,x3inv
  write(flog,*) x1inv,x2inv,x3inv
  write(flog,*) ' enter the model n1inv,n2inv,n3inv '
  write(*,*) ' enter the model n1inv,n2inv,n3inv '
  read(*,*) n1inv,n2inv,n3inv
  write(flog,*) n1inv,n2inv,n3inv
  write(flog,*) ' enter the model h1inv,h2inv,h3inv '
  write(*,*) ' enter the model h1inv,h2inv,h3inv '
  read(*,*) h1inv,h2inv,h3inv
  write(flog,*)  h1inv,h2inv,h3inv
  write(flog,*) ' enter (x1f,x2f,x3f) origin of the forward grid '
  write(*,*) ' enter (x1f,x2f,x3f) origin of the forward grid '
  read(*,*) x1_fwd,x2_fwd,x3_fwd
  write(flog,*)  x1_fwd,x2_fwd,x3_fwd
  write(flog,*) ' enter the forward model n1,n2,n3 '
  write(*,*) ' enter the forward model n1,n2,n3 '
  read(*,*) n1,n2,n3
  write(flog,*) n1,n2,n3
  write(flog,*) ' enter the forward model h1,h2,h3 h1=h2=h3'
  write(*,*) ' enter the forward model h1,h2,h3 : no test '
  read(*,*) h1,h2,h3
  write(flog,*) h1,h2,h3
  ! #####################################  arrays
  ALLOCATE(modinvp(n1inv,n2inv,n3inv))
  ALLOCATE(modforwp(n1,n2,n3))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     ALLOCATE(modinvs(n1inv,n2inv,n3inv))
     ALLOCATE(modforws(n1,n2,n3))
  endif

  ! #####################################  io
  open(7,file=name_modelp,access='direct',recl=4*n1inv*n2inv*n3inv)
  open(8,file=name_forwardp,access='direct',recl=4*n1*n2*n3)
  !------------ if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     open(9,file=name_models,access='direct',recl=4*n1inv*n2inv*n3inv)
     open(10,file=name_forwards,access='direct',recl=4*n1*n2*n3)
  endif

  ! #####################################  reading the inverse model
  read(7,rec=1) modinvp             ! true grid in P
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     read(9,rec=1) modinvs             ! true grid in S
  endif

  ! ####################################  make the interpolation from inversion to FMM
  write(*,*) ' P wave forward grid computation '
  call subinterpol(modinvp,modforwp,x1inv,x2inv,x3inv,n1inv,n2inv,n3inv, &
       h1inv,h2inv,h3inv,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,iopt)
  write(*,*) ' P min/max velocity values ',vmin,vmax
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     write(*,*) ' S wave forward grid computation '
     call subinterpol(modinvs,modforws,x1inv,x2inv,x3inv,n1inv,n2inv,n3inv, &
          h1inv,h2inv,h3inv,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,iopt)
     write(*,*) ' S min/max velocity values ',vmin,vmax
  endif

  ! #####################################  writing the forward model
  write(8,rec=1) modforwp
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     write(10,rec=1) modforws
  endif
  close(7)
  close(9)
  close(8)
  close(10)
  stop
end program model_interpol






