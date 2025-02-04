!###############################################################################|
! build-up synthetic model and input parameters file                            |
!                                                                               |
! based on the initial model of Potin (PhD thesis p80)                          |
!                                                                               |
! could be used for initial model design as well                                |
!                                                                               |
!===============================================================================|
!                                                                               |
!###############################################################################|
program model_layer
  implicit none

  real(kind=4)             :: xlength,ylength,zlength,dx,dz,dy
  integer(kind=4)          :: nx,ny,nz,ix,iy,iz
  real(kind=4)             :: xorg,yorg,zorg,xtotal,ytotal,ztotal
  integer(kind=4)          :: nxtotal,nytotal,nztotal
  real(kind=4),allocatable,dimension(:,:,:) :: vel3D

  real(kind=4)             :: df,d_ray
  real(kind=4)             :: vel_depth


  real(kind=4),allocatable,dimension(:) :: prof,speed
  integer(kind=4) :: nprof,iprof
  real(kind=4) :: dprof

  real(kind=4) :: zdepth

  open(77,file='Potin_layer.txt',status='old')
  read(77,*) nprof
  allocate(prof(nprof),speed(nprof))
  do iprof=1,nprof
     read(77,*) prof(iprof),speed(iprof)
  enddo
  close(77)
  dprof=(prof(nprof)-prof(1))/float(nprof-1)

  write(*,*) ' end of reading the layer structure ',dprof
  
  write(*,*) ' enter the dimensions of the model (meters) '
  read(*,*) xlength,ylength,zlength
  write(*,*) ' enter the expected sampling in x, y and z (meters) '
  read(*,*) dx,dy,dz
  write(*,*) ' enter xorg,yorg,zorg '
  read(*,*) xorg,yorg,zorg
  nx=int(xlength/dx)+1
  ny=int(ylength/dy)+1  
  nz=int(zlength/dz)+1

  nxtotal=nx! +2      ! one additional point on both sides
  nytotal=ny! +2      ! one additional point on both sides
  nztotal=nz! +2      ! be careful about the positions of sources/receivers
  !=========================
  xlength=(nx-1)*dx ! rounding
  ylength=(ny-1)*dy
  zlength=(nz-1)*dz
  write(*,*) ' we assume an origin at (-dx,-dy,-dz) '
  write(*,*) ' sources and receivers should not be exactly at borders '
  
 
  xtotal=xorg+float(nxtotal-1)*dx
  ytotal=yorg+float(nytotal-1)*dy
  ztotal=zorg+float(nztotal-1)*dz
  write(*,*) ' computationnal dimensions '
  write(*,*) ' origin:',xorg,yorg,zorg
  write(*,*) ' last point:',xtotal,ytotal,ztotal
  write(*,*) ' active dimensions '
  write(*,*) ' origin:',0.,0.,0.
  write(*,*) ' last point:',xlength,ylength,zlength
  write(*,*) ' number of nodes in each dimension ',nxtotal,nytotal,nztotal
  write(*,*) ' number of internal nodes ',nx,ny,nz

  !======================================================
  !
  ! velocity from the layer structure
  !
  !======================================================

  allocate(vel3D(nxtotal,nytotal,nztotal))
  
  vel_depth=0.
  
  do iz=1,nztotal
     zdepth=zorg+float(iz-1)*dz  ! temporary for getting the iprof index
     iprof=1+int(zdepth/dprof)
     if(iprof<1) iprof=1
     if(iprof > nprof-1) iprof=nprof-1
     vel_depth=speed(iprof) + (speed(iprof+1)-speed(iprof))*(zdepth-prof(iprof))/(prof(iprof+1)-prof(iprof))
!     if(zdepth <= prof(1)) vel_depth=speed(1)
!     if(zdepth >= prof(nprof)) vel_depth=speed(nprof)
!     write(99,*) zdepth,vel_depth
     do iy=1,nytotal
        do ix=1,nxtotal
           vel3D(ix,iy,iz)=vel_depth
        enddo
     enddo
  enddo

  open(7,file='model.ini',access='direct',recl=4*nxtotal*nytotal*nztotal)
  write(7,rec=1) vel3D
  close(7)

  deallocate(vel3D)

  !=========================================== input parameters

  write(*,*) ' discretization of the forward problem h (meters) '
  read(*,*) df
  write(*,*) ' discretization of the ray sampling (meters) '
  read(*,*) d_ray

  open(7,file='model.head') 
  write(7,'(a)') '# P (1) or P+1 (2) to be defined'
  write(7,*) 2
  write(7,'(a)') '# origin of the grid '
  write(7,*) xorg,yorg,zorg
  write(7,'(a)') '# dimensions nx,ny,nz '
  write(7,*) nxtotal,nytotal,nztotal
  write(7,'(a)') '# samplings along x,y,z '
  write(7,*) dx,dy,dz
  write(7,'(a)') '# forward problem sampling (cube)'
  write(7,*) df
  write(7,'(a)') '# ray sampling'
  write(7,*) d_ray
  write(7,'(a)') ' starting node '
  write(7,*) 1,1,1
  write(7,'(a)') ' ending node '
  write(7,*) nxtotal,nytotal,nztotal
  close(7)

  open(7,file='inversion.head')
  write(7,'(a)') '# residual weighting '
  write(7,*)  3.0000000, 4.0000000  
  write(7,'(a)') '# ray weighting '
  write(7,*) 10000000.000, 20000000.000  
  write(7,'(a)') '# activate the smoothing (0=KO,1=OK)'
  write(7,*) 1
  write(7,'(a)') '# Laplacian option 1 to 3 '
  write(7,*) 2
  write(7,'(a)') '# weight along x '
  write(7,*) 5.
  write(7,'(a)') '# weight along y '
  write(7,*) 5.
  write(7,'(a)') '# weight along z '
  write(7,*) 5.
  write(7,'(a)') '# trigger reading for residues weights (0 set weight to 1 and 1 read file fwei) and others'
  write(7,*) 1,0
  write(7,'(a)') '# ratio for P wave derivative, ratio for S wave'
  write(7,*) 1.,1. 
  write(7,'(a)') '# ratio for hypocenter, ratio for origin time'
  write(7,*) 1.,1.
  write(7,'(a)') '# damping parameter for lsqr '
  write(7,*) 0.05
  write(7,'(a)') '# maximum number of lsqr iterations '
  write(7,*) 500
  write(7,'(a)') '# VP perturbation maximum (m/s)'
  write(7,*) 800.
  write(7,'(a)') '# VS perturbation maximum (m/s)'
  write(7,*) 600.
  write(7,'(a)') '# Horizontal shift maximum (m)'
  write(7,*) 1500.
  write(7,'(a)') '# vertical shift maximum (m)'
  write(7,*) 500.
  write(7,'(a)') '# original time shift maximum (s)'
  write(7,*) 1.5
  close(7)

!!! creation du fichier smoothing.head active par le fichier regular.head
!!! qui active le smoothing et/ou le denoising ...  
  open(7,file='smoothing.head')
  write(7,*)  6,nxtotal/3,nytotal/3,nztotal/3
  close(7)
  
  stop
end program model_layer
