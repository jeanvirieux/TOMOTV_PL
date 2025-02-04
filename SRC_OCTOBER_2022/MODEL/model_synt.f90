!###############################################################################|
! build-up synthetic model and input parameters file                            |
!                                                                               |
! based on a gradient and on a sum of exponential balls                         |
!                                                                               |
! could be used for initial model design as well                                |
!                                                                               |
!===============================================================================|
!                                                                               |
! modified for WINSISM  2015                                                    |
!###############################################################################|
program model_synt
  implicit none

  real(kind=4)             :: xlength,ylength,zlength,dx,dz,dy
  integer(kind=4)          :: nx,ny,nz,ix,iy,iz
  real(kind=4)             :: xorg,yorg,zorg,xtotal,ytotal,ztotal
  integer(kind=4)          :: nxtotal,nytotal,nztotal
  real(kind=4),allocatable,dimension(:,:,:) :: vel3D
  real(kind=4)             :: vel_background,gamx,gamy,gamz
  real(kind=4)             :: xbou,ybou,zbou,xbou_len,ybou_len,zbou_len
  real(kind=4)             :: distx,disty,distz,vel_ano
  integer(kind=4)          :: ibou

  real(kind=4)             :: df,d_ray
  real(kind=4)             :: depth,offsetx,offsety,vel_depth,vel_offsetx,vel_offsety

  write(*,*) ' enter the dimensions of the model (meters) '
  read(*,*) xlength,ylength,zlength
  write(*,*) ' enter the expected sampling in x, y and z (meters) '
  read(*,*) dx,dy,dz
  nx=int(xlength/dx)+1
  ny=int(ylength/dy)+1  
  nz=int(zlength/dz)+1

  nxtotal=nx+2      ! one additional point on both sides
  nytotal=ny+2
  nztotal=nz+2      ! be careful about the positions of sources/receivers
  !=========================
  xlength=(nx-1)*dx ! rounding
  ylength=(ny-1)*dy
  zlength=(nz-1)*dz
  write(*,*) ' we assume an origin at (-dx,-dy,-dz) '
  write(*,*) ' sources and receivers should not be exactly at borders '
  !================= xtotal=xlength+2*dx=(nx-1)*dx+2*dx=(nx+1)*dx=(nxtotal-1)*dx
  xorg=-dx
  yorg=-dy
  zorg=-dz
  xtotal=(nxtotal-1)*dx
  ytotal=(nytotal-1)*dy
  ztotal=(nztotal-1)*dz
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
  ! medium
  !
  !======================================================

  write(*,*) ' enter top velocity '
  read(*,*) vel_background
  write(*,*) ' enter depth and specified velocity '
  read(*,*) depth, vel_depth
  write(*,*) ' enter offsetx and specified velocity at this offset'
  read(*,*) offsetx,vel_offsetx
  write(*,*) ' enter offsety and specified velocity at this offset'
  read(*,*) offsety,vel_offsety
  !
  !  write(*,*) ' enter horizontal and vertical gradients '
  !  read(*,*) gamx,gamz
  gamz= (vel_depth-vel_background)/depth
  gamx= (vel_offsetx-vel_background)/offsetx
  gamy= (vel_offsety-vel_background)/offsety

  write(*,*) ' vertical and horizontal velocity gradients '
  write(*,*) ' vertical one: ',gamz
  write(*,*) ' horizontal one along x: ',gamx
  write(*,*) ' horizontal one along y: ',gamy

  allocate(vel3D(nxtotal,nytotal,nztotal))

  do iz=1,nztotal
     do iy=1,nytotal
        do ix=1,nxtotal
           vel3D(ix,iy,iz)=vel_background+(ix-1)*dx*gamx+(iy-1)*dy*gamy+(iz-1)*dz*gamz
        enddo
     enddo
  enddo

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ adding exponential ball
  ibou=1

  do while (ibou > 0)

     write(*,*) ' enter ball option <0 means stop '
     read(*,*) ibou
     write(*,*) ' enter position of the center of the ball '
     read(*,*) xbou,ybou,zbou
     write(*,*) ' enter horizontal radius and vertical radius '
     read(*,*) xbou_len,ybou_len,zbou_len
     write(*,*) ' enter velocity anomaly wrt background velocity '
     read(*,*) vel_ano

     do iz=1,nztotal
        distz=((iz-1)*dz-zbou)/zbou_len
        distz=distz*distz
        do iy=1,nytotal
           disty=((iy-1)*dy-ybou)/ybou_len
           disty=disty*disty
           do ix=1,nxtotal
              distx=((ix-1)*dx-xbou)/xbou_len
              distx=distx*distx
              vel3D(ix,iy,iz)=vel3D(ix,iy,iz)+vel_ano*exp(-distx-disty-distz)
           enddo
        enddo
     enddo

  enddo

  open(7,file='modelP.ini',access='direct',recl=4*nxtotal*nytotal*nztotal)
  write(7,rec=1) vel3D
  close(7)

  deallocate(vel3D)

  !=========================================== input parameters

  write(*,*) ' discretization of the forward problem h (meters) '
  read(*,*) df
  write(*,*) ' discretization of the ray sampling (meters) '
  read(*,*) d_ray

  open(7,file='modelP.head')

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
  close(7)

  open(7,file='inversion.head')
  write(7,'(a)') '# residual weighting '
  write(7,*) 20.0,40.0
  write(7,'(a)') '# ray weighting '
  write(7,*) 100000.,200000.
  write(7,'(a)') '# activate the smoothing '
  write(7,*) 1
  write(7,'(a)') '# Laplacian option 1 to 3 '
  write(7,*) 3
  write(7,'(a)') '# weight along x '
  write(7,*) 1.
  write(7,'(a)') '# weight along y '
  write(7,*) 1.
  write(7,'(a)') '# weight along z '
  write(7,*) 1.
  write(7,'(a)') '# damping parameter for lsqr '
  write(7,*) 0.5
  write(7,'(a)') '# maximum number of lsqr iterations '
  write(7,*) 500
  write(7,'(a)') '# VP perturbation maximum (m/s)'
  write(7,*) 800.
  write(7,'(a)') '# VS perturbation maximum (m/s)'
  write(7,*) 600.
  write(7,'(a)') '# Horizontal shift maximum (m)'
  write(7,*) 150.
  write(7,'(a)') '# vertical shift maximum (m)'
  write(7,*) 50.
  write(7,'(a)') '# original time shift maximum (s)'
  write(7,*) 0.5
  close(7)

  stop
end program model_synt
