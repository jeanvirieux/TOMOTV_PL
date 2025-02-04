!***********************************************************************
!   moving from a 2D model towards a 3D model
!   adding extra boundary nodes
!   computing the origin point of the grid while keeping the same 
!   coordinate system for stations and shots
!
! cautious : permutation from z as the fastest index to x 
!-----------------------------------------------------------------------
! 2014
!***********************************************************************
program model2Dto3D
implicit none

REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modv_3d
REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:), TARGET :: modv_2d
!----------------------------- variables

character(len=132) :: name_model_2d                 ! filename for model and FMM for P wave
character(len=132) :: name_model_3d                 ! filename for model and FMM for S wave

integer(kind=4) :: nx,ny,nz,ix,iy,iz,flog,ixc,izc
real(kind=4) :: vitesse,scale
real(kind=4) :: dx,dy,dz,vmin,vmax

flog=51
open(flog,file='flog.2Dto3D')

! ####################################  input
write(flog,*) ' enter the 2d model name '
read(*,'(a)') name_model_2d
write(flog,'(a)') name_model_2d
write(flog,*) ' enter the 3d model name '
read(*,'(a)') name_model_3d
write(flog,'(a)') name_model_3d
!
!
!
write(*,*) ' enter nx,nz '
read(*,*) nx,nz
write(*,*) ' enter dx,dz and dy '
read(*,*) dx,dz,dy
ny=4
write(flog,*) ' final nx,ny,nz ',nx+2,ny,nz+2
write(*,*) ' final nx,ny,nz ',nx+2,ny,nz+2
scale=1.
write(*,*) ' enter scale if units have to be changed (1.00 will be neutral) ',scale
read(*,*) scale

! #####################################  arrays
ALLOCATE(modv_2d(nx,nz))                      ! x first ... 
ALLOCATE(modv_3d(nx+2,ny,nz+2))             ! x first ... add two points

! #####################################  io
open(7,file=name_model_2d,access='direct',recl=4*nx*nz)
open(8,file=name_model_3d,access='direct',recl=4*(nx+2)*ny*(nz+2))

! #####################################  reading the 2D model
read(7,rec=1) modv_2d             ! true grid

! ##################################### extension to 3D and permutation

do iz=1,nz+2
izc=iz
if(iz == 1) izc=2
if(iz == nz+2) izc=nz+1

do ix=1,nx+2
ixc=ix
if(ix == 1) ixc=2
if(ix == nx+2) ixc=nx+1

vitesse=modv_2d(ixc-1,izc-1)*scale

modv_3d(ix,1,iz)=vitesse
modv_3d(ix,2,iz)=vitesse
modv_3d(ix,3,iz)=vitesse
modv_3d(ix,4,iz)=vitesse
enddo
enddo

vmin=1.e+29
vmax=-1.e+29

do iz=1,nz+2
do iy=1,4
do ix=1,nx+2
if(modv_3d(ix,iy,iz) > vmax) vmax= modv_3d(ix,iy,iz)
if(modv_3d(ix,iy,iz) < vmin) vmin= modv_3d(ix,iy,iz)
enddo
enddo
enddo

write(*,*) ' vmin, vmax ',vmin, vmax
write(*,*) ' origin ',-dx,-dy,-dz

! #####################################  writing the 3D model
write(8,rec=1) modv_3d

close(7)
close(8)
DEALLOCATE(modv_2d)
DEALLOCATE(modv_3d)
stop
endprogram model2Dto3D






