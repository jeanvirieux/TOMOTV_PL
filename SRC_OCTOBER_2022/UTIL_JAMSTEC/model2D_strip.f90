!***********************************************************************
!   read output model3D from tomoTV
!   reduce it to a 2D model and remove boundaries
! cautious : permutation from z as the fastest index to x 
!-----------------------------------------------------------------------
! 2016
!***********************************************************************
program model3Dto2D_strip
implicit none

REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modv_3d
REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:), TARGET :: modv_2d
!----------------------------- variables

character(len=132) :: name_model_2d                 ! filename for model and FMM for P wave
character(len=132) :: name_model_3d                 ! filename for model and FMM for S wave

integer(kind=4) :: nx,ny,nz,ix,iz,flog,ixc,izc
real(kind=4) :: scale
real(kind=4) :: dx,dy,dz,vmin,vmax

flog=51
open(flog,file='flog.2D_strip')

! ####################################  input
write(flog,*) ' enter the 3d model name '
read(*,'(a)') name_model_3d
write(flog,'(a)') name_model_3d
write(flog,*) ' enter the 2d model name '
read(*,'(a)') name_model_2d
write(flog,'(a)') name_model_2d
!
!
!
write(*,*) ' enter nx,ny,nz '
read(*,*) nx,ny,nz
write(*,*) ' enter dx,dy,dz '
read(*,*) dx,dy,dz
ny=4
write(flog,*) ' strip 2D nx,nz ',nx-2,nz-2
write(*,*) ' strip 2D nx,nz ',nx-2,nz-2
scale=1.
write(*,*) ' enter scale if units have to be changed (1.00 will be neutral) ',scale
read(*,*) scale

! #####################################  arrays
ALLOCATE(modv_2d(nx-2,nz-2))                 ! x first ... 
ALLOCATE(modv_3d(nx,ny,nz))                  ! x first ... add two points

! #####################################  io
open(7,file=name_model_3d,access='direct',recl=4*nx*ny*nz)
open(8,file=name_model_2d,access='direct',recl=4*(nx-2)*(nz-2))

! #####################################  reading the 3D model
read(7,rec=1) modv_3d             ! true grid

! ##################################### reduction to 2D without boundaries

do iz=1,nz-2
izc=iz+1
do ix=1,nx-2
ixc=ix+1
modv_2d(ix,iz)=0.5*scale*(modv_3d(ixc,2,izc)+modv_3d(ixc,3,izc))
enddo
enddo

vmin=1.e+29
vmax=-1.e+29
do iz=1,nz-2
do ix=1,nx-2
if(modv_2d(ix,iz) > vmax) vmax= modv_2d(ix,iz)
if(modv_2d(ix,iz) < vmin) vmin= modv_2d(ix,iz)
enddo
enddo

write(*,*) ' vmin, vmax ',vmin, vmax

! #####################################  writing the 3D model
write(8,rec=1) modv_2d

close(7)
close(8)
DEALLOCATE(modv_2d)
DEALLOCATE(modv_3d)
stop
end program model3Dto2D_strip






