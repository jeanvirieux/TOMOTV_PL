!***********************************************************************
!   Making a 3D model into a 2.5D model
!
!=======================================================================
! valid only if ny=4
!-----------------------------------------------------------------------
! 2014
!***********************************************************************
program model3Dto2_5D
implicit none

REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modv_3d
REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:), TARGET :: modv_2d
REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modv_2_5d

!----------------------------- variables
real(kind=4) :: vitesse
integer(kind=4) :: nx,ny,nz,ix,iz,flog
character(len=1) :: carac

flog=51
open(flog,file='flog.model2D-3D')
!
!
!
open(49,file='model.head',status='old',err=1000)
read(49,'(a)') carac
read(49,'(a)') carac
read(49,'(a)') carac
read(49,'(a)') carac
read(49,'(a)') carac
read(49,*) nx,ny,nz
if( ny /= 4) then
  write(flog,*) ' error in the conversion from 3D to 2.5D ny=4 ',ny
  write(*,*) ' error in the conversion from 3D to 2.5D ny=4 ',ny
  stop
endif
write(flog,*) ' dimension of the 3D model nx,ny,nz ',nx,ny,nz

! #####################################  arrays
ALLOCATE(modv_3d(nx,ny,nz))
ALLOCATE(modv_2_5d(nx,ny,nz))
ALLOCATE(modv_2d(nx,nz))

! #####################################  reading the 3D model
open(7,file='modelP',access='direct',recl=4*nx*ny*nz)
read(7,rec=1) modv_3d             ! true grid
close(7)

! ##################################### making the 2.5d medium 
do iz=1,nz
do ix=1,nx
vitesse=0.5*(modv_3d(ix,2,iz)+modv_3d(ix,3,iz))
modv_2_5d(ix,1,iz)=vitesse
modv_2_5d(ix,2,iz)=vitesse
modv_2_5d(ix,3,iz)=vitesse
modv_2_5d(ix,4,iz)=vitesse
modv_2d(ix,iz)=vitesse
enddo
enddo

! #####################################  writing the 3D model
open(8,file='modelP',access='direct',recl=4*nx*ny*nz)
write(8,rec=1) modv_2_5d
close(8)

!======================== for plotting with ximage
open(9,file='modelP.2d',access='direct',recl=4*nx*nz)
write(9,rec=1) ((modv_2d(ix,iz),iz=1,nz),ix=1,nx)
close(9)

! ##################################### making the 2.5d medium 
! ##################################### ascii format
open(10,file='modelP.asc')
do iz=2,nz-1
do ix=2,nx-1
write(10,*) modv_2d(ix,iz)
enddo
enddo

DEALLOCATE(modv_3d)
DEALLOCATE(modv_2_5d)
DEALLOCATE(modv_2d)

write(flog,*) ' end of the 2.5 conversion and ascii output '
stop

1000 continue
write(flog,*) ' missing the file model.head'
stop
end program model3Dto2_5D
