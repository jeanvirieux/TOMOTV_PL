!***********************************************************************
!   moving from a 2D model towards a 3D model
!
! cautious : permutation from z as the fastest index to x 
!-----------------------------------------------------------------------
! 2016 : a dummy version
!***********************************************************************
program model2Dto3D
implicit none

REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modv_3d
REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:), TARGET :: modv_2d
!----------------------------- variables

character(len=1) :: string

integer(kind=4) :: nx,ny,nz,ix,iz,flog
real(kind=4) :: vitesse

flog=51
open(flog,file='flog.2Dto3D')

! ####################################  input
!
!  lecture dans model.head
!
open(9,file='model.head',status='old')
read(9,*) string
read(9,*) string
read(9,*) string
read(9,*) string
read(9,*) string
read(9,*) nx,ny,nz

if(ny /= 4) then
write(*,*) ' ny different from 4 '
stop
endif
write(flog,*) ' nx,ny,nz ',nx,ny,nz
write(*,*) ' nx,ny,nz ',nx,ny,nz

! #####################################  arrays
ALLOCATE(modv_2d(nz,nx))                ! z first
ALLOCATE(modv_3d(nx,ny,nz))             ! x first

! #####################################  io
open(7,file='modelP.2d',access='direct',recl=4*nx*nz)
open(8,file='modelP',access='direct',recl=4*nx*ny*nz)

! #####################################  reading the 2D model
read(7,rec=1) modv_2d             ! true grid

! ##################################### extension to 3D and permutation

do iz=1,nz
do ix=1,nx
vitesse=modv_2d(iz,ix)
modv_3d(ix,1,iz)=vitesse
modv_3d(ix,2,iz)=vitesse
modv_3d(ix,3,iz)=vitesse
modv_3d(ix,4,iz)=vitesse
enddo
enddo

! #####################################  writing the 3D model
write(8,rec=1) modv_3d

close(7)
close(8)
DEALLOCATE(modv_2d)
DEALLOCATE(modv_3d)
stop
endprogram model2Dto3D






