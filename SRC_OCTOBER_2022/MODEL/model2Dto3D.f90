!***********************************************************************
!   moving from a 2D model towards a 3D model
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

integer(kind=4) :: nx,ny,nz,ix,iz,flog
real(kind=4) :: vitesse,scale

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
ny=4
write(flog,*) ' nx,ny,nz ',nx,ny,nz
write(*,*) ' nx,ny,nz ',nx,ny,nz
scale=1.
write(*,*) ' enter scale if units have to be changed ',scale
read(*,*) scale

! #####################################  arrays
ALLOCATE(modv_2d(nz,nx))                ! z first
ALLOCATE(modv_3d(nx,ny,nz))             ! x first

! #####################################  io
open(7,file=name_model_2d,access='direct',recl=4*nx*nz)
open(8,file=name_model_3d,access='direct',recl=4*nx*ny*nz)

! #####################################  reading the 2D model
read(7,rec=1) modv_2d             ! true grid

! ##################################### extension to 3D and permutation

do iz=1,nz
do ix=1,nx
vitesse=modv_2d(iz,ix)*scale
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






