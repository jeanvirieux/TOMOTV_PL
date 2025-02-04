!***********************************************************************
! module : model.cross     XXXXXXX  PASSAGE EN KM/S  XXXXXXX
! ==============
! purpose : interpolation the inversion model into a new model
!           with different spatial discretization
!
! ==== and extraction of horizontal sections put into the XYZ
!           with correct file names for Anne's script
!
! the log file name is "flog.model.tomo.cross" 
!
!-----------------------------------------------------------------------
!
!***********************************************************************
program model_cross
use s_interpol_mod
implicit none

REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinv,modnew
REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:), TARGET :: modcross

!----------------------------- variables
real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
real(kind=4) h1,h2,h3,x1_fwd,x2_fwd,x3_fwd       ! spatial step for FMM and origin position
integer(kind=4) n1inv,n2inv,n3inv,n1,n2,n3       ! dimensions for inversion grid and for FMM grid
integer(kind=4) flog                             ! output flux

real(kind=4) :: vmin,vmax
integer(kind=4) :: ierr,ilen,idepth,jdepth,ndepth

real(kind=4) :: depth_min,depth_max,depth_dz,depth

character(len=2) :: carac

integer(kind=4) :: iopt=0      ! harmonic
real(kind=4) :: scale=1000.          ! 1000.   m/s >>> km/s

flog=51
open(flog,file='flog.model.cross')
write(*,*) ' m or km scale '
read(*,*) scale
if(abs(scale-1000.) < 0.1) iopt=1    ! harmonic si conversion m/s into km/s
if(abs(scale+1000.) < 0.1) then
   iopt=0
   scale=-scale
endif
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
write(flog,*) ' enter the new model h1,h2,h3 '
write(*,*) ' enter the forward model h1,h2,h3 : no test '
read(*,*) h1,h2,h3
write(flog,*) h1,h2,h3
n1=int(((x1inv+(n1inv-1)*h1inv)-x1_fwd)/h1)+1 
n2=int(((x2inv+(n2inv-1)*h2inv)-x2_fwd)/h2)+1 
n3=int(((x3inv+(n3inv-1)*h3inv)-x3_fwd)/h3)+1 
write(flog,*) ' (n1,n2,n3) extended grid ',n1,n2,n3

! #####################################  arrays
ALLOCATE(modinv(n1inv,n2inv,n3inv))
ALLOCATE(modnew(n1,n2,n3))

! #####################################  io
open(7,file='model.old',access='direct',recl=4*n1inv*n2inv*n3inv)
open(8,file='model.new',access='direct',recl=4*n1*n2*n3)

! #####################################  reading the inverse model
read(7,rec=1) modinv             ! true grid in P

modinv(:,:,:)=modinv(:,:,:)/scale    !!!! m/s >>>> km/s

! ####################################  make the interpolation from inversion to FMM
write(*,*) ' P wave forward grid computation '
call subinterpol(modinv,modnew,x1inv,x2inv,x3inv,n1inv,n2inv,n3inv, &
                 h1inv,h2inv,h3inv,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,iopt)
! #####################################  writing the forward model
write(8,rec=1) modnew
close(7)
close(8)

deallocate(modinv)

allocate(modcross(n1,n2))

write(*,*) ' enter depth min,  depth max and depth step '
read(*,*) depth_min,depth_max,depth_dz

ndepth=1+int((depth_max-depth_min)/depth_dz)

do idepth=1,ndepth

   depth=depth_min+float(idepth-1)*depth_dz
   jdepth=depth/1000.
   write(carac,'(I2.2)') jdepth

do ilen=1,len(carac)
   if(carac(ilen:ilen) == ' ') carac(ilen:ilen)='0'
enddo   

call subcross_z(depth,modnew,modcross,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)

if(ierr == 0) then
open(8,file='XYZ/vp_lay'//carac//'.bin',access='direct',recl=4*n1*n2)
write(8,rec=1) modcross
close(8)
endif

enddo

write(*,*) ' Extraction of horizontal sections has ended',n1,n2,vmin,vmax
write(*,*) ' [depth min,depth max] with depth step ',depth_min,depth_max,depth_dz

deallocate(modcross)
deallocate(modnew)

contains

  subroutine subcross_z(depth,v_model,v_cross,x1orig,x2orig,x3orig,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)
    implicit none
    integer(kind=4) :: n1,n2,n3
    real(kind=4) :: v_model(n1,n2,n3)
    real(kind=4) :: v_cross(n1,n2)
    real(kind=4) :: h1,h2,h3
    real(kind=4) :: x1orig,x2orig,x3orig
    !
    real(kind=4) :: depth,zz
    integer(kind=4) :: idepth

    integer(kind=4) :: i,j,ierr
    real(kind=4) :: vmin,vmax

    ierr=999
    vmin=1.e+29; vmax=-1.e29

    idepth=1+int((depth-x3orig)/h3)
    if(idepth > n3) idepth=n3-1
    zz=(depth-x3orig)/h3-float(idepth-1)
    if(zz < 0.) then
       write(*,*) ' depth',depth,'zz',zz,'should positive'
       return
    endif
    if(zz > 1.) then
       write(*,*) ' depth',depth,'zz',zz,'should < 1.'
       return
    endif
    do j=1,n2
       do i=1,n1
          v_cross(i,j)=v_model(i,j,idepth)+(v_model(i,j,idepth+1)-v_model(i,j,idepth))*zz
          if(v_cross(i,j) > vmax) vmax=v_cross(i,j)
          if(v_cross(i,j) < vmin) vmin=v_cross(i,j)
       enddo
    enddo
    ierr=0
    
  end subroutine subcross_z
  
end program model_cross






