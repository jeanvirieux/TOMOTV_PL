!***********************************************************************
! module : model.coupe     XXXXXXX  PASSAGE EN KM/S  XXXXXXX
! ==============
! purpose : interpolation the inversion model into a new model
!           with different spatial discretization
!
! ==== and extraction of horizontal sections
!
! the log file name is "flog.model.tomo.coupe" 
!
!-----------------------------------------------------------------------
!
!***********************************************************************
program model_coupe
use s_interpol_mod
implicit none

REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinv,modnew
REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:), TARGET :: modcoupe

!----------------------------- variables
real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
real(kind=4) h1,h2,h3,x1_fwd,x2_fwd,x3_fwd       ! spatial step for FMM and origin position
integer(kind=4) n1inv,n2inv,n3inv,n1,n2,n3       ! dimensions for inversion grid and for FMM grid
integer(kind=4) flog                             ! output flux

real(kind=4) :: vmin,vmax,depth,ysec,xsec
integer(kind=4) :: ierr,ilen

real(kind=4) :: scale=1000.

character(len=5) :: carac

integer(kind=4) :: iopt=0

flog=51
open(flog,file='flog.model.coupe')
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
write(*,*) ' (n1,n2,n3) extended grid ',n1,n2,n3

! #####################################  arrays
ALLOCATE(modinv(n1inv,n2inv,n3inv))
ALLOCATE(modnew(n1,n2,n3))

! #####################################  io
open(7,file='model.old',access='direct',recl=4*n1inv*n2inv*n3inv)
open(8,file='model.new',access='direct',recl=4*n1*n2*n3)

! #####################################  reading the inverse model
read(7,rec=1) modinv             ! true grid in P

modinv(:,:,:)=modinv(:,:,:)/scale    !!!! m/s >>>> km/s

write(*,*) ' min/max modinv ',minval(modinv),maxval(modinv)

! ####################################  make the interpolation from inversion 
write(*,*) ' P wave forward grid computation '
call subinterpol(modinv,modnew,x1inv,x2inv,x3inv,n1inv,n2inv,n3inv, &
                 h1inv,h2inv,h3inv,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,iopt)
! #####################################  writing the forward model
write(8,rec=1) modnew
close(7)
close(8)

write(*,*) ' vmin,vmax ',vmin,vmax

write(*,*) ' min/max modnew ',minval(modnew),maxval(modnew)

deallocate(modinv)



allocate(modcoupe(n1,n2))

1000 continue

write(*,*) ' enter the depth of the section (>100000 = stop) '
read(*,*) depth
write(carac,'(F5.0)') depth/1000.

do ilen=1,len(carac)
   if(carac(ilen:ilen) == ' ') carac(ilen:ilen)='0'
enddo   

if(depth > 100000.) goto 2000

call subcoupe_z(depth,modnew,modcoupe,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)

if(ierr == 0) then
open(8,file='hori_section.'//carac//'bin',access='direct',recl=4*n1*n2)
write(8,rec=1) modcoupe
close(8)
endif

write(*,*) 'depth',depth/1000.,' min/max modcoupe ',minval(modcoupe),maxval(modcoupe)

goto 1000
2000 continue

write(*,*) ' end of section extraction ',n1,n2,vmin,vmax
deallocate(modcoupe)



allocate(modcoupe(n1,n3))
4000 continue
write(*,*) ' enter the constant y position of the vertical section (>100000 = stop) '
read(*,*) ysec
write(carac,'(F5.0)') ysec/1000.
do ilen=1,len(carac)
   if(carac(ilen:ilen) == ' ') carac(ilen:ilen)='0'
enddo   

if(ysec > 1000000.) goto 3000

call subcoupe_x(ysec,modnew,modcoupe,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)

if(ierr == 0) then
open(8,file='vert_x_section.'//carac//'bin',access='direct',recl=4*n1*n3)
write(8,rec=1) modcoupe
close(8)
endif

write(*,*) 'ysec',ysec/1000.,' min/max coupe ',minval(modcoupe),maxval(modcoupe)

goto 4000

3000 continue
write(*,*) ' end of vertical section extraction y=cte ',n1,n3,vmin,vmax

deallocate(modcoupe)

allocate(modcoupe(n2,n3))
6000 continue
write(*,*) ' enter the constant x position of the vertical section (>100000 = stop) '
read(*,*) xsec
write(carac,'(F5.0)') xsec/1000.
do ilen=1,len(carac)
   if(carac(ilen:ilen) == ' ') carac(ilen:ilen)='0'
enddo   

if(xsec > 1000000.) goto 5000

call subcoupe_y(xsec,modnew,modcoupe,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)

if(ierr == 0) then
open(8,file='vert_y_section.'//carac//'bin',access='direct',recl=4*n2*n3)
write(8,rec=1) modcoupe
close(8)
endif

write(*,*) 'xsec',xsec/1000.,' min/max coupe ',minval(modcoupe),maxval(modcoupe)

goto 6000
5000 continue
write(*,*) ' end of vertical section extraction x=cte ',n2,n3,vmin,vmax

deallocate(modcoupe)



deallocate(modnew)

contains

  subroutine subcoupe_z(depth,v_model,v_coupe,x1orig,x2orig,x3orig,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)
    implicit none
    integer(kind=4) :: n1,n2,n3
    real(kind=4) :: v_model(n1,n2,n3)
    real(kind=4) :: v_coupe(n1,n2)
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

    !    write(*,*) 'depth,zorig,hz,idepth,zz',depth,x3orig,h3,idepth,zz
    
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
          v_coupe(i,j)=v_model(i,j,idepth)+(v_model(i,j,idepth+1)-v_model(i,j,idepth))*zz
          if(v_coupe(i,j) > vmax) vmax=v_coupe(i,j)
          if(v_coupe(i,j) < vmin) vmin=v_coupe(i,j)
       enddo
    enddo
    ierr=0
    
  end subroutine subcoupe_z
  
  subroutine subcoupe_x(ysec,v_model,v_coupe,x1orig,x2orig,x3orig,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)
    implicit none
    integer(kind=4) :: n1,n2,n3
    real(kind=4) :: v_model(n1,n2,n3)
    real(kind=4) :: v_coupe(n1,n3)
    real(kind=4) :: h1,h2,h3
    real(kind=4) :: x1orig,x2orig,x3orig
    !
    real(kind=4) :: ysec,yy
    integer(kind=4) :: iy

    integer(kind=4) :: i,k,ierr
    real(kind=4) :: vmin,vmax

    ierr=999
    vmin=1.e+29; vmax=-1.e29

    iy=1+int((ysec-x2orig)/h2)
    if(iy > n2) iy=n2-1
    yy=(ysec-x2orig)/h2-float(iy-1)
    write(*,*) ' yy ',yy,iy
    if(yy < 0.) then
       write(*,*) ' ysec',ysec,'yy',yy,'should positive'
       return
    endif
    if(yy > 1.) then
       write(*,*) ' ysec',ysec,'yy',yy,'should < 1.'
       return
    endif
    do k=1,n3
       do i=1,n1
          v_coupe(i,k)=v_model(i,iy,k)+(v_model(i,iy+1,k)-v_model(i,iy,k))*yy
          if(v_coupe(i,k) > vmax) vmax=v_coupe(i,k)
          if(v_coupe(i,k) < vmin) vmin=v_coupe(i,k)
       enddo
    enddo
    ierr=0
    
  end subroutine subcoupe_x
  
  subroutine subcoupe_y(xsec,v_model,v_coupe,x1orig,x2orig,x3orig,n1,n2,n3,h1,h2,h3,vmin,vmax,ierr)
    implicit none
    integer(kind=4) :: n1,n2,n3
    real(kind=4) :: v_model(n1,n2,n3)
    real(kind=4) :: v_coupe(n2,n3)
    real(kind=4) :: h1,h2,h3
    real(kind=4) :: x1orig,x2orig,x3orig
    !
    real(kind=4) :: xsec,xx
    integer(kind=4) :: ix

    integer(kind=4) :: j,k,ierr
    real(kind=4) :: vmin,vmax

    ierr=999
    vmin=1.e+29; vmax=-1.e29

    ix=1+int((xsec-x1orig)/h1)
    if(ix > n1) ix=n1-1
    xx=(xsec-x1orig)/h1-float(ix-1)
    write(*,*) ' xx ',xx,ix
    if(xx < 0.) then
       write(*,*) ' xsec',xsec,'xx',xx,'should positive'
       return
    endif
    if(xx > 1.) then
       write(*,*) ' xsec',xsec,'xx',xx,'should < 1.'
       return
    endif
    do k=1,n3
       do j=1,n2
          v_coupe(j,k)=v_model(ix,j,k)+(v_model(ix+1,j,k)-v_model(ix,j,k))*xx
          if(v_coupe(j,k) > vmax) vmax=v_coupe(j,k)
          if(v_coupe(j,k) < vmin) vmin=v_coupe(j,k)
       enddo
    enddo
    ierr=0
    
  end subroutine subcoupe_y
  
end program model_coupe






