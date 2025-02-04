!###############################################################################|
! PROGRAM fsta2a     fsta >>> ASCII
!
!       input (binary file)
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!                                                                               |
!
!     output (ascii file)
!       fsta.asc
!
!
!###############################################################################|
program fsta2a
implicit none

integer(kind=4) :: ista,nsta,id_sta,i1,i2,i3,i4,i5,i6
real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta
character(len=4) :: name

real(kind=4) :: xmin,xmax,ymin,ymax,zmin,zmax

xmin=1.e29
xmax=-1.e29
ymin=1.e29
ymax=-1.e29
zmin=1.e29
zmax=-1.e29

open(10,file='fsta',access='direct',recl=7*4) ! we open the stream on file fsta
open(11,file='fsta.asc',status='unknown')

read(10,rec=1) nsta,i1,i2,i3,i4,i5,i6
write(11,'(7i8)') nsta,i1,i2,i3,i4,i5,i6

write(*,*) ' conversion of fsta into ascii file fsta.asc '
write(*,*) ' number of stations ',nsta

do ista=1,nsta
  read(10,rec=1+ista) id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
  write(11,'(I8,5F25.6,2x,A4)') id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
write(*,'(I8,5F25.6,2x,A4)') id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
  if(xmin >= x_sta) xmin=x_sta
  if(ymin >= y_sta) ymin=y_sta
  if(zmin >= z_sta) zmin=z_sta
  if(xmax <= x_sta) xmax=x_sta
  if(ymax <= y_sta) ymax=y_sta
  if(zmax <= z_sta) zmax=z_sta
enddo

write(*,*) ' extension of the box for stations '
write(*,*) ' xmin,xmax ',xmin,xmax
write(*,*) ' ymin,ymax ',ymin,ymax
write(*,*) ' zmin,zmax ',zmin,zmax
write(*,*) ' ================================ '

close(10)
close(11)

stop
end program fsta2a
