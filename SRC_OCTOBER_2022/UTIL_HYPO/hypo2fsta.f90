!###############################################################################|
! PROGRAM hypo2fsta     HYPO71 >>> fsta
!	  input
!       file.inp   
!     output (binary file)
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!
!
!###############################################################################|
program hypo2fsta
implicit none

integer(kind=4) :: RefEllipsoid
character(len=1) ::  UTMZone
integer(kind=4) :: ZoneNumber

integer(kind=4) :: ista,nsta,id_sta,i1,i2,i3,i4,i5,i6
real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta
character(len=4) :: name
character(len=132) :: string
character(len=132) :: name_hypo     ! file name for HYPO71 input file
integer(kind=4) :: iempty

open(55,file='check_utm')
call constants    ! set constants for UTM projection
RefEllipsoid = 24 ! 24 is WGS-84. 

write(*,*) ' enter the input HYPO71 file *.inp'
read(*,'(a)') name_hypo
open(11,file=name_hypo,status='old')
open(10,file='fsta',access='direct',recl=7*4) ! we open the stream on file fsta

nsta=0;i1=0;i2=0;i3=0;i4=0;i5=0;i6=0
write(10,rec=1) nsta,i1,i2,i3,i4,i5,i6   ! dummy line upto now

!
! we assume one blank line before stations and one blank line after stations
!

1000 continue
read(11,'(a)') string

! if for detecting empty
iempty=0
call test_empty(string,iempty)
if(iempty == 1) goto 1000 ! if it is not blank I loop to 1000

! first blank line has been detected

ista=0
2000  continue
read(11,'(a132)') string
! is it a blank line if yes jump
iempty=0
call test_empty(string,iempty)
if(iempty == 0) goto 3000  ! if blank jump to 3000

ista=ista+1

! extract positions lat lon ... do utm conversion
call read_string(string,x_sta,y_sta,z_sta,name,RefEllipsoid,UTMZone,Zonenumber) 

dtp_sta=0.
dts_sta=0.
id_sta=ista

write(*,*) ' ************* '
write(*,*) ' station ',id_sta
write(*,*) x_sta,y_sta,z_sta
write(*,'(a)') name
write(*,*) ' ************* '

write(10,rec=1+ista) id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name

goto 2000

3000 continue
nsta=ista
write(10,rec=1) nsta,i1,i2,i3,i4,i5,i6   ! dummy line upto now
write(*,*) ' number of detected stations ', nsta
close(10)
close(11)

stop
end program hypo2fsta

subroutine test_empty(string,iempty)
!
!  use of the function len_trim   iempty=0 means the line is blank
!
character*(*) :: string
integer(kind=4) :: iempty

iempty=1   ! I assume that the line is not blank a priori
if(len_trim(string) == 0 ) then
  iempty=0
endif
return
end subroutine

subroutine read_string(string,x_sta,y_sta,z_sta,name,RefEllipsoid,UTMZone,Zonenumber)

integer(kind=4) :: RefEllipsoid
character(len=1) ::  UTMZone
integer(kind=4) :: ZoneNumber

character*(*) :: string
character*(*) :: name
real(kind=4) :: x_sta,y_sta,z_sta
real(kind=8) :: xlat,xlon,x_sta_d,y_sta_d

integer(kind=4) :: lat,lon,iz_sta
character(len=1) :: conv_lat,conv_lon

ioffset=3
name=string(ioffset:ioffset+3)  ! extraction of the name of the station

read(string,'(6x,i2,f5.2,a1)') lat,xlat,conv_lat          !   we assume convention NORTH
if(conv_lat /= 'N') then
  write(*,*) ' error in the North convention '
  stop 
endif
xlat=lat+xlat/60.   ! conversion in decimal degree
read(string,'(15x,i2,f5.2,a1)') lon,xlon,conv_lon          !   we assume convention EAST
if(conv_lon /= 'E') then
  write(*,*) ' error in the East convention '
  stop
endif
xlon=lon+xlon/60.   ! conversion in decimal degree

call lltoutm(RefEllipsoid,xlat,xlon,y_sta_d,x_sta_d,UTMZone,ZoneNumber)
!         X along Easting
!         Y along Northing

write(*,*) ' #################### '
write(*,*) 'xlat,xlon',xlat,xlon
write(*,*) 'x_sta,y_sta',x_sta_d,y_sta_d
write(*,*) ' #################### '

write(55,*) UTMZone,ZoneNumber


x_sta=sngl(x_sta_d)
y_sta=sngl(y_sta_d)

read(string,'(23x,i4)') iz_sta   ! altitude as an integer

z_sta=-iz_sta   ! depth in meter

return
end subroutine
