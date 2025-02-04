!###############################################################################|
! PROGRAM sta2fsta    conversion from a station network description to fsta     |
!                                                                               |
!       input (ascii file)                                                      |
!       sta.txt    format is case dependent and should vary                     |
!       ASCII FORMAT TO BE ADAPTED    <<<<<<<<<<<<                              |
!                                                                               |
!     output (binary file)                                                      |
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!                                                                               |
!     output (ascii file)
!       network.txt
!                                                                               |
!###############################################################################|
program sta2fsta
implicit none

integer(kind=4) :: ista,nsta,id_sta,i1,i2,i3,i4,i5,i6
real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta
real(kind=4) :: x_orig,y_orig
character(len=4) :: name
character(len=8) :: nom_station
character(len=162) :: string

write(*,*) ' enter origin point (x0,y0) '
read(*,*) x_orig,y_orig

!===================================== find the number of stations (one per line)
open(11,file='sta.txt',status='old')
nsta=1
100 continue
read(11,*,end=200) name
nsta=nsta+1
goto 100
200 continue
nsta=nsta-1
close(11)

write(*,*) ' number of stations nsta:',nsta

open(11,file='sta.txt',status='old')
open(10,file='fsta',access='direct',recl=7*4) ! we open the stream on file fsta
open(12,file='network.txt',status='unknown')
 i1=0;i2=0;i3=0;i4=0;i5=0;i6=0
write(10,rec=1) nsta,i1,i2,i3,i4,i5,i6

write(*,*) ' conversion of ascii file sta.txt into binary file fsta'
write(*,*) ' writing a file network.txt where is the relation between id_sta and station name'

!=================== set statics to zero
dtp_sta=0.
dts_sta=0.

do ista=1,nsta
read(11,'(A162)') string

!============================================================ depend on the file sta.txt
!============================================================ should be tuned for different formats ...
nom_station='        '
if(string(9:9) == ' ') then
  read(string,'(9X,F6.0,1X,F7.0,1X,F4.0)') x_sta,y_sta,z_sta
  read(string,'(A8)') nom_station
elseif(string(8:8) == ' ') then
  read(string,'(8X,F6.0,1X,F7.0,1X,F4.0)') x_sta,y_sta,z_sta
  read(string,'(A7)') nom_station
elseif(string(7:7) == ' ') then
  read(string,'(7X,F6.0,1X,F7.0,1X,F4.0)') x_sta,y_sta,z_sta
  read(string,'(A6)') nom_station
elseif(string(6:6) == ' ') then
  read(string,'(6X,F6.0,1X,F7.0,1X,F4.0)') x_sta,y_sta,z_sta
  read(string,'(A5)') nom_station
else
 write (*,*) ' unknown format '
 stop
endif

write(name,'(I4.4)') ista   ! setup this name for tomoTV

!============================================================

x_sta=x_sta-x_orig
y_sta=y_sta-y_orig
z_sta=-z_sta                ! altitude is given

write(10,rec=1+ista) ista,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
write(12,'(A8,1x,I6,1x,A4)') nom_station,ista,name
!write(12,*) nom_station,ista,name

enddo

close(10)
close(11)
close(12)

stop
end program sta2fsta
