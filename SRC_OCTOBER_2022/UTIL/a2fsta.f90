!###############################################################################|
! PROGRAM a2fsta     ASCII >>> fsta
!
!       input (ascii file)
!       fsta.asc                                                                |
! ASCII  FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!                                                                               |
!
!     output (binary file)
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!
!
!###############################################################################|
program a2fsta
implicit none

integer(kind=4) :: ista,nsta,id_sta,i1,i2,i3,i4,i5,i6
real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta
character(len=4) :: name
character(len=162) :: string

open(11,file='fsta.asc',status='old')
open(10,file='fsta',access='direct',recl=7*4) ! we open the stream on file fsta

read(11,*) nsta,i1,i2,i3,i4,i5,i6
write(10,rec=1) nsta,i1,i2,i3,i4,i5,i6

write(*,*) ' conversion of ascii file fsta.asc into binary file fsta'
write(*,*) ' number of stations ',nsta

do ista=1,nsta
read(11,'(A162)') string
write(*,*) string
read(string,*) id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta
! il faut prendre les 4 derniers caracteres pour le nom
name(1:4)=string(len_trim(string)-3:len_trim(string))
write(10,rec=1+ista) id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name

enddo

close(10)
close(11)

stop
end program a2fsta
