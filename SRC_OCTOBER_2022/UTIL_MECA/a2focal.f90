!###############################################################################|
! program a2focal    ASCII >> focal                                             |
!                                                                               |
!  alias focal                                                                  |
!  alias focal                                                                  |
!                                                                               |
!   input (ascii file)                                                          |
!                                                                               |
!       focal.asc  (sources)    10*4 = 40 bytes for each record                 |
! BINARY FORMAT     nt,ntp,nts                                                  |
!   record 1+irec   id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase  |
!                                                                               | 
!   output (binary file)                                                        |
!                                                                               |
!       focal                                                                   |
!                                                                               |
!###############################################################################|
program a2focal

implicit none
integer(kind=4) :: id_dat,lut_src,lut_sta,impulse,iphase,nt,ntp,nts
real(kind=4) :: px,py,pz,azimuth,dip

integer(kind=4) :: i4,i5,i6,i7,i8,i9,i10,irec

i4=0;i5=0;i6=0;i7=0;i8=0;i9=0;i10=0
write(*,*) ' file conversion focal into focal.asc '

open(10,file='focal.asc',status='unknown')
open(11,file='focal',access='direct',recl=10*4)

read(10,*) nt,ntp,nts,i4,i5,i6,i7,i8,i9,i10
write(11,rec=1) nt,ntp,nts,i4,i5,i6,i7,i8,i9,i10

write(*,*) '======================================== '
write(*,*)  'number of picks ',nt
write(*,*)  'number of P picks ',ntp
write(*,*)  'number of S picks ',nts
write(*,*) '======================================== '

do irec=1,nt
  read(10,*) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
  write(11,rec=irec+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
enddo

close(unit=10)
close(unit=11)

end program a2focal
