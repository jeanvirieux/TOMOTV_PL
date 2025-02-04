!###############################################################################|
! program fsphe2a    fsphe >> ASCII                                             |
!                                                                               |
!  alias fsphe                                                                  |
!  alias focal                                                                  |
!                                                                               |
!   input (binary file)                                                         |
!                                                                               |
!       fsphe      (sources)    10*4 = 40 bytes for each record                 |
! BINARY FORMAT     nt,ntp,nts                                                  |
!   record 1+irec   id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase  |
!                                                                               | 
!   output (ascii file)                                                         |
!                                                                               |
!       fsphe.asc                                                               |
!                                                                               |
!###############################################################################|
program fsphe2a

implicit none
integer(kind=4) :: id_dat,lut_src,lut_sta,impulse,iphase,nt,ntp,nts
real(kind=4) :: px,py,pz,azimuth,dip

integer(kind=4) :: i4,i5,i6,i7,i8,i9,i10,irec

i4=0;i5=0;i6=0;i7=0;i8=0;i9=0;i10=0
write(*,*) ' file conversion fsphe into fsphe.asc '

open(10,file='fsphe',access='direct',recl=10*4)
open(11,file='fsphe.asc',status='unknown')

read(10,rec=1) nt,ntp,nts,i4,i5,i6,i7,i8,i9,i10
write(11,*) nt,ntp,nts,i4,i5,i6,i7,i8,i9,i10

write(*,*) '======================================== '
write(*,*)  'number of picks ',nt
write(*,*)  'number of P picks ',ntp
write(*,*)  'number of S picks ',nts
write(*,*) '======================================== '



do irec=1,nt
  read(10,rec=irec+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
  write(11,'(3i8,5f25.6,2i8)') id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
enddo

close(unit=10)
close(unit=11)

end program fsphe2a
