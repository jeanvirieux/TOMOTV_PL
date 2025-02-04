!###############################################################################|
! program a2fsrc    ASCII >>> fsrc                                              |
!
!   input (ascii file) 
!
!      fsrc.dat
!
!  output (binary file)                                                         |
!                                                                               |
!       fsrc        (sources)    8*4 = 32 bytes for each record                 |
! BINARY FORMAT     nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7                  |
!   record 1+irec   sr_x(irec),sr_y(irec),sr_z(irec),sr_to(irec),   &           |
!                   ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)             |
! 
!###############################################################################|

program a2fsrc
implicit none
integer(kind=4) :: nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7,irec
integer(kind=4) :: ki_pos,ki_to,ki_m,ki_id
real(kind=4) :: sr_x,sr_y,sr_z,sr_to
integer(kind=4) :: iopt

open(11,file='fsrc.asc',status='old')
open(10,file='fsrc',access='direct',recl=8*4)

write(*,*) ' file conversion from fsrc.asc into fsrc '
write(*,*) ' input iopt (0= keep event situation) (1= reset event situation)'
read(*,*) iopt

read(11,*) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7
if(iopt == 1) then
   neqks=neqks+neqks_out
   neqks_out=0
endif   
write(10,rec=1) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7

write(*,*) '========================================== '
write(*,*) 'number of sources ',nsrc
write(*,*) ' number of quakes inside the domain ', neqks
write(*,*) ' number of quakes outside the domain ', neqks_out
if(iopt == 1) write(*,*) ' if reset, should be zero '
write(*,*) ' number of quakes ',neqks+neqks_out
write(*,*) ' number of shots ',nshot
write(*,*) ' number of blasts ',nblast
write(*,*) '========================================== '

do irec=1,nsrc
  read(11,*) sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
  !  write(*,*) 'source ',irec,sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
  if(iopt==1) ki_m=1
  write(10,rec=irec+1) sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
enddo

close(unit=10)
close(unit=11)

end program a2fsrc
