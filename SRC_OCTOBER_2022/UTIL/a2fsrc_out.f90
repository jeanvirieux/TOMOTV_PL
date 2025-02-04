!###############################################################################|
! program a2fsrc_out    ASCII >>> fsrc                                          |
!                                                                               |
!   input (ascii file)                                                          |
!                                                                               |
!      fsrc.asc                                                                 |
!                                                                               |
!   remove events outside the box given by model.head (only quakes)             |
!                                                                               |
!  output (binary file)                                                         |
!                                                                               |
!       fsrc        (sources)    8*4 = 32 bytes for each record                 |
! BINARY FORMAT     nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7                  |
!   record 1+irec   sr_x(irec),sr_y(irec),sr_z(irec),sr_to(irec),   &           |
!                   ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)             |
!                                                                               |
!###############################################################################|

program a2fsrc_out
implicit none
integer(kind=4) :: nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7,irec
integer(kind=4) :: ki_pos,ki_to,ki_m,ki_id
real(kind=4) :: sr_x,sr_y,sr_z,sr_to
integer(kind=4) :: irecord,iflag

character(len=50) :: carac

real(kind=4) :: xmin,xmax,ymin,ymax,zmin,zmax
real(kind=4) :: xorg,yorg,zorg,dx,dy,dz
integer(kind=4) :: nx,ny,nz

open(11,file='fsrc.asc',status='old')
open(10,file='fsrc',access='direct',recl=8*4)

open(7,file='model.head',status='old')
read(7,*) carac
read(7,*) carac
read(7,*) carac
read(7,*) xorg,yorg,zorg
read(7,*) carac
read(7,*) nx,ny,nz
read(7,*) carac
read(7,*) dx,dy,dz
close(7)

open(8,file='fsrc_id.out',status='unknown')

xmin=xorg+dx
xmax=xorg+float(nx-2)*dx
ymin=yorg+dy
ymax=yorg+float(ny-2)*dy
zmin=zorg+dz
zmax=zorg+float(nz-2)*dz

write(*,*) ' xmin,xmax ',xmin,xmax
write(*,*) ' ymin,ymax ',ymin,ymax
write(*,*) ' zmin,zmax ',zmin,zmax


read(11,*) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7
write(10,rec=1) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7

write(*,*) '========================================== '
write(*,*) 'number of sources ',nsrc
write(*,*) ' number of quakes inside the domain ', neqks
write(*,*) ' number of quakes outside the domain ', neqks_out
write(*,*) ' number of quakes ',neqks+neqks_out
write(*,*) ' number of shots ',nshot
write(*,*) ' number of blasts ',nblast
write(*,*) '========================================== '

irecord=0
do irec=1,nsrc
  read(11,*) sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
  !  write(*,*) 'source ',irec,sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
  iflag=1  ! selection
  if(ki_m == 0) then
     write(8,*) ki_id
     iflag=0
  endif   
  if(sr_x < xmin .or. sr_x > xmax .or. sr_y < ymin .or. sr_y > ymax .or. &
       sr_z < zmin .or. sr_z > zmax) then
     ki_m=0
     neqks=neqks-1
     neqks_out=neqks_out+1
     write(8,*) ki_id
     iflag=0
  endif
  if(iflag == 1) then
     irecord=irecord+1
     write(10,rec=irecord+1) sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
  endif   
enddo

nsrc=irecord
neqks_out=0
neqks=irecord    !only quakes
neqks_out=0

write(10,rec=1) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7

close(10)
close(11)
close(8)

write(*,*) '========================================== '
write(*,*) 'number of sources ',nsrc
write(*,*) ' number of quakes inside the domain ', neqks
write(*,*) ' number of quakes outside the domain ', neqks_out
write(*,*) ' number of quakes ',neqks+neqks_out
write(*,*) ' number of blasts ',nblast
write(*,*) '========================================== '
end program a2fsrc_out
