! ###############################################################################     
!    reading of the source file fsrc
!
!    input   : usrc,nsrc
!   usrc     : unite logique du fichier fsrc
!   nsrc     : nombre de sources
!
!    output  : following values 
!
!   sr_x(:)  : coordonnees des sources + temps origine     
!   sr_y(:)  :
!   sr_z(:)  :
!   sr_to(:) :
!
!   ki_pos(:): flag for inverting positions   1 = yes
!   ki_to(:) : flag for inverting origin time 1 = yes
!   ki_m(:)  : flag for the source inside the box 1=in
!   ki_id(:) : ID of the source (should be set once for an experiment)
!
!------------------------------------------------------------------------------- 
MODULE s_read_src
contains

subroutine read_fsrc(usrc,sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id,nsrc)
implicit none

integer(kind=4) :: nsrc,irec,usrc
integer(kind=4), dimension(:) :: ki_pos(nsrc),ki_to(nsrc),ki_m(nsrc),ki_id(nsrc)
real(kind=4), dimension(:) :: sr_x(nsrc),sr_y(nsrc),sr_z(nsrc),sr_to(nsrc)

do irec=1,nsrc
  read(usrc,rec=irec+1) sr_x(irec),sr_y(irec),sr_z(irec),sr_to(irec), &
                        ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)
enddo

return
end subroutine read_fsrc

END MODULE s_read_src




