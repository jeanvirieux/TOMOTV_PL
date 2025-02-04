!###############################################################################|
! program fsrc2in    fsrc >>> fsrc.ins    put back events into the list         |
!
!   input (binary file)                                                         |
!                                                                               |
!       fsrc        (sources)    8*4 = 32 bytes for each record                 |
! BINARY FORMAT     nsrc,neqks,nshot,nblast,i4,i5,i6,i7                         |
!   record 1+irec   sr_x(irec),sr_y(irec),sr_z(irec),sr_to(irec),   &           |
!                   ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)             |
! 
!   output (ascii file)
!
!       fsrc.asc
!
!
!
! RAPPEL voir raytrace.tomo.f90
!
!---fsrc (binaire rec=32) fichier decrivant les sources-------------------------|  
!                                                                               |
!             +---------------------------------------------------------------+ |
! rec=1     : |nsr   |neqks |nshot |nblast |neqks_out|   0   |   0    |  0    | |
! ...         +------+------+------+-------+---------+-------+--------+-------| |
! rec=i+1   : | x(i) | y(i) | z(i) | to(i) | kp(i)   | kt(i) | in(i)  |id(i)  | |
! ...         +------+------+------+-------+---------+-------+--------+-------| |
! rec=nsr+1 : |x(nsr)|y(nsr)|z(nsr)|to(nsr)| kp(nsr) |kt(nsr)| in(nsr)|id(nsr)| |
!             +------+------+------+-------+---------+-------+--------+-------| |
!                                                                               |
! nsr       integer(kind=4): nombre total de sources                            |
! neqks     integer(kind=4): nbre de seismes(position et temps origine inconnus)|
! nshot     integer(kind=4): nbre de tir(position et temps origine connus)      |
! nblast    integer(kind=4): nbre d'explosions(temps origine inconnu)           |
! neqks_out integer(kind=4): nbre de seisme hors du domaine de calcul           |
!                                                                               |
!    x(j),y(j),z(j) real(kind=4): coordonnee de la source j (j<src+1)           |
!    to(j)          real(kind=4): temps origine de la source j                  | 
!    kp(j)       integer(kind=4): index sur la position des sources             |
!                                                                               | 
!           kp est une fonction de l'ensemble des sources a valeurs entieres    |
!           qui verifie :                                                       |
!                                                                               |
!                      (i)  kp(j)=0 sauf pour les seismes.                      |  
!                      (ii) la restriction de kp a l'ensemble des seismes       |
!                           est une fonction strictement croissante avec des    |
!                           valeurs dans {1,2,..,eqks}, ce qui permet de donner |
!                           la position dans la matrice des derivees            |
!                                                                               |
!    kt(j)     integer(kind=4)  : index sur le temps origine des sources        |
!                                                                               | 
!           kt est une fonction de l'ensemble des sources a valeurs entieres    |
!           qui verifie :                                                       |
!                                                                               |
!                      (i)  kt(j)=0 sauf pour les seismes et explosions.        |  
!                      (ii) la restriction de kp a l'ensemble des seismes       |
!                           et explosions est une fonction strictement          |
!                            croissante, a valeurs dans {1,2,..,eqks+shot}      |
!                           kt(j) est une fonction croissante qui               |
!                           numerote les seismes et les explosions              |
!                                                                               |
!    in(j)     integer(kind=4)  : (in(j)=1=>source j dans le domaine)           |
!                           (in(j)=0=>source j hors du domaine)                 |
!                                                                               |
!    id(j=)    integer(kind=4)  : un ID unique par source                       |
!
!  kp == ki_pos;  kt == ki_to;    in == ki_m;    id == ki_id
!###############################################################################|
program fsrc2in

implicit none
integer(kind=4) :: nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7
integer(kind=4) :: ki_pos,ki_to,ki_m,ki_id,irec,icount
real(kind=4) :: sr_x,sr_y,sr_z,sr_to
real(kind=4) :: xmin,xmax,ymin,ymax,zmin,zmax

xmin=1.e29
xmax=-1.e29
ymin=1.e29
ymax=-1.e29
zmin=1.e29
zmax=-1.e29

write(*,*) ' file conversion fsrc into fsrc.asc '

open(10,file='fsrc',access='direct',recl=8*4)
open(11,file='fsrc.ins',status='unknown')

read(10,rec=1) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7

write(*,*) '======================================== '
write(*,*)  'number of sources ',nsrc
write(*,*) ' number of quakes inside the domain ', neqks
write(*,*) ' number of quakes outside the domain ', neqks_out
write(*,*) ' number of quakes ',neqks+neqks_out
write(*,*) ' number of shots ',nshot
write(*,*) ' number of blasts ',nblast
write(*,*) '======================================== '

neqks=neqks+neqks_out
neqks_out=0
write(11,*) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7

icount=0
do irec=1,nsrc
   read(10,rec=irec+1) sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
   if(ki_m == 0) icount=icount+1
!!!!!!!!!!!!!!!! put the event into the active list    ki_m=1
   write(11,'(4f25.6,4i8)') sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,1,ki_id
  if(xmin >= sr_x) xmin=sr_x
  if(ymin >= sr_y) ymin=sr_y
  if(zmin >= sr_z) zmin=sr_z
  if(xmax <= sr_x) xmax=sr_x
  if(ymax <= sr_y) ymax=sr_y
  if(zmax <= sr_z) zmax=sr_z
enddo

write(*,*) ' extension of the box for sources '
write(*,*) ' xmin,xmax ',xmin,xmax
write(*,*) ' ymin,ymax ',ymin,ymax
write(*,*) ' zmin,zmax ',zmin,zmax
write(*,*) ' ================================ '
write(*,*) ' number of frozen events put back into the active list',icount
write(*,*) ' ================================ '

close(unit=10)
close(unit=11)

end program fsrc2in
