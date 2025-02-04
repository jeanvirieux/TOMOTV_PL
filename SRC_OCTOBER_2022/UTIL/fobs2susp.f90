!###############################################################################|
! Program fobs2a      fobs >>> ASCII
!
!       input file (binary)
!       fobs        (data arrival times) 6*4 = 24 bytes for each record         |
! BINARY FORMAT     nt,ntp,nts,i1,i2,i3                                         |
!   record 1+irec   id_dat(irec),temps(irec),dtemps(irec), &                    |
!                   lut_src(irec),lut_sta(irec),lut_ray(irec)                   |
!                                                                               |
!---fobs (binaire rec=24) fichier contenant les temps d'arrivee-----------------| 
!             +--------+---------+---------+-----------+-----------+-----------+|
! rec=1    :  | nt     | ntp     | nts     |   0       |     0     |     0     ||
!   ...       +--------+---------+---------+-----------+-----------+-----------+| 
! rec=i+1  :  | id(i)  |  t(i)   | dt(i)   |lut_src(i) |lut_sta(i) |lut_ray(i) || 
!   ...       +--------+---------+---------+-----------+-----------+-----------+|
! rec=nt+1 :  | id(nt) |  t(nt)  | dt(nt)  |lut_src(nt)|lut_sta(nt)|lut_ray(nt)||
!             +--------+---------+---------+-----------+-----------+-----------+|
!       output file (ASCII)
!       fobs.asc
!
!###############################################################################|
program fobs2susp
implicit none

integer(kind=4) :: irec,nt,ntp,nts,ndum1,ndum2,ndum3
real(kind=4) :: t_lu,dt_lu
integer(kind=4) :: id_dat,lut_src,lut_sta,lut_ray

integer(kind=4) :: nwrong,iwrong
integer(kind=4),dimension(10000) :: id_susp
real(kind=4), dimension(10000) :: t_obs,t_syn,res

open(8,file='suspicious.dat',status='old')

iwrong=0
1000 continue
iwrong=iwrong+1
read(8,*,end=2000) id_susp(iwrong),t_obs(iwrong),t_syn(iwrong),res(iwrong)
goto 1000
2000 continue
nwrong=iwrong-1

write(*,*) ' entree suspicieuse ',nwrong

write(*,*) ' file conversion fobs into fobs.asc with correction 60 or 120 s'
write(*,*) ' adding a fwei equal to 1 '

open(10,file='fobs',access='direct',recl=24)
open(11,file='fobs.asc',status='unknown')
open(12,file='fwei',access='direct',recl=4)

read(10,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3
write(11,'(6I8)') nt,ntp,nts,ndum1,ndum2,ndum3

do irec=1,nt ! lecture 
   read(10,rec=irec+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray

   do iwrong=1,nwrong
      if(id_dat == id_susp(iwrong)) then
         if(res(iwrong) < 75.) t_lu=t_lu-60.
         if(res(iwrong) > 75.) t_lu=t_lu-120.
         goto 4000
      endif
   enddo
   
4000 continue
  write(11,'(I8,2F25.6,3I8)') id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
  write(12,rec=irec) 1.
enddo 

close(10)
close(11)
close(12)

stop
end program fobs2susp
