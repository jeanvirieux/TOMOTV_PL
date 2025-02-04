!###############################################################################|
! Program fobs2tt      fobs >>> ASCII
!
!          ===> correction of time if around (+-) 60 or (+-) 120 and so on
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
program fobs2tt
  implicit none

  integer(kind=4) :: irec,nt,ntp,nts,ndum1,ndum2,ndum3
  real(kind=4) :: t_lu,dt_lu
  integer(kind=4) :: id_dat,lut_src,lut_sta,lut_ray

  integer(kind=4) :: id_cor1,icor,ncor
  real(kind=4) :: tobs1,tsyn1,tdif1

  integer(kind=4),allocatable,dimension(:) :: id_cor
  real(kind=4),allocatable,dimension(:)  :: tobs,tsyn,tdif

  real(kind=4) :: seuil
  integer(kind=4) :: nt_cor,ntp_cor,nts_cor

  write(*,*) ' file conversion fobs into fobs.ttt '
  write(*,*) ' an attempt for correcting wrong time conversion '
  write(*,*) ' around multiples of 60 sec'

  write(*,*) ' enter threshold (sec) '
  read(*,*) seuil

!!!====================================== reading the file with odd times
  open(9,file='fobs.cor',status='old')
  icor=0
100 continue
  read(9,*,end=200) id_cor1,tobs1,tsyn1,tdif1
  icor=icor+1
  goto 100
200 continue
  close(9)
  ncor=icor
  write(*,*) ' number of wrong pickings ',ncor

  allocate(id_cor(ncor))
  allocate(tobs(ncor))
  allocate(tsyn(ncor))
  allocate(tdif(ncor))
  open(9,file='fobs.cor',status='old')
  do icor=1,ncor
     read(9,*) id_cor(icor),tobs(icor),tsyn(icor),tdif(icor)
  enddo
  close(9)

  open(10,file='fobs',access='direct',recl=24)
  open(11,file='junk.ttt',status='unknown')

  ntp_cor=0 
  nts_cor=0

  read(10,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3
  write(11,'(6I8)') nt,ntp,nts,ndum1,ndum2,ndum3
  do irec=1,ntp ! lecture P
     read(10,rec=irec+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
     ntp_cor=ntp_cor+1
     do icor=1,ncor
        if(id_dat == id_cor(icor)) then
           if(abs(tdif(icor)-60.) < seuil) then
              t_lu=t_lu-60.;goto 1
           endif
           if(abs(tdif(icor)+60.) < seuil) then
              t_lu=t_lu+60.;goto 1
           endif
           if(abs(tdif(icor)-120.) < seuil) then
              t_lu=t_lu-120.;goto 1
           endif
           if(abs(tdif(icor)+120.) < seuil) then
              t_lu=t_lu+120.;goto 1
           endif
           if(abs(tdif(icor)-180.) < seuil) then
              t_lu=t_lu-180.;goto 1
           endif
           if(abs(tdif(icor)+180.) < seuil) then
              t_lu=t_lu+180.;goto 1
           endif
           if(abs(tdif(icor)-240.) < seuil) then
              t_lu=t_lu-240.;goto 1
           endif
           if(abs(tdif(icor)+240.) < seuil) then
              t_lu=t_lu+240.;goto 1
           endif
           if(abs(tdif(icor)-300.) < seuil) then
              t_lu=t_lu-300.;goto 1
           endif
           if(abs(tdif(icor)+300.) < seuil) then
              t_lu=t_lu+300.;goto 1
           endif
           if(abs(tdif(icor)-360.) < seuil) then
              t_lu=t_lu-360.;goto 1
           endif
           if(abs(tdif(icor)+360.) < seuil) then
              t_lu=t_lu+360.;goto 1
           endif
           if(abs(tdif(icor)-420.) < seuil) then
              t_lu=t_lu-420.;goto 1
           endif
           if(abs(tdif(icor)+420.) < seuil) then
              t_lu=t_lu+420.;goto 1
           endif
           id_dat=0   ! really wrong remove
           ntp_cor=ntp_cor-1
           goto 1
        endif
     enddo
1    continue  
     if(id_dat /= 0) write(11,'(I8,2F25.15,3I8)') id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
  enddo
  do irec=ntp+1,nt ! lecture S
     read(10,rec=irec+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
     nts_cor=nts_cor+1
     do icor=1,ncor
        if(id_dat == id_cor(icor)) then
           if(abs(tdif(icor)-60.) < seuil) then
              t_lu=t_lu-60.;goto 2
           endif
           if(abs(tdif(icor)+60.) < seuil) then
              t_lu=t_lu+60.;goto 2
           endif
           if(abs(tdif(icor)-120.) < seuil) then
              t_lu=t_lu-120.;goto 2
           endif
           if(abs(tdif(icor)+120.) < seuil) then
              t_lu=t_lu+120.;goto 2
           endif
           if(abs(tdif(icor)-180.) < seuil) then
              t_lu=t_lu-180.;goto 2
           endif
           if(abs(tdif(icor)+180.) < seuil) then
              t_lu=t_lu+180.;goto 2
           endif
           if(abs(tdif(icor)-240.) < seuil) then
              t_lu=t_lu-240.;goto 2
           endif
           if(abs(tdif(icor)+240.) < seuil) then
              t_lu=t_lu+240.;goto 2
           endif
           if(abs(tdif(icor)-300.) < seuil) then
              t_lu=t_lu-300.;goto 2
           endif
           if(abs(tdif(icor)+300.) < seuil) then
              t_lu=t_lu+300.;goto 2
           endif
           if(abs(tdif(icor)-360.) < seuil) then
              t_lu=t_lu-360.;goto 2
           endif
           if(abs(tdif(icor)+360.) < seuil) then
              t_lu=t_lu+360.;goto 2
           endif
           if(abs(tdif(icor)-420.) < seuil) then
              t_lu=t_lu-420.;goto 2
           endif
           if(abs(tdif(icor)+420.) < seuil) then
              t_lu=t_lu+420.;goto 2
           endif
           id_dat=0   ! really wrong remove it
           nts_cor=nts_cor-1
        endif
     enddo
2    continue
     if(id_dat /= 0) write(11,'(I8,2F25.15,3I8)') id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
  enddo

  nt_cor=ntp_cor+nts_cor

  write(*,*) '============================'
  write(*,*) ' previous fobs nt,ntp,nts ',nt,ntp,nts
  write(*,*) ' acccepted fobs.ttt nt,ntp,nts',nt_cor,ntp_cor,nts_cor
  write(*,*) '============================'
  write(*,*) ' final rejected picks ',nt-nt_cor
  write(*,*) ' re-integrated picks ',ncor-(nt-nt_cor)
  write(*,*) '============================'

  close(10)
  close(11)

  open(10,file='junk.ttt',status='old')
  open(11,file='fobs.ttt',status='unknown')

  read(10,*) nt,ntp,nts,ndum1,ndum2,ndum3
  write(11,'(6I8)') nt_cor,ntp_cor,nts_cor,ndum1,ndum2,ndum3

3000 continue
  read(10,*,end=4000) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
  write(11,'(I8,2F25.15,3I8)') id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
  goto 3000
4000 continue
  close(10)
  close(11)

  !call system('/bin/rm junk.ttt')
  deallocate(id_cor)
  deallocate(tobs)
  deallocate(tsyn)
  deallocate(tdif)

  stop
end program fobs2tt
