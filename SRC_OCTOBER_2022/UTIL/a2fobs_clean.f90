!###############################################################################|
! Program a2fobs_clean      ASCII >>> fobs 
!
!       input file (ASCII)
!       fobs.asc
!
!       output file (binary)
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
!###############################################################################|
program a2fobs_clean
  implicit none

  integer(kind=4) :: irec,nt,ntp,nts,ndum1,ndum2,ndum3
  real(kind=4) :: t_lu,dt_lu
  integer(kind=4) :: id_dat,lut_src,lut_sta,lut_ray

  integer(kind=4),parameter :: ibadp=5000
  integer(kind=4), dimension(ibadp) :: id_out
  integer(kind=4) :: ibad,isrc,iflag,irecord

  open(11,file='fobs.asc',status='old')
  open(10,file='fobs',access='direct',recl=24)
  open(12,file='fwei',access='direct',recl=4)

  open(8,file='fsrc_id.out',status='old')
  ibad=1
1000 continue
  read(8,*,end=2000) id_out(ibad)
  ibad=ibad+1
  if(ibad > 5000) then
     write(*,*) ' increase the expected number of rejected events in the code', ibad,ibadp
     stop
  endif   
  goto 1000
2000 continue
  close(8)

  write(*,*) ' number of events to be removed ',ibad

  write(*,*) ' file conversion from fobs.asc into fobs '
  write(*,*) ' adding a weight equal to one : should use dt ...'

  read(11,*) nt,ntp,nts,ndum1,ndum2,ndum3
  write(*,*) ' nt,ntp,nts ',nt,ntp,nts
  write(10,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3

  irecord=0

  do irec=1,nt 
     read(11,*) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
     iflag=1
     do isrc=1,ibad
        if(id_out(isrc) == lut_src) then
           if(irec <= ntp) then
              iflag=0
              ntp=ntp-1
           else
              iflag=0
              nts=nts-1
           endif
           goto 500
        endif
     enddo
500  continue
     
     if(iflag == 1) then
        irecord=irecord+1
        write(10,rec=irecord+1) irecord,t_lu,dt_lu,lut_src,lut_sta,lut_ray
        write(12,rec=irecord) 1.
     endif
  enddo

  nt=ntp+nts
  write(*,*) ' number of recorded data ',irecord
  write(*,*) ' number of events after substration ',nt
  write(*,*) ' summary '
  write(*,*) ' nt,ntp,nts ',nt,ntp,nts,'ntp+nts=',ntp+nts, ' = nt?'
  
  write(10,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3

  close(10)
  close(11)
  close(12)

  stop
end program a2fobs_clean
