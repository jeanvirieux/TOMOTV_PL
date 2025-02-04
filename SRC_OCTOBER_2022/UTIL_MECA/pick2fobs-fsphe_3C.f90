!###############################################################################|
! Program obs2fobs-fsph      
!        
!
!       input file (ASCII)
!       picking files from GEOPSY
!       network.txt
!       fdate.txt
!       pick_event.txt
!
!       output file (binary)
!       fobs        (data arrival times) 6*4 = 24 bytes for each record         |
! BINARY FORMAT     nt,ntp,nts,i1,i2,i3                                         |
!   record 1+irec   id_dat(irec),temps(irec),dtemps(irec), &                    |
!                   lut_src(irec),lut_sta(irec),lut_ray(irec)                   |
!                                                                               |
!---fobs (binaire rec=6*4) file with travel times and uncertainties ------------| 
!             +--------+---------+---------+-----------+-----------+-----------+|
! rec=1    :  | nt     | ntp     | nts     |   0       |     0     |     0     ||
!   ...       +--------+---------+---------+-----------+-----------+-----------+| 
! rec=i+1  :  | id(i)  |  t(i)   | dt(i)   |lut_src(i) |lut_sta(i) |lut_ray(i) || 
!   ...       +--------+---------+---------+-----------+-----------+-----------+|
! rec=nt+1 :  | id(nt) |  t(nt)  | dt(nt)  |lut_src(nt)|lut_sta(nt)|lut_ray(nt)||
!             +--------+---------+---------+-----------+-----------+-----------+|
!                                                                               |
!---fsphe (binary rec=10*4)                                                     |
!
!
!
!
!###############################################################################|
program obs2fobs
  implicit none

  integer, parameter :: neventp=10000,nstatp=1000

  integer(kind=4) :: irec,nt,ntp,nts,ndum1,ndum2,ndum3,icount
  real(kind=4) :: t_lu,dt_lu
  integer(kind=4) :: id_dat,lut_src,lut_sta,lut_ray

  character(len=256) :: string,substring(14)
  character(len=256) :: string_v,substring_v(14)
  character(len=256) :: string_e,substring_e(14)
  character(len=256) :: string_n,substring_n(14)


  integer(kind=4),dimension(neventp) :: id_src
  !  integer(kind=4),dimension(nstatp) :: id_sta                  ! for fsta not used replaced by network.txt (long name)
  integer(kind=4),dimension(nstatp) :: id_sta_ref

  character(len=8),dimension(nstatp) :: name_sta_long
  !  character(len=4),dimension(nstatp) :: name_sta               ! for fsta not used replaced by network.txt (long name)
  character(len=8) :: name_sta_cur

  character(len=256),dimension(neventp) :: event_lock_key
  integer(kind=4) :: nsrc_located

  character(len=256),dimension(neventp) :: event_pick_loc,event_pick_pick
  character(len=256) :: event_loc,event_pick
  integer(kind=4) :: nsrc_event_pick

  character(len=23) :: date_string
  integer(kind=4) :: ndata,idata,n_pick,n_nopick

  integer(kind=4) :: i1,i2,i3,i4,i5,i6,i7
  integer(kind=4) :: irecp,irecs,isrc,ista
  integer(kind=4) :: ilut,lut_found

  real(kind=8) :: utime_event,utime_pick
  real(kind=8),dimension(neventp) :: utime_src_ref,utime_src
  integer(kind=4),dimension(neventp) :: id_src_ref
  character(len=256),dimension(neventp) :: event_loc_key

  real(kind=4) :: azimuth,dip,px,py,pz
  !  real(kind=4) :: dtp_sta,dts_sta                     !  for fsta not used replaced by network.txt (long name)
  real(kind=4) :: sr_x,sr_y,sr_z,sr_to
  real(kind=4) :: x_sta,y_sta,z_sta
  real(kind=4) :: time_pick

  integer(kind=4) :: ki_m,ki_pos,ki_to
  integer(kind=4) :: impulse,iphase
  integer(kind=4) :: nblast,neqks,neqks_out,nshot,nsrc,nsta,nsta_network

  !============================== set the format for date string
  date_string(5:5)='.'
  date_string(8:8)='.'
  date_string(11:11)='-'
  date_string(14:14)='.'
  date_string(17:17)='.'
  date_string(20:23)='.000'

  !=========================== consider two temporary files

  open(11,file='ftem',access='direct',recl=6*4)
  open(12,file='fpro',access='direct',recl=10*4)

  write(*,*) ' file conversion from picking files  into fobs '
  write(*,*) ' adding a weight equal to one : should use dt when available ...'

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! check the different links   network.txt for stations
  !                             fdate.txt for event location ... find the utime for fsrc
  !                             pic_event.txt for link between a located event and the pick file
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !=========================================== build ids of sources and stations
  !====================== check stations: get IDs from network.txt file where long names are
  open(10,file='network.txt',status='old')
  nsta_network=0
100 continue
  read(10,'(a)',end=101) string
  call pick_line_space(string,substring)
  nsta_network=nsta_network+1
  name_sta_long(nsta_network)=trim(substring(1))
  read(substring(2),*) id_sta_ref(nsta_network)            ! get the series of IDs of stations
  goto 100
101 continue
  close(10)
  write(*,*) ' end of reading network.txt '
  write(*,*) ' number of referenced stations ', nsta_network
  !=========================================================== pas forcement necessaire sauf coherence
  !  open(10,file='fsta',access='direct',recl=7*4)
  ! read(10,rec=1) nsta,i1,i2,i3,i4,i5,i6
  ! do ista=1,nsta
  !   read(10,rec=1+ista) id_sta(ista),x_sta,y_sta,z_sta,dtp_sta,dts_sta,name_sta(ista)
  ! enddo
  ! close(10)
  !======================= done for stations
  !======================= check events
  open(10,file='fdate.txt',status='old')
  nsrc_located=0
102 continue
  read(10,'(a)',end=103) string
  call pick_line_space(string,substring)
  nsrc_located=nsrc_located+1
  read(substring(3),*) utime_src_ref(nsrc_located)
  event_loc_key(nsrc_located)=trim(substring(5))        ! store Location.20160822.114158 related to an id_src_ref
  read(substring(1),*) id_src_ref(nsrc_located)
  goto 102
103 continue
  close(10)
  write(*,*) ' end of reading fdate.txt '
  open(10,file='fsrc',access='direct',recl=8*4)
  read(10,rec=1) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7
  do isrc=1,nsrc
     read(10,rec=isrc+1) sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,id_src(isrc)
     do irec=1,nsrc_located
        if(id_src(isrc) ==  id_src_ref(irec) ) then
           utime_src(isrc)=utime_src_ref(irec)
           goto 104
        endif
     enddo ! irec
     write(*,*) ' error in detecting the ID of the event in the reference list of located events fdate.txt '
104  continue
  enddo ! isrc
  close(10)
  write(*,*) ' end of reading fsrc '
  !====================== done for events
  open(10,file='pick_event.txt',status='old')
  nsrc_event_pick=0
105 continue
  read(10,'(a)',end=106) string
  call pick_line_space(string,substring)
  nsrc_event_pick=nsrc_event_pick+1
  event_pick_loc(nsrc_event_pick)=trim(substring(1))               ! get Location.20160822.114158
  event_pick_pick(nsrc_event_pick)=trim(substring(2))              ! in relation with picking file
  goto 105
106 continue
  close(10)
  write(*,*) ' end of reading pick_event.txt '
  write(*,*) ' end of reading files for linkage '

  !=========================================== id_sta_ref and id_src and id_src_ref tables are defined
  !=========================================== from fsrc
  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
  !=========================================== fake values of number of picking
  nt=0
  ntp=0
  nts=0
  ndum1=0;ndum2=0;ndum3=0
  !=================== write the number of picking   into fobs and into fmeca
  write(11,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3
  write(12,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3,ndum1,ndum2,ndum3,ndum1
  id_dat=0   !====================== zero picking   start counting with new IDs ...
  irec=0      !====================== start writing at rec=2
  lut_ray=0     !====================== ray not yet specified
  px=-999.
  py=-999.
  pz=-999.
  azimuth=-999.
  dip=-999.
  impulse=0

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !=========================================== super loop over events
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  nsrc=1
  do isrc=1,nsrc

     write(*,*) ' ======================================== '
     write(*,*) '    Event: looking at its picking file for event ',isrc,' over total number ',nsrc
     write(*,*) ' ======================================== '
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !############################which event
     !  call find_id_event()   already defined by id_src
     !############################ we know which event and its origin time
     lut_src=id_src(isrc)
     utime_event=utime_src(isrc)
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     !##################################################################
     ! looking for a pick file
     !##################################################################

     do ilut=1,nsrc_located                         ! Is it a located event ?
        if(id_src(isrc) == id_src_ref(ilut)) then
           event_loc=event_loc_key(ilut)
           goto 500
        endif
     enddo
     write(*,*) ' unable to find the related event in fdate.txt '
     write(*,*) ' go to next event '
     goto 10000

500  continue

     do ilut=1,nsrc_event_pick                       ! OK then find the corresponding pick file
        if(event_loc == event_pick_loc(ilut)) then
           event_pick=event_pick_pick(ilut)
           goto 600
        endif
     enddo
     write(*,*) ' unable to find the related pick to this event in pick_event.txt '
     write(*,*) ' go to next event '
     goto 10000

600  continue

     write(*,*) ' scanning the pick file ',trim(event_pick)

     !############################scan the pick file for this event
     open(21,file=trim(event_pick),status='old')
     !=========================================== reading to detect the ID line
1000 continue
     read(21,'(a)',end=2000) string
     if(string(1:2) == 'ID') goto 2100
     goto 1000
2000 continue
     write(*,*) ' error in reading picking file for station: please check',isrc
     stop
     !========================================== ID has been detected
2100 continue
     ndata=0
     n_pick=0
     n_nopick=0
2110 continue
!=====================================================================================
!   be careful ... we assume three lines par picking
!=====================================================================================
     read(21,'(a)',end=2200) string_e               ! separation is TAB from GEOPSY
     read(21,'(a)',end=2200) string_n
     read(21,'(a)',end=2200) string_v
     call pick_line_tab(string_e,substring_e)
     call pick_line_tab(string_n,substring_n)
     call pick_line_tab(string_v,substring_v)
     ! 
     !  extension _e for EAST  _n for NORTH _v VERTICAL
     !substring(9) picked time for P up
     !substring(10) picked time for P down
     !substring(11) picked time for P w/o sign
     !substring(12) picked time for S up
     !substring(13) picked time for S down
     !substring(14) picked time for S w/o sign
     !
     !
     !========================================== have you got a picking (no 's' as criteria for picking)
     string=substring_e(9)                                              ! remove the TAB
     substring_e(9)=trim(string(1:len(trim(string))-1))
     string=substring_e(10)
     substring_e(10)=trim(string(1:len(trim(string))-1))
     string=substring_e(11)
     substring_e(11)=trim(string(1:len(trim(string))-1))
     string=substring_e(12)
     substring_e(12)=trim(string(1:len(trim(string))-1))
     string=substring_e(13)
     substring_e(13)=trim(string(1:len(trim(string))-1))
     string=substring_e(14)                                          
     substring_e(14)=trim(string(1:len(trim(string))-1))

     string=substring_n(9)                                              ! remove the TAB
     substring_n(9)=trim(string(1:len(trim(string))-1))
     string=substring_n(10)
     substring_n(10)=trim(string(1:len(trim(string))-1))
     string=substring_n(11)
     substring_n(11)=trim(string(1:len(trim(string))-1))
     string=substring_n(12)
     substring_n(12)=trim(string(1:len(trim(string))-1))
     string=substring_n(13)
     substring_n(13)=trim(string(1:len(trim(string))-1))
     string=substring_n(14)                                          
     substring_n(14)=trim(string(1:len(trim(string))-1))

     string=substring_v(9)                                              ! remove the TAB
     substring_v(9)=trim(string(1:len(trim(string))-1))
     string=substring_v(10)
     substring_v(10)=trim(string(1:len(trim(string))-1))
     string=substring_v(11)
     substring_v(11)=trim(string(1:len(trim(string))-1))
     string=substring_v(12)
     substring_v(12)=trim(string(1:len(trim(string))-1))
     string=substring_v(13)
     substring_v(13)=trim(string(1:len(trim(string))-1))
     string=substring_v(14)                                          
     substring_v(14)=trim(string(1:len(trim(string))-1))

     if( trim(substring_e(9)) /= '0s' .or. trim(substring_e(10)) /= '0s' .or. trim(substring_e(11)) /= '0s'  .or. &
          trim(substring_e(12)) /= '0s' .or. trim(substring_e(13)) /= '0s' .or. trim(substring_e(14)) /= '0s') then

        n_pick=n_pick+1
     write(76,*) '!',substring(9),'!'
     write(76,*) '!',substring(10),'!'
     write(76,*) '!',substring(11),'!'
     write(76,*) '!',substring(12),'!'
     write(76,*) '!',substring(13),'!'
     write(76,*) '!',substring(14),'!'
    else
     write(77,*) '!',substring(9),'!'
     write(77,*) '!',substring(10),'!'
     write(77,*) '!',substring(11),'!'
     write(77,*) '!',substring(12),'!'
     write(77,*) '!',substring(13),'!'
     write(77,*) '!',substring(14),'!'
        n_nopick=n_nopick+1
     endif
     ndata=ndata+1
     goto 2110
2200 continue
     write(*,*) ' for this event, the number of picking lines is ',ndata
     write(*,*) ' for this event, the number of picked  times is ',n_pick
     write(*,*) ' for this event, the number of unpicked times is ',n_nopick
     if(ndata == 0) goto 10000    ! move to next event as we have zero picking
     rewind(21)
     !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ RESTART
     !=========================================== reading to detect the ID line
1001 continue
     read(21,'(a)',end=2002) string
     if(string(1:2) == 'ID') goto 2001
     goto 1001
2002 write(*,*) ' error in reading this picking file ',trim(event_pick)
     goto 10000   ! next event
2001 continue
     !####################################################################
     !============================================== loop over data for this event
     !####################################################################
     do idata=1,ndata
        read(21,'(a)') string      ! should have a correct number

        write(*,*) idata
        write(*,*) string

        call pick_line_tab(string,substring)
        string=substring(9)                                            ! remove the TAB
        substring(9)=trim(string(1:len(trim(string))-1))
        string=substring(10)
        substring(10)=trim(string(1:len(trim(string))-1))
        string=substring(11)
        substring(11)=trim(string(1:len(trim(string))-1))
        string=substring(12)
        substring(12)=trim(string(1:len(trim(string))-1))
        string=substring(13)
        substring(13)=trim(string(1:len(trim(string))-1))
        !string=substring(14)                                          ! LAST without a TAB
        !substring(14)=trim(string(1:len(trim(string))-1))
        ! substring(9) picked time for P up  !substring(10) picked time for P down  !substring(11) picked time for P w/o sign
        ! substring(12) picked time for S up !substring(13) picked time for S down  !substring(14) picked time for S w/o sign
        !========================================== have you got a picking (no 's' as criteria for picking)






        if( trim(substring(9)) /= '0s' .or. trim(substring(10)) /= '0s' .or. trim(substring(11)) /= '0s'  .or. &
             trim(substring(12)) /= '0s' .or. trim(substring(13)) /= '0s' .or. trim(substring(14)) /= '0s') then
           !%%%%%%%%%%%%%%%%% which station?
           name_sta_cur=trim(substring(2))
           call find_id_station(name_sta_cur,name_sta_long,id_sta_ref,nsta_network,lut_found)
           if(lut_found < 0) then
              write(*,*) ' unable to find the station in the list '
              goto 9999    ! next data line
           endif
           lut_sta=lut_found
           !%%%%%%%%%%%%%%%%% which event?  already specified
           !################# which data?
           !     read(substring(1),*) id_data        ! local ID ... probably local to this picking ... not good for general ID of data
           string=substring(4)                 ! 23/08/2016   ! French style
           date_string(9:10)=trim(string(1:2)) !len(trim(string))-8) ! day
           date_string(6:7)=trim(string(4:5))  !len(trim(string))-5)  ! month
           date_string(1:4)=trim(string(7:10))  !len(trim(string)))    ! year
           string=substring(5)                 ! 21:04:08
           date_string(12:13)=trim(string(1:2)) !len(trim(string))-6)    ! hour
           date_string(15:16)=trim(string(4:5)) !len(trim(string))-4)    ! mn
           date_string(18:19)=trim(string(7:8)) !len(trim(string)))     ! up to the second
           call str2utime1(date_string,utime_pick)   ! conversion in utime upto the second of starting the recording
           !=========================================== which picked time?
           ! substring(9) picked time for P up  !substring(10) picked time for P down  !substring(11) picked time for P w/o sign
           ! substring(12) picked time for S up !substring(13) picked time for S down  !substring(14) picked time for S w/o sign
           icount=0
           string=substring(9)
           if( trim(string) /= '0s' ) then	
              substring(9)=trim(string(1:len(trim(string))-1))
              read(substring(9),*) time_pick
              impulse=1
              iphase=1
              ntp=ntp+1
              icount=1
           endif
           string=substring(10)
           if( trim(string) /= '0s' ) then
              substring(10)=trim(string(1:len(trim(string))-1))
              read(substring(10),*) time_pick
              impulse=-1
              iphase=1
              ntp=ntp+1
              icount=1
           endif
           string=substring(11)
           if( trim(string) /= '0s' ) then
              substring(11)=trim(string(1:len(trim(string))-1))
              read(substring(11),*) time_pick
              impulse=0
              iphase=1
              ntp=ntp+1
              icount=1
           endif
           !============================================================ S picking ... should be done in the second pass
           string=substring(12)
           if( trim(string) /= '0s' ) then
              substring(12)=trim(string(1:len(trim(string))-1))
              read(substring(12),*) time_pick
              impulse=1
              iphase=2
              nts=nts+1
              icount=1
           endif
           string=substring(13)
           if( trim(string) /= '0s' ) then
              substring(13)=trim(string(1:len(trim(string))-1))
              read(substring(13),*) time_pick
              impulse=-1
              iphase=2
              nts=nts+1
              icount=1
           endif
           string=substring(14)
           if( trim(string) /= '0s' ) then
              substring(14)=trim(string(1:len(trim(string))-1))
           write(*,*) 'substring 14 ',substring(14)
              read(substring(14),*) time_pick
              impulse=0
              iphase=2
              nts=nts+1
              icount=1
           endif
           if(icount == 0) goto 9999
           !============================================================
           utime_pick=utime_pick+time_pick         ! add these seconds to utime
           t_lu=utime_pick-utime_event             ! substract the origin time
           dt_lu=0.
!           write(*,*) nt,ntp,nts,t_lu,time_pick,utime_event,utime_pick
           !=================== write data in fobs      increment the id_dat
           id_dat=id_dat+1

           irec=irec+1         ! initial value =0, so we have 1 and write the record 2 initially
           write(11,rec=irec+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
           write(12,rec=irec+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase

        endif ! end of analysis of this pick time

9999    continue
        !####################################################################
        !============================================== loop over picking for this event
        !####################################################################
     enddo ! idata

10000 continue
     close(21)      ! close the pick file and move to the next event
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !=========================================== super loop over events
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  enddo    ! isrc

  nt=ntp+nts
  !=================== write the final number of pickings
  !      attention P ou S ...
  write(11,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3
  write(12,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3,ndum1,ndum2,ndum3,ndum1
  close(11)     ! ftem
  close(12)     ! fpro

  write(*,*) ' end of reading the pick file '
  write(*,*) ' reordering the pick file with phase '

  open(11,file='ftem',access='direct',recl=6*4)
  open(12,file='fpro',access='direct',recl=10*4)
  open(13,file='fwei',access='direct',recl=4)

  open(14,file='fobs',access='direct',recl=6*4)
  open(15,file='fsphe',access='direct',recl=10*4)

!================================================== first line

  write(14,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3
  write(15,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3,ndum1,ndum2,ndum3,ndum1

  irec=0
  irecp=0
  irecs=0
  do idata=1,ntp
     irec=irec+1
     read(11,rec=irec+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
     !     write(*,*) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
     read(12,rec=irec+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
     !     write(*,*) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
     if(iphase == 1) then                ! phase P
        irecp=irecp+1
        write(14,rec=irecp+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
        write(15,rec=irecp+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
        write(13,rec=irecp) 1.
     endif
  enddo
  !==================  write the following record for phase S
  do idata=1,nts
     irec=irec+1
     read(11,rec=irec+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
     read(12,rec=irec+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
     if(iphase == 2) then                ! phase S
        irecs=irecs+1
        write(14,rec=irecp+irecs+1) id_dat,t_lu,dt_lu,lut_src,lut_sta,lut_ray
        write(15,rec=irecp+irecs+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase
        write(13,rec=irecp+irecs) 1.
     endif
  enddo
  write(*,*) ' total records ',irec,' P records ',irecp,' S records ',irecs

  close(11)     ! ftem
  close(12)     ! fpro
  close(13)     ! fwei
  close(14)     ! fobs
  close(15)     ! fsphe

  call system('rm ftem')
  call system('rm fpro')
  stop
end program obs2fobs

subroutine find_id_station(name_sta_cur,name_sta_long,id_sta_ref,nsta,lut_found)
  integer, parameter :: neventp=10000,nstatp=1000
  character(len=8),dimension(nstatp) :: name_sta_long
  character(len=4),dimension(nstatp) :: name_sta
  integer(kind=4),dimension(nstatp) :: id_sta_ref
  character(len=8) :: name_sta_cur
  integer(kind=4) :: nsta,lut_found

  integer(kind=4) :: ista

  do ista=1,nsta
    if(trim(name_sta_cur) == trim(name_sta_long(ista))) then
      lut_found=id_sta_ref(ista)    ! get the ID
      return
    endif
  enddo
  lut_found=-1
  return
end subroutine find_id_station
