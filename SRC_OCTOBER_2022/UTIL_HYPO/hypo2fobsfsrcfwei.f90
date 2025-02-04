!###############################################################################|
! PROGRAM hypo2fobsfsrcfwei     HYPO71 >>> fobs
!     input  
!       file.prt where are locations AND stations data                          !
!     output (binary file)
!       fsrc        (events)        8*4 = 32 bytes for each record
!       fobs        (observation)   6*4 = 24 bytes for each record              |
!       fwei       (weight)        1*4
!
!   7 December 2010 Ortensia: ATTENTION fwei is WRONG!!  
!                                          it contains one more value    
!
! BINARY FORMAT     nt,ntp,nts,i1,12,13                                         |                 
!   record 1+irec   id_obs,t_lu,dt_lu,id_src,id_sta,id_ray                      |
!
!
!   MODIFICATION BY Ortensia AMOROSO    December 2010
!
!###############################################################################|
program hypo2fobsfsrcfwei
implicit none

integer(kind=4) i1,i2,i3,iempty,ierr

character(len=132) :: name_hypo
character(len=132) :: string
character(len=4) :: name
character(len=1) :: pickp

character(len=4), allocatable, dimension(:), target :: sta_name
integer(kind=4), allocatable, dimension(:), target :: sta_id

integer(kind=4) id_obs,id_src,id_sta,id_ray


integer(kind=4) :: nsrc,nshots,neqks,neqks_out,nblasts
integer(kind=4) :: ista,nsta,ista_sel,iobs,isrc

integer(kind=4) :: iquality,iweightp,weightp,iweights,weights

integer(kind=4) :: kp_src,kt_src,in_src

integer(kind=4) :: nt,ntp,nts
real(kind=4) :: t_lu,dt_lu,dt,dtp_sta,dts_sta
real(kind=4) :: x_src,y_src,z_src,t0_src
real(kind=4) :: x_sta,y_sta,z_sta,t_sta
real(kind=4) :: weip,weis
integer(kind=4) :: yr,mt,dy,hr,mn
real(kind=4) :: sec
real(kind=8) :: t0_src_d

integer(kind=4) :: RefEllipsoid
character(len=1) ::  UTMZone
integer(kind=4) :: ZoneNumber

open(55,file='check_utm')
call constants    ! set constants for UTM projection
RefEllipsoid = 24 ! 24 is WGS-84. 

!*********************************************
! we need fsta file and we need to create a lut for id_sta and name
!*********************************************
open(10,file='fsta',access='direct',recl=7*4) ! we open the stream 10 on file fsta

read(10,rec=1) nsta
write(*,*) ' ==========================='
write(*,*) ' number of stations ', nsta
write(*,*) ' ==========================='
! allocate memory space for the LUT station
allocate(sta_id(nsta))
allocate(sta_name(nsta))

do ista=1,nsta
read(10,rec=1+ista) id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
sta_name(ista)=name
sta_id(ista)=id_sta
enddo
close(10)

!*********************************************
! read the output file of HYPO : ****.prt
!*********************************************
write(*,*) ' enter the output HYPO71 file *.prt'
read(*,'(a)') name_hypo
open(9,file=name_hypo,status='old')

write(*,*) ' enter the quality of location '
write(*,*) ' 1 means only A '
write(*,*) ' 2 means only A+B '
write(*,*) ' 3 means only A+B+C '
write(*,*) ' 4 means only A+B+C+D '
read(*,*) iquality
! ORTENSIA peso pick P
write(*,*) ' enter the quality of picking P '
write(*,*) ' 0 means take only high quality 0 '
write(*,*) ' 1 means take only 0+1 '
write(*,*) ' 2 means take only 0+1+2 '
write(*,*) ' 3 means take only 0+1+2+3 '
write(*,*) ' 4 means take every pick data 0+1+2+3+4 '
read(*,*) iweightp
! ORTENSIA peso pick S
write(*,*) ' enter the quality of picking S '
write(*,*) ' 0 means take only high quality 0 '
write(*,*) ' 1 means take only 0+1 '
write(*,*) ' 2 means take only 0+1+2 '
write(*,*) ' 3 means take only 0+1+2+3 '
write(*,*) ' 4 means take every pick data 0+1+2+3+4 '
read(*,*) iweights

open(11,file='fsrc',access='direct',recl=8*4) ! we open the stream 11 on file fsrc
open(12,file='fobs',access='direct',recl=6*4) ! we open the stream 12 on file fobs
open(13,file='fwei',access='direct',recl=4) ! we open the stream 13 on file fwei


nsrc=0;neqks=0;nshots=0;nblasts=0;neqks_out=0;i1=0;i2=0;i3=0
isrc=0
write(11,rec=1) nsrc,neqks,nshots,nblasts,neqks_out,i1,i2,i3  ! dummy line upto now

nt=0;ntp=0;nts=0;i1=0;i2=0;i3=0
iobs=0
write(12,rec=1) nt,ntp,nts,i1,12,13 ! dummy line upto now for fobs


id_ray=0   ! no information yet about the ray
t0_src=0.  ! always set to zero for each earthquake

!====================================
!   LOOP ON EVENTS
!====================================
2000  continue
read(9,'(a132)',end=5000) string

if(string(3:7) /= 'DATE ') goto 2000

!write(*,*) ' @@@@@@@@@@@@@@@@@@@ '
!write(*,'(a)') string
!
!  ONE EVENT !
!
isrc=isrc+1
read(9,'(a132)') string
!write(*,*) ' EVENT ', string(80:80)
!write(*,*) ' reading P ',iquality,string(80:80)
!================ if this quake has not the required quality, jump to another one
if(iquality == 1 .and. string(80:80) /= 'A') then
  goto 2000
elseif(iquality == 2 .and. (string(80:80) /= 'A' .and. string(80:80) /= 'B') ) then
  goto 2000
elseif(iquality == 3 .and. (string(80:80) /= 'A' .and. string(80:80) /= 'B' .and. string(80:80) /= 'C') ) then
  goto 2000
endif
!write(*,*) ' reading P ',iquality,string(80:80)

! extract positions lat lon ... do utm conversion
call read_string(string,x_src,y_src,z_src,t0_src_d,RefEllipsoid,UTMZone,Zonenumber)


kp_src=isrc    ! counter ... not the ID which is supposed to be defined by the data acquisition system
               ! because it is an earthqake (we have to find x_src,y_src,z_src
kt_src=isrc    ! for earthquakes        ( # between shots & blasts)
in_src=1       ! the quake is in the box !
id_src=isrc    ! id of the earthquake - should be provided by the database .. if not use the counter

!write(*,*) ' ************* '
!write(*,*) ' source ',id_src,kp_src,kt_src,in_src
!write(*,*)  x_src,y_src,z_src,t0_src_d 
!write(*,*) ' ************* '

write(11,rec=1+isrc) x_src,y_src,z_src,t0_src,kp_src,kt_src,in_src,id_src

!*****************************
!   LOOP ON RECEIVER FOR THE ISRC EVENT
!*****************************

3000 continue
read(9,'(a132)') string

if(string(3:6) /= 'STN ') goto 3000

!
! ONE RECEIVER !
!
4000 continue

read(9,'(a132)') string
! is it a blank line if yes jump
iempty=0
call test_empty(string,iempty)
if(iempty == 0) goto 2000  ! if blank jump to 2000 TO NEXT EVENT

! extract date or utime, weip

call read_fobs_P(string,t_sta,dt,t0_src_d,name,ierr,iweightp,weip)
! check the selection of the data
if(ierr /= 0) goto 4000 ! do not consider this picked time

!write(*,*) ' ############## '
!write(*,'(a)') string
!write(*,*) ' ############## '

iobs=iobs+1
t_lu=t_sta ! data traveltime as event starts at t=0.
dt_lu=dt

do ista=1,nsta
ista_sel=ista
if(name == sta_name(ista)) goto 3100
enddo
!
!   wrong fsta because the station name is missing
!
write(*,*) ' error in scanning the fsta file station name NOT FOUND'
stop

3100 continue
id_sta=sta_id(ista)   ! get the id of the station

id_obs=iobs           ! set the id of the time data

write(12,rec=1+iobs) id_obs,t_lu,dt_lu,id_src,id_sta,id_ray
write(13,rec=iobs) weip

!write(*,*) '-----------------------ORTENSIA------------------------'
!write(*,*) ' DATA P ', id_obs,t_lu,dt_lu,id_src,id_sta,id_ray
!write(*,*) ' weight ',weip 
!write(*,*) '-------------------------------------------------'
goto 4000

5000 continue

ntp=iobs   ! number of P picked times

!///////////////////////////////////////
!    END OF P TIME READING
!//////////////////////////////////////
nsrc=isrc

write(*,*) ' number of valid events ',nsrc
neqks=nsrc
write(11,rec=1) nsrc,neqks,nshots,nblasts,neqks_out,i1,i2,i3 
close(11)

!///////////////////////////////////////
!  DO NOW THE S READING 
!///////////////////////////////////////
rewind(9)   

isrc=0
!====================================
!   LOOP ON EVENTS
!====================================
2001  continue
read(9,'(a132)',end=5001) string

if(string(3:7) /= 'DATE ') goto 2001

!write(*,*) ' @@@@@@@@@@@@@@@@@@@ '
!write(*,'(a)') string
!
!  ONE EVENT !   KEEP THE COUNTER .... ALWAYS THE SAME WHATEVER IS YOUR SELECTION
!
isrc=isrc+1
read(9,'(a132)') string
!================ if this quake has not the required quality, jump to another one
if(iquality == 1 .and. string(80:80) /= 'A') then
  goto 2001
elseif(iquality == 2 .and. (string(80:80) /= 'A' .and. string(80:80) /= 'B') ) then
  goto 2001
elseif(iquality == 3 .and. (string(80:80) /= 'A' .and. string(80:80) /= 'B' .and. string(80:80) /= 'C') ) then
  goto 2001
endif

id_src=isrc  ! id of the earthquake - provided by the db.. here use the counter

!*****************************
!   LOOP ON RECEIVER FOR THE ISRC EVENT
!*****************************

3001 continue
read(9,'(a132)') string

if(string(3:6) /= 'STN ') goto 3001

!
! ONE RECEIVER !
!
4001 continue

read(9,'(a132)') string

! is it a blank line if yes jump
iempty=0
call test_empty(string,iempty)
if(iempty == 0) goto 2001  ! if blank jump to 2000 TO NEXT EVENT

! extract date or utime, iquality,
call read_fobs_S(string,t_sta,dt,t0_src_d,name,ierr,iweights,weis)
! check the selection of the data
if(ierr /= 0) goto 4001 ! do not consider this picked time

!write(*,*) ' ############## '
!write(*,'(a)') string

iobs=iobs+1
t_lu=t_sta ! data traveltime  as event start at t=0.
dt_lu=dt   ! mv estimation of dt into the dt_lu

do ista=1,nsta
ista_sel=ista
if(name == sta_name(ista)) goto 3101
enddo
!
!   wrong fsta because the station name is missing
!
write(*,*) ' as we have checked for P times, we should not arrive here ! '
stop

3101 continue
id_sta=sta_id(ista_sel)   ! get the id of the station

id_obs=iobs           ! set the id of the time data

write(12,rec=1+iobs) id_obs,t_lu,dt_lu,id_src,id_sta,id_ray
write(13,rec=iobs) weis

!write(*,*) '-----------------------ORTENSIA------------------------'
!write(*,*) ' DATA S ', id_obs,t_lu,dt_lu,id_src,id_sta,id_ray
!write(*,*) ' WEIGHT S ',weis 
!write(*,*) '-------------------------------------------------'

goto 4001

5001 continue

!////////////////////////////////////////////
!    END OF S TIME READING
!///////////////////////////////////////////

nt=iobs
nts=nt-ntp

write(12,rec=1) nt,ntp,nts,0.,0.,0.

close(9)
close(12)
close(13)

stop
end program hypo2fobsfsrcfwei

subroutine test_empty(string,iempty)
!
!  use of the function len_trim   iempty=0 means the line is blank
!
character*(*) :: string
integer(kind=4) :: iempty

iempty=1   ! I assume that the line is not blank a priori
if(len_trim(string) == 0 ) then
  iempty=0
endif
return
end subroutine

subroutine read_string(string,x_src,y_src,z_src,t0_src_d,RefEllipsoid,UTMZone,Zonenumber)
character*(*) :: string
character(len=23) :: date
real(kind=4) :: x_src,y_src,z_src
real(kind=8) :: xlat_src_d,xlon_src_d,x_src_d,y_src_d

real(kind=8) :: t0_src_d

integer(kind=4) :: lat_src,lon_src,iz_src,yr,mt,dy,hr,mn
real(kind=4) :: sec

integer(kind=4) :: RefEllipsoid
character(len=1) ::  UTMZone
integer(kind=4) :: ZoneNumber

!@@@@@@@@@@@@@@@@@@@@@@@@@@
! LATLONG   TO UTM
!@@@@@@@@@@@@@@@@@@@@@@@@@@

read(string,'(19x,i2,1x,f5.2)') lat_src,xlat_src_d  
xlat_src_d=lat_src+xlat_src_d/60.   ! conversion in decimal degree 
read(string,'(29x,i2,1x,f5.2)') lon_src,xlon_src_d  
xlon_src_d=lon_src+xlon_src_d/60.   ! conversion in decimal degree


call lltoutm(RefEllipsoid,xlat_src_d, xlon_src_d, y_src_d, x_src_d,UTMZone,ZoneNumber)
!         X along Easting
!         Y along Northing

!write(*,*) ' #################### '
!write(*,*) 'xlat,xlon',xlat_src_d,xlon_src_d
!write(*,*) 'x_src,y_src',x_src_d,y_src_d
!write(*,*) ' #################### '

x_src=x_src_d
y_src=y_src_d

write(55,*) UTMZone,ZoneNumber

read(string,'(39x,f5.2)') z_src   ! depth in km
! convertion in meters
z_src=z_src*1000. 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! convert the string from HYPO71 into a string yyyy.mm.dd-hh:mn:ss.sss
!
! HYPO71  '  50817  840 57.97'
!          xxxxxxxxxxxxxxxxxx
!          123456789012345678
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@

date(1:2)='20'
read(string,'(1x,i2)') yr   ! turn around the problem of writing the year as a number that might be ' 5' instead of '05'
write(date(3:4),'(i2.2)') yr
date(5:5)='.'
date(6:7)=string(4:5)       ! month
date(8:8)='.'
date(9:10)=string(6:7)     ! day
date(11:11)='-'
read(string,'(8x,i2)') hr   ! turn around blank character problem instead of zero
write(date(12:13),'(i2.2)') hr
date(14:14)=':'
date(15:16)=string(11:12)    ! mn
date(17:17)=':'
date(18:22)=string(14:18)   ! sec in format xx.xx    ! format HYPO71
date(23:23)='0'             ! we put the 0 by hand  
call str2utime1(date,t0_src_d)   ! deduce from the date of the event the universal time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! in fact not really needed from the output of HYPO71

!write(*,*) '==============='
!write(*,*) date
!write(*,*) t0_src_d
!write(*,*) '==============='

return
end subroutine

subroutine read_fobs_P(string,t_sta,dt,t0_src_d,name,ierr,iweightp,weip)

integer(kind=4) :: ierr,iweightp,weightp
character*(*) :: string
real(kind=4) :: t_sta,weip
character*(*) :: name
character(len=1) :: pickp
character(len=6) :: appo
character(len=23) :: date
real(kind=4) :: dt_tab(5),wei_tab0(5),wei_tab1(5),wei_tab2(5),wei_tab3(5),wei_tab4(5)

real(kind=8) :: t0_src_d,t_sta_d
!============ tabulation
!  of time uncertainties ... NOT USED
!  of weights between 1. down to 0.
!============
data dt_tab/0.05,0.075,0.15,0.35,0.75/

data wei_tab0/1.00,0.75,0.50,0.25,0.00/   ! should be red from a file
data wei_tab1/1.00,1.00,0.75,0.50,0.25/
data wei_tab2/1.00,1.00,1.00,0.75,0.50/
data wei_tab3/1.00,1.00,1.00,1.00,0.00/
data wei_tab4/1.00,1.00,1.00,1.00,1.00/

name(1:4)=string(2:5) ! name of the station

ierr=0 ! 0 means we select the data
t_sta=0.    ! default value

!
! check existence of the picking
! read the quality of the picking ... make a corresponding dt
read(string,'(21x,a1)') pickp
if(pickp /= 'P') then          ! this picking is not a P
  ierr=1
  goto 6000
endif
read(string,'(23x,i1)') weightp  ! read the quality of the picking 
!=============================== high quality
if(iweightp == 0 .and. weightp /= 0) then     ! not enough good for what it is asked
  ierr=1
  goto 6000
else
  dt=dt_tab(weightp+1)
  weip=wei_tab0(weightp+1)
  goto 6010
endif
!============================== going down
if(iweightp == 1 .and. (weightp /= 0 .and. weightp /= 1) ) then
  ierr=1
  goto 6000
else
  dt=dt_tab(weightp+1)
  weip=wei_tab1(weightp+1)
  goto 6010
endif
!============================== down
if(iweightp == 2 .and. (weightp /= 0 .and. weightp /= 1 .and. weightp /= 2) ) then
  ierr=1
  goto 6000
else
  dt=dt_tab(weightp+1)
  weip=wei_tab2(weightp+1)
  goto 6010
endif
!============================== down
if(iweightp == 3 .and. (weightp /= 0 .and. weightp /= 1 .and. weightp /= 2 .and. weightp /= 3) ) then
  ierr=1
  goto 6000
else
  dt=dt_tab(weightp+1)
  weip=wei_tab3(weightp+1)
  goto 6010
endif
!============================== down
if(iweightp == 4 .and. (weightp /= 0 .and. weightp /= 1 .and. weightp /= 2 .and. weightp /= 3 .and. weightp /=4) ) then
  ierr=1
  goto 6000
else
  dt=dt_tab(weightp+1)
  weip=wei_tab4(weightp+1)
  goto 6010
endif

if(iweightp > 4) then    ! remove other values in case of an error
  ierr=1
  goto 6000
endif
! ORTENSIA test per l'esistenza del pick P
6010 continue
read(string,'(35x,a6)') appo
if(appo /= '******') then
  read(string,'(36x,f5.2)') t_sta! time since the origin time of the earthquake ... IT IS PROPAGATING TIME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! we do not need to go back to universal time because we shall assume that
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! quakes start at zero time (at least before the inversion)
  t_sta_d=t0_src_d+t_sta         ! just for checking  should be commented when sure that it works
  call utime2str1(t_sta_d,date)
endif
6000 continue
return
end subroutine read_fobs_P

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    S reading
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_fobs_S(string,t_sta,dt,t0_src_d,name,ierr,iweights,weis)

integer(kind=4) :: ierr,weights,iweights
character*(*) :: string
real(kind=4) :: t_sta,dt,weis
character*(*) :: name
character(len=6) :: appo
character(len=1) :: picks
character(len=23) :: date
real(kind=4) :: dt_tab(5),wei_tab0(5),wei_tab1(5),wei_tab2(5),wei_tab3(5),wei_tab4(5)

real(kind=8) :: t0_src_d,t_sta_d

!============ tabulation
!  of time uncertainties ... NOT USED
!  of weights between 1. down to 0.
!============
data dt_tab/0.05,0.075,0.15,0.35,0.75/

data wei_tab0/1.00,0.75,0.50,0.25,0.00/   ! should be red from a file
data wei_tab1/1.00,1.00,0.75,0.50,0.25/
data wei_tab2/1.00,1.00,1.00,0.75,0.50/
data wei_tab3/1.00,1.00,1.00,1.00,0.00/
data wei_tab4/1.00,1.00,1.00,1.00,1.00/

name(1:4)=string(2:5) ! name of the station

ierr=0 ! 0 means we select the data
t_sta=0.    ! default value

! check existence of the picking
! read the type of the picking ... make a corresponding dt
read(string,'(100x,a1)') picks
if(picks /= 'S') then
  ierr=1
  goto 7000
endif
! read the quality of the picking
read(string,'(102x,i1)') weights
!=============================== high quality
if(iweights == 0 .and. weights /= 0) then     ! not enough good for what it is asked
  ierr=1
  goto 7000
else
  dt=dt_tab(weights+1)
  weis=wei_tab0(weights+1)
  goto 7010
endif
!============================== going down
if(iweights == 1 .and. (weights /= 0 .and. weights /= 1) ) then
  ierr=1
  goto 7000
else
  dt=dt_tab(weights+1)
  weis=wei_tab1(weights+1)
  goto 7010
endif
!============================== down
if(iweights == 2 .and. (weights /= 0 .and. weights /= 1 .and. weights /= 2) ) then
  ierr=1
  goto 7000
else
  dt=dt_tab(weights+1)
  weis=wei_tab2(weights+1)
  goto 7010
endif
!============================== down
if(iweights == 3 .and. (weights /= 0 .and. weights /= 1 .and. weights /= 2 .and. weights /= 3) ) then
  ierr=1
  goto 7000
else
  dt=dt_tab(weights+1)
  weis=wei_tab3(weights+1)
  goto 7010
endif
!============================== down
if(iweights == 4 .and. (weights /= 0 .and. weights /= 1 .and. weights /= 2 .and. weights /= 3 .and. weights /=4) ) then
  ierr=1
  goto 7000
else
  dt=dt_tab(weights+1)
  weis=wei_tab4(weights+1)
  goto 7010
endif

if(iweights > 4) then    ! remove other values in case of an error
  ierr=1
  goto 7000
endif

! ORTENSIA test per l'esistenza del pick S
7010 continue
read(string,'(109x,a6)') appo
if(appo /= '******') then
  read(string,'(110x,f5.2)') t_sta ! time since the origin time of the earthquake ... IT IS PROPAGATING TIME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! we do not need to go back to universal time because we shall assume that
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! quakes start at zero time (at least before the inversion)
  t_sta_d=t0_src_d+t_sta           ! just for checking  should be commented when sure that it works
  call utime2str1(t_sta_d,date)
endif

7000 continue
return
end subroutine read_fobs_S

