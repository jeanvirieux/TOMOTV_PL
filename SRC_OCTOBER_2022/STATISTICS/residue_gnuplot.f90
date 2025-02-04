!###############################################################################|
! PROGRAM FOR STATISTICS ON RESIDUAL                                            | 
!                                                                               |
!  COMPUTATION OF RESIDUALS DISTRIBUTION FOR EACH STATION                       |
!     input                                                                     |
!                                                                               |
!       fobs        (data arrival times) 6*4 = 24 bytes for each record         |
! BINARY FORMAT     nt,ntp,nts,i1,i2,i3                                         |
!   record 1+ista   id_dat(ista),temps(ista),dtemps(ista), &                    |
!                   lut_src(ista),lut_sta(ista),lut_ray(ista)                   |
!                                                                               |
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+ista   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!                                                                               |
!       fresi       UNWEIGHTED residues for each observation                    |
!                                                                               |
!                                                                               |
!                  version aout 2011                                            |
!                  Jean Virieux                                                 |
!     split between P and S residues
!###############################################################################|
program residue_analysis
use s_read_obs


implicit none
character(len=4) :: name
character(len=132) :: name_res

character(len=4), allocatable, dimension(:), target :: sta_name
integer(kind=4), allocatable, dimension(:), target :: sta_id
integer(kind=4), allocatable, dimension(:), target :: sta_num
real(kind=4), allocatable, dimension(:), target :: sta_dtp,sta_dts

integer(kind=4) :: nsta,ista
integer(kind=4), allocatable, dimension(:), target :: id_dat
integer(kind=4) :: id_sta0
real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta

real(kind=4), allocatable, dimension(:), target :: temps,dtemps,res
integer(kind=4), allocatable, dimension(:), target :: lut_ray,lut_sta,lut_src
integer(kind=4) :: nt,ntp,nts

integer(kind=4) :: i1,i2,i3,iunit,idata

real(kind=4), allocatable, dimension(:), target :: dt_step_P, dt_step_S
integer(kind=4), allocatable, dimension(:), target :: nt_step_P, nt_step_S
integer(kind=4) :: it0,it,nt_hist_P,nt_hist_S,icountp,icounts,chx
real(kind=4) :: dtmin_P,dtmax_P,dt_hist_P
real(kind=4) :: dtmin_S,dtmax_S,dt_hist_S

real(kind=4), allocatable, dimension(:,:), target :: res_sta

write(*,*) ' P(0) or S (1) or P+S (2)  enter '
read(*,*) chx
write(*,*) ' enter the name of the residues file fresi, fresu, fresv, fresw '
read(*,'(a)') name_res
write(*,*) 'enter dt_hist_P and nt_hist_P (odd) '
read(*,*) dt_hist_P, nt_hist_P
write(*,*) 'enter dt_hist_S and nt_hist_S (odd) '
read(*,*) dt_hist_S, nt_hist_S

!
!  compute structure of the histogram both for P and S
!
!****************************************** P
if(chx == 0 .or. chx == 2) then
dtmin_P=-(nt_hist_P/2)*dt_hist_P
dtmax_P=(nt_hist_P/2)*dt_hist_P
write(*,*) ' dtmin_P, dtmax_P for P ',dtmin_P,dtmax_P
allocate(dt_step_P(nt_hist_P))
allocate(nt_step_P(nt_hist_P))

dt_step_P(1)=dtmin_P
do it=2,nt_hist_P
dt_step_P(it)=dt_step_P(it-1)+dt_hist_P
enddo 
write(*,*) ' dtmax_P, dt_step_P_max ',dtmax_P,dt_step_P(nt_hist_P)
endif

!****************************************** S
if(chx == 1 .or. chx == 2) then
dtmin_S=-(nt_hist_S/2)*dt_hist_S
dtmax_S=(nt_hist_S/2)*dt_hist_S
write(*,*) ' dtmin_S, dtmax_S for S ',dtmin_S,dtmax_S
allocate(dt_step_S(nt_hist_S))
allocate(nt_step_S(nt_hist_S))

dt_step_S(1)=dtmin_S
do it=2,nt_hist_S
dt_step_S(it)=dt_step_S(it-1)+dt_hist_S
enddo 
write(*,*) ' dtmax_S, dt_step_S_max ',dtmax_S,dt_step_S(nt_hist_S)
endif

!*********************************************
! we need fsta file and we need to save a lut for id_sta and name
!*********************************************
open(10,file='fsta',access='direct',recl=7*4) ! we open the stream 10 on file fsta
read(10,rec=1) nsta
write(*,*) ' ==========================='
write(*,*) ' number of stations ', nsta
write(*,*) ' ==========================='
! allocate memory space for the LUT station
allocate(sta_id(nsta))
allocate(sta_name(nsta))
allocate(sta_dtp(nsta))
allocate(sta_dts(nsta))
allocate(sta_num(nsta))
do ista=1,nsta
  read(10,rec=1+ista) id_sta0,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
  sta_name(ista)=name
!  write(*,'(a)') sta_name(ista)
  sta_id(ista)=id_sta0
enddo
close(10)

!*************************************************************************
! we read now the fobs file for picked time
! and         the fresi file for residues at all observations
!*************************************************************************

! fobs : id_dat,temps,dtemps,lut_src,lut_sta,lut_ray
iunit=10
open(iunit,file='fobs',access='direct',recl=6*4)
read(iunit,rec=1) nt,ntp,nts,i1,i2,i3   ! read the total number of data (split into P & S data)

write(*,*) ' time data nt, ntp, nts ',nt,ntp,nts

allocate(id_dat(nt))     ! id of the data
allocate(temps(nt))      ! both P and S waves   ntp first and then nts
allocate(dtemps(nt))
allocate(lut_src(nt))    ! id of the source
allocate(lut_sta(nt))    ! id of the station
allocate(lut_ray(nt))    ! lut for each data (=0 if no ray or iray if one ray)

call read_fobs(iunit,id_dat,temps,dtemps,lut_src,lut_sta,lut_ray,nt,ntp,nts)
close(iunit)

allocate(res(nt))        !  residues at each observation   -999. magic number for rejected data
open(iunit,file=name_res,access='direct',recl=4*nt)     ! write per float  residues
read(iunit,rec=1) res
close(iunit)

write(*,*) ' end of reading of fsta, fobs and fresi '

!####################################
! OPEN the shell script for gnuplot
!####################################
open(9,file='STAT_run_stat')
!,access='append')

!**************************************** P
if(chx == 0 .or. chx == 2) then
!
!  make the global histogramme in P
!
nt_step_P(:)=0
do idata=1,ntp
  if(res(idata) >= dtmin_P .and. res(idata) <= dtmax_P) then
    it0=1
    do it=2,nt_hist_P
      if(res(idata) > dt_step_P(it)) then
        it0=it
      endif
    enddo  
    nt_step_P(it0)=nt_step_P(it0)+1
  endif ! inside the histogram
enddo
open(iunit,file='STAT_histo_all_P')
do it=1,nt_hist_P
write(iunit,*) dt_step_P(it)+0.5*dt_hist_P,nt_step_P(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "STAT_histo_res_all_P.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"STAT_histo_all_P"'//' linecolor 1 with boxes'
write(9,'(a)') 'EOD'
endif

!********************************** S
if(chx == 1 .or. chx == 2) then
!
!  make the global histogram in S
!
nt_step_S(:)=0
do idata=ntp+1,ntp+nts
  if(res(idata) >= dtmin_S .and. res(idata) <= dtmax_S) then
    it0=1
    do it=2,nt_hist_S
      if(res(idata) > dt_step_S(it)) then
        it0=it
      endif
    enddo  
    nt_step_S(it0)=nt_step_S(it0)+1
  endif ! inside the histogram
enddo
open(iunit,file='STAT_histo_all_S')
do it=1,nt_hist_S
write(iunit,*) dt_step_S(it)+0.5*dt_hist_S,nt_step_S(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "STAT_histo_res_all_S.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"STAT_histo_all_S"'//' linecolor 3 with boxes'
write(9,'(a)') 'EOD'
endif

allocate(res_sta(nt,nsta))     ! save the entire set of residues over stations
res_sta(:,:)=-999.

!********************************* P
if(chx == 0 .or. chx == 2) then
!=============================================================
!   loop over data   P first
!=============================================================
sta_dtp(:)=0.    ! residues at each station is set to zero
sta_num(:)=0
do idata=1,ntp
if(res(idata) /= -999.) then
  do ista=1,nsta
  if(lut_sta(idata) == sta_id(ista)) then
    sta_num(ista)=sta_num(ista)+1
    sta_dtp(ista)=sta_dtp(ista)+res(idata)
    res_sta(idata,ista)=res(idata)
  endif
  enddo   ! end on station ista
endif
enddo ! end on data   in P
! make the average for each station and deduce the time static equal to minus ...
do ista=1,nsta
  if(sta_num(ista) /= 0) then
    sta_dtp(ista)=-sta_dtp(ista)/float(sta_num(ista))   ! static delay at the station
  else
    sta_dtp(ista)=-999.
  endif
enddo  ! end on ista
endif

!********************************* S
if(chx == 1 .or. chx == 2) then
!=======================================================
!   loop over S data
!=======================================================
sta_dts(:)=0.    ! residues at each station is set to zero
sta_num(:)=0

do idata=ntp+1,ntp+nts
if(res(idata) /= -999.) then   ! if we have estimated a residue for this data
  do ista=1,nsta
  if(lut_sta(idata) == sta_id(ista)) then
    sta_num(ista)=sta_num(ista)+1
    sta_dts(ista)=sta_dts(ista)+res(idata)
    res_sta(idata,ista)=res(idata)
  endif
  enddo   ! end on station ista
endif
enddo ! end on data   in S
! make the average for each station and deduce the time static equal to minus
do ista=1,nsta
  if(sta_num(ista) /= 0) then
    sta_dts(ista)=-sta_dts(ista)/float(sta_num(ista))   ! static delay at the station
  else
    sta_dts(ista)=-999.
  endif
enddo ! end on ista
endif

!
!  write now the new fsta with estimated statics
!
open(10,file='fsta',access='direct',recl=7*4) ! we open again the stream 10 on file fsta
open(11,file='fsta.dt_static',access='direct',recl=7*4) ! we open the stream 11 for writing
read(10,rec=1) nsta
write(11,rec=1) nsta
do ista=1,nsta
  read(10,rec=1+ista) id_sta0,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
  write(*,*) ' id station ',name,' static p',sta_dtp(ista),' static s',sta_dts(ista)
  write(11,rec=1+ista) id_sta0,x_sta,y_sta,z_sta,sta_dtp(ista),sta_dts(ista),name ! new statics
enddo
close(10)   ! end of the fsta
close(11)   ! end of the new fsta

!
!   LOOK OVER STATIONS   FOR INDIVUAL HISTOGRAMS
!
do ista=1,nsta

!************************************** P
if(chx == 0 .or. chx == 2) then
!
!  Make histograms for P at each station
!
nt_step_P(:)=0
icountp=0
do idata=1,ntp
  if(res_sta(idata,ista) >= dtmin_P .and. res_sta(idata,ista) <= dtmax_P) then
    icountp=icountp+1
    it0=1
    do it=2,nt_hist_P
      if(res_sta(idata,ista) > dt_step_P(it)) then
        it0=it
      endif
    enddo  
    nt_step_P(it0)=nt_step_P(it0)+1
  endif ! inside the histogram
enddo
open(iunit,file='STAT_histo_P'//'.'//sta_name(ista))
do it=1,nt_hist_P
write(iunit,*) dt_step_P(it)+0.5*dt_hist_P,nt_step_P(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "STAT_histo_res_P'//'.'//sta_name(ista)//'.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"STAT_histo_P'//'.'//sta_name(ista)//'" linecolor 1 with boxes'
write(9,'(a)') 'EOD'
endif

!*************************************** S
if(chx == 1 .or. chx == 2) then
!
!  make the histogram in S  for one station
!
nt_step_S(:)=0
icounts=0
do idata=ntp+1,ntp+nts
  if(res_sta(idata,ista) >= dtmin_S .and. res_sta(idata,ista) <= dtmax_S) then
    icounts=icounts+1
    it0=1
    do it=2,nt_hist_S
      if(res_sta(idata,ista) > dt_step_S(it)) then
        it0=it
      endif
    enddo  
    nt_step_S(it0)=nt_step_S(it0)+1
  endif ! inside the histogram
enddo
open(iunit,file='STAT_histo_S'//'.'//sta_name(ista))
do it=1,nt_hist_S
write(iunit,*) dt_step_S(it)+0.5*dt_hist_S,nt_step_S(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "STAT_histo_res_S'//'.'//sta_name(ista)//'.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"STAT_histo_S'//'.'//sta_name(ista)//'" linecolor 3 with boxes'
write(9,'(a)') 'EOD'
endif


write(*,*) ' for station ',sta_name(ista),' with id', sta_id(ista)
write(*,*) ' number of picked P ',icountp, 'number of picked S ',icounts

enddo ! over station ista

close(9)
stop
end program residue_analysis

