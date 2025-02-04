!###############################################################################|
! PROGRAM FOR STATISTICS ON STATIONS                                            | 
!                       traveltimes                                             |
!  COMPUTATION OF THE AVERAGE RESIDUE AT EACH STATION AND MAKE NEW fsta FILE    |
!     input                                                                     |
!                                                                               |
!       fobs        (data arrival times) 6*4 = 24 bytes for each record         |
! BINARY FORMAT     nt,ntp,nts,i1,i2,i3                                         |
!   record 1+irec   id_dat(irec),temps(irec),dtemps(irec), &                    |
!                   lut_src(irec),lut_sta(irec),lut_ray(irec)                   |
!                                                                               |
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!                                                                               |
!       fresi       residues for each observation   UNWEIGHTED                  |
!       fresu       residues for each observation   WEIGHTED                    |
!                                                                               |
!                                                                               |
!                                                                               |
!###############################################################################|
program station_res
  use s_read_obs


  implicit none
  character(len=4) :: name

  character(len=4), allocatable, dimension(:), target :: sta_name
  integer(kind=4), allocatable, dimension(:), target :: sta_id
  integer(kind=4), allocatable, dimension(:), target :: sta_num
  real(kind=4), allocatable, dimension(:), target :: sta_dtp,sta_dts

  integer(kind=4) :: nsta,ista
  integer(kind=4), allocatable, dimension(:), target :: id_dat
  integer(kind=4) :: id_sta0
  real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta

  real(kind=4) :: rms,rms_w
  real(kind=4), allocatable, dimension(:), target :: temps,dtemps,res,res_w
  integer(kind=4), allocatable, dimension(:), target :: lut_ray,lut_sta,lut_src
  integer(kind=4) :: nt,ntp,nts,ndp,nds,nd

  integer(kind=4) :: i1,i2,i3,iunit,idata

  open(7,file='out_P',status='unknown')
  open(8,file='out_S',status='unknown')

  open(20,file='residu_P',status='unknown')
  open(21,file='residu_S',status='unknown')

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
     !     write(*,'(a)') sta_name(ista)
     sta_id(ista)=id_sta0
  enddo
  close(10)

  !*************************************************************************
  ! we read now the fobs file for picked time
  ! and         the fresi file for WEIGHTED residues at all observations
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
  allocate(res_w(nt))      ! weighted
  open(iunit,file='fresi',access='direct',recl=4*nt)     ! write per float  residues
  read(iunit,rec=1) res
  close(iunit)
  open(iunit,file='fresu',access='direct',recl=4*nt)     ! write per float  residues
  read(iunit,rec=1) res_w
  close(iunit)

  write(*,*) ' end of reading of both fsta and fobs '

  rms=0.
  rms_w=0.
  ndp=0
  nds=0
  nd=0
  !=============================================================
  !   loop over data   P first
  !=============================================================
  sta_dtp(:)=0.    ! residues at each station is set to zero
  sta_num(:)=0
  do idata=1,ntp
     if(abs(res(idata)+999.) > 1.) then
        if(abs(res(idata)) < 10.) then
           ndp=ndp+1     ! one decent P residue 
           do ista=1,nsta
              if(lut_sta(idata) == sta_id(ista)) then
                 sta_num(ista)=sta_num(ista)+1
                 sta_dtp(ista)=sta_dtp(ista)+res(idata)
              endif
           enddo   ! end on station ista
           rms=rms+res(idata)**2    ! computation of unweighted residuals
           rms_w=rms_w+res_w(idata)**2    ! computation of weighted residuals

        else
           write(20,*) ' too big P residue ',idata,res(idata)
        endif      ! test over residus
        ! else
        !    write(*,*) ' missing P data ',res(idata)
     endif
  enddo ! end on data   in P
  ! make the average for each station and deduce the time static equal to minus ...
  do ista=1,nsta
     if(sta_num(ista) /= 0) then
        sta_dtp(ista)=-sta_dtp(ista)/float(sta_num(ista))   ! static delay at the station
        write(7,*) ista,sta_dtp(ista)
     else
        sta_dtp(ista)=-999.    !!! correction JEAN flag -999. est active dans derive_slowness.tomo
     endif
  enddo

  !=======================================================
  !   loop over S data
  !=======================================================
  sta_dts(:)=0.    ! residues at each station is set to zero
  sta_num(:)=0

  do idata=ntp+1,ntp+nts
     if(abs(res(idata)+999.) > 1.) then   ! if we have estimated a residue for this data

        if(abs(res(idata)) < 10.) then
           nds=nds+1   ! one decent S residu
           do ista=1,nsta
              if(lut_sta(idata) == sta_id(ista)) then
                 sta_num(ista)=sta_num(ista)+1
                 sta_dts(ista)=sta_dts(ista)+res(idata)
              endif
           enddo   ! end on station ista
           rms=rms+res(idata)**2    ! computation of unweighted residuals
           rms_w=rms_w+res_w(idata)**2    ! computation of weighted residuals

        else
           write(21,*) ' too big P residue ',idata,res(idata) 
        endif ! test sur residus
        !  else
        !     write(*,*) ' missing S data ',res(idata)
     endif


  enddo ! end on data   in S
  ! make the average for each station and deduce the time static equal to minus
  do ista=1,nsta
     if(sta_num(ista) /= 0) then
        sta_dts(ista)=-sta_dts(ista)/float(sta_num(ista))   ! static delay at the station
        write(8,*) ista,sta_dts(ista)
     else
        sta_dts(ista)=-999.     !!! correction JEAN flag -999. est active dans derive_slowness.tomo
     endif
  enddo


  !
  !   RMS FOR CHECKING
  !
  nd=ndp+nds    !

  write(*,*) ' total of data ',nt
  write(*,*) ' P data validated ',ndp
  write(*,*) ' S data validated ',nds
  write(*,*) ' total data validated ',nd

  if(nt /= 0 .or. nd /= 0) then
     write(*,*) 'variance ',rms/float(nd)
     rms=sqrt(rms/float(nd))
     write(*,*) ' UNWEIGHTED RMS FOR RESIDUALS (from fresi computed in DERIVE) '
     write(*,*) ' RMS ',rms, ' for all validated data ',nd
     write(*,*) '=================='
     rms_w=sqrt(rms_w/float(nd))
     write(*,*) ' WEIGHTED RMS FOR RESIDUALS (from fresu computed in DERIVE) '
     write(*,*) ' RMS ',rms_w, ' for all validated data  ',nd
     write(*,*) '=================='
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

  close(20)   ! residu_P
  close(21)   ! residu_S
  close(7)    ! out_P
  close(8)    ! out_S
  stop
end program station_res

