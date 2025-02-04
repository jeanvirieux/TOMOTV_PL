! ###############################################################################
!                                                                               |
!       Compute slowness and hypocentral derivative travel times                |
!                                                                               |
!                                                                               |
!   input (ascii)                                                               |
!                                                                               |
!       inversion.par                                                           |
!                                                                               |
!                                                                               |
!   input (binary)                                                              |
!                                                                               |
!       fsrc        (sources)    8*4 = 32 bytes for each record                 |
! BINARY FORMAT     nsrc,neqks,nshot,nblast,i4,i5,i6,i7                         |
!   record 1+irec   sr_x(irec),sr_y(irec),sr_z(irec),sr_to(irec),   &           |
!                   ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)             |
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
!       modelP      (P velocity)                                                |
!       [modelS]    (S velocity)                                                |
!                                                                               !
!       frays.coor  (rays coordinates)                                          | 
!       frays.tabl  (rays pointers)                                             |
!                                                                               |
!   output (ascii)                                                              |
!                                                                               |
!       matrice.par                                                             |
!       flog.derive                                                             |
!                                                                               |
!   output (binary)                                                             |
!                                                                               |
!       fmat.x     (derivative matrix)                                          |
!       fmat.ic    (index colunmn matrix)                                       |
!       fmat.id    (pointer line matrix)                                        |
!       fmat.pt    (pointer line matrix)                  OBSOLETE              |
!       fdift      (residuals)                                                  |
!       fcal       (calculated travel times)                                    |
!                                                                               |
!      -----------------input files description------------------------         |
!                                                                               |
!---inversion.par---------------------------------------------------------------| 
!                                                                               |
!      chx             (integer) if chx=1=> P waves, if chx=2 => S waves        |
!      xo yo zo        (real) model origin                                      |
!      n1 n2 n3        (integer) nodes number                                   |
!      h1 h2 h3        (real) mesh spacing                                      |
!      modelP          (character) P model file                                 |
!      ixo iyo izo     (integer) inversion grid origin                          |
!      nix niy niz     (integer) inversible nodes numbers                       |
!      [modelS]        (character) if chx=2 S model file                        |
!                      -----------                                              |
!     ixo,iyo,izo             : origin of inversible domain                     |
!     ixo+ni1,iyo+ni2,izo+ni3 : end of inversible domain                        |
!     ni1 < 0 , ni2 < 0 , ni3 <0 => velocity model not inversible               |
!                                                                               | 
!                                                                               |
!                                  MARCH 2010                                   |
!================================================================================
! modif : introduction of the residues for each data (nt values) not used for
!         the inversion but for further analysis
! introduction of output file fresi      which is residues without any weight
! introduction of output file fresu      with weight  
! introduction of output file fidnl
! introduction of input file  fwei      user defined weight
!################################################################################
!                                  AUGUST 2012                                  |
!  backraytracing may fail: some rays are missing : only one fake node stored   |
!  during raytracing and we include here the way to manage it                   |
!  which was not done before although included in the raytracing module         |
!                                                                               |
!  synthetic data should have negative times ....                               |
!================================================================================
! APRIL 2017
! focal mechanism
! ajout d'un test sur dt_static pour la valeur fake -999. si pas d'estimation
!================================================================================
program derivee

  use s_interpol3ds
  use s_read_obs
  use s_read_par
  use s_read_src
  use s_read_fray
  use s_submain
  use s_valu_lent
  use s_store_line

  implicit none
  integer(kind=4) :: nnz,nu,iunit,flog

  ! ray definition
  integer(kind=4) :: iray,rmax
  integer(kind=4) :: nl,nlp,nls                ! number of data with rays
  integer(kind=4) :: nt,ntp,nts                ! number of data
  integer(kind=4) :: nd,ndp,nds,ndw            ! number of data used in the inversion
  real(kind=4) :: l_ray,tps_ray

  ! grid definition  
  integer(kind=4) :: n1,n2,n3,i1,i2,i3,i4,i5,i6,i7
  real(kind=4)    :: h1,h2,h3,xo,yo,zo
  integer(kind=4) :: ixo,iyo,izo,ni1,ni2,ni3,ni11,ni22,ni33
  integer(kind=4) :: nvit,chx

  !     units files 
  integer(kind=4) :: umat,uic,uid,urai,udif,ural,uvp,uvs,ucal,ures,uwes,umeca,uidn
  !     inversion files
  character(len=132) :: fmat_x,fmat_ic,fmat_id

  ! ------------------- sources

  real(kind=4), allocatable,dimension(:) :: sr_x,sr_y,sr_z,sr_to
  real(kind=4) :: sr_to_cur               ! backup it for 'tres' estimation
  integer(kind=4), allocatable,dimension(:) :: ki_pos,ki_to,ki_m,ki_id
  real(kind=4), allocatable,dimension(:) :: plent,slent
  real(kind=4) :: plent_cur,slent_cur
  integer(kind=4) :: ki_pos_cur,ki_to_cur,ki_m_cur
  integer(kind=4) :: idsrc,nsrc,neqks,nshot,nblast
  real(kind=4) :: dtx,dty,dtz            ! slowness vector
  real(kind=4) :: azimuth,dip            ! for focal mecanisms ...

  ! ------------------- stations

  integer(kind=4) :: idsta,ista,nsta,iter,i
  integer(kind=4), allocatable,dimension(:) :: id_sta
  real(kind=4), allocatable,dimension(:) :: dtp_sta,dts_sta
  real(kind=4) :: xdum,ydum,zdum,dt_static
  character(len=4) :: namedum

  ! ------------------- data

  real(kind=4) :: sum_weight
  real(kind=4), allocatable,dimension(:) :: temps,dtemps,weight
  real(kind=4) :: rms,rms_weight,tdat
  integer(kind=4), allocatable,dimension(:) :: lut_src,lut_sta,id_dat,lut_ray
  integer(kind=4) :: idata,noffset,iloop

  ! ------------------- P model or S model

  real(kind=4), allocatable,dimension(:,:,:) :: vit
  character(len=132) :: fvp,fvs

  ! ------------------- inversion matrices for one DATA

  real(kind=4), allocatable,dimension(:,:,:) :: der   ! 3D full for velocity
  !real(kind=4), allocatable,dimension(:) :: der_line  ! 1D full for velocity, positions and origin times 
  !integer(kind=4), allocatable,dimension(:) :: der_ig ! 1D index
  integer(kind=4) :: kptm,kptl,nc,nc_src

  ! ------------------- sparse matrices    for all DATA  (should be cancelled out ... writing it directly)
  !real(kind=4), allocatable,dimension(:) :: mat
  !integer(kind=4), allocatable,dimension(:) :: ic,id

  ! ------------------- rayons

  integer(kind=4), allocatable,dimension(:) :: pt_ray    ! numero du rayon
  integer(kind=4), allocatable,dimension(:,:) :: ind_ray    ! ID couple (SRC,STA) for each ray
  real(kind=4), allocatable,dimension(:,:) :: coord_ray
  !real(kind=4), allocatable,dimension(:) :: ray_length
  integer(kind=4) :: incr,incr_t

  ! ------------------- synthetic data

  real(kind=4), allocatable,dimension(:) :: tcal,tfake

  ! ------------------- statistics to be build-up

  !integer(kind=4), allocatable,dimension(:) :: ihisto
  !      i(28)=maloc(nclass)            ! histogramme
  !real(kind=4), allocatable,dimension(:) :: poids_stat
  !      i(29)=maloc(nl)                ! poids_stat

  iunit=10           ! we start at ten
  flog=iunit         ! log file on derivative
  iunit=iunit+1
  ucal=iunit         ! ucal=11
  iunit=iunit+1
  uvp=iunit          ! uvp=12
  iunit=iunit+1
  uvs=iunit          ! uvs=13   set to zero if only P
  iunit=iunit+1

  open(flog,file='flog.derive')

  ! =================================
  ! ---------------------------------
  !
  !                   READ INPUTS
  !
  ! ---------------------------------
  ! =================================
  write(flog,*) ' reading inversion.par -------------------------'

  call read_par(iunit,chx,fvp,fvs,xo,yo,zo,n1,n2,n3,h1,h2,h3,    &
       ixo,iyo,izo,ni1,ni2,ni3,flog)

  if(chx == 2) then
     write(flog,*) ' both P and S velocity inversion '
  elseif(chx == 1) then
     uvs=0
     write(flog,*) ' only P velocity inversion '
  else
     write(flog,*) ' error in choice P or P&S derive_V stop '
     write(*,*) ' error in choice P or P&S derive_V stop '
     stop
  endif

  ! =================================
  !  allocation de pt_ray et de ind_ray se fait dans la subroutine
  !  ouverture de urai relie a frays.coord et lecture/fermeture de frays.tabl
  ! =================================
  write(flog,*) ' working on rays files ---------------------'

  urai=iunit    ! expected 14
  iunit=iunit+1 ! will be used and freed in read_fray

  write(flog,'(10x,a,i5,a11)') 'open ',iunit,'frays.tabl'
  open(iunit,file='frays.tabl',access= 'direct',recl=3*4)

  read(iunit,rec=1) nl,nlp,nls    ! read number of computed rays split into P and S rays
  ! numbers here for allocating memory stacks for pt_ray and ind_ray
  allocate(pt_ray(nl+1))          ! global index for each ray (2:nl+1) allowing direct-access reading
  allocate(ind_ray(2,nl))         ! indexes of source and station related to this ray NOT YET IDS !!!!

  call read_fray(iunit,urai,pt_ray,ind_ray,rmax,nl,nlp,nls,flog)

  write(flog,*) ' maximum number of points for each ray rmax ',rmax
  write(*,*) ' MAXIMUM number of points for each ray rmax ',rmax

  ! =================================
  !  lecture des sources
  ! =================================
  write(flog,*) ' working on sources files-------------------------------'
  open(iunit,file='fsrc',access='direct',recl=8*4)  ! eight datas of 4 bytes 
  read(iunit,rec=1) nsrc,neqks,nshot,nblast,i4,i5,i6,i7    ! lecture du nbre de seismes et du nbre de blasts

  nc_src=4*neqks+nblast    ! nbre de degres de liberte for sources

  allocate(sr_x(nsrc))      ! x coordinate of source general frame
  allocate(sr_y(nsrc))      ! y
  allocate(sr_z(nsrc))      ! z
  allocate(sr_to(nsrc))     ! t0 origin time
  allocate(ki_pos(nsrc))    ! flag on position
  allocate(ki_to(nsrc))     ! flag on origin time
  allocate(ki_m(nsrc))      ! flag on model
  allocate(ki_id(nsrc))     ! id of the source

  call read_fsrc(iunit,sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id,nsrc)
  close(iunit)

  write(flog,*) 'evaluate P & S hypocenters slowness -----------------------' 
  ! allocation of plent and slent in this subroutine   some compilers do not want it : allocate here
  allocate(plent(nsrc))
  if(uvs /= 0) then
     allocate(slent(nsrc))
  endif
  call valu_lent(uvp,uvs,fvp,fvs,sr_x,sr_y,sr_z,ki_pos,nsrc,xo,yo,zo,n1,n2,n3,h1,h2,h3, &
       plent,slent,flog)

  write(*,*) ' end of the allocation for events '

  ! =================================
  !  lecture des stations  : important for P & S static delay times at stations
  ! =================================

  !
  ! fsta file : format id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name_sta (we only need number of stations)
  !

  open(iunit,file='fsta',access='direct',recl=7*4) ! we open the stream on file fsta
  read(iunit,rec=1) nsta,i1,i2,i3,i4,i5,i6
  allocate(id_sta(nsta))     ! id of the station
  allocate(dtp_sta(nsta))     ! delay time at the station
  allocate(dts_sta(nsta))     ! delay time at the station
  do ista=1,nsta
     read(iunit,rec=1+ista) id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),namedum
  enddo
  close(iunit)

  write(*,*) ' end of the allocation for stations '

  !     workspace
  !if(ni1.gt.0.and.ni2.gt.0.and.ni3.gt.0) then
  !  nvit=ni1*ni2*ni3
  !else
  !  nvit=0
  !  ni1=0
  !  ni2=0
  !  ni3=0
  !endif

  !allocate(ihisto(nclass))
  !allocate(poids_stat(nl))
  !ihisto(:)=0
  !poids_stat(:)=0.0


  ! ================================================================
  ! reading travel time data
  ! ================================================================

  write(flog,*) 'read recorded travel time------------------------'

  ! opening/closing files is done here
  ! fobs : id_dat,temps,dtemps,lut_src,lut_sta,lut_ray
  open(iunit,file='fobs',access='direct',recl=6*4)
  read(iunit,rec=1) nt,ntp,nts,i1,i2,i3   ! read the total number of data (split into P & S data)

  allocate(id_dat(nt))    ! id of the data
  allocate(temps(nt))      ! both P and S waves   ntp first and then nts
  allocate(dtemps(nt))
  allocate(lut_src(nt))    ! id of the source
  allocate(lut_sta(nt))    ! id of the station
  allocate(lut_ray(nt))    ! lut for each data (=0 if no ray or iray if one ray)

  call read_fobs(iunit,id_dat,temps,dtemps,lut_src,lut_sta,lut_ray,nt,ntp,nts)

  write(flog,*) ' total number of data ',nt
  write(flog,*) ' P data ',ntp
  write(flog,*) ' S data ',nts
  close(iunit)

  !========================= reading user weight
  allocate(weight(nt))
  open(iunit,file='fwei',access='direct',recl=4*nt)
  read(iunit,rec=1) weight
  close(iunit)

  write(*,*) ' end of the allocation for data '

  !=================================================================
  ! writing the sparse sensivity matrix per line directly to FILES
  !=================================================================

  fmat_x='fmat.x'
  fmat_ic='fmat.ic'
  fmat_id='fmat.id'

  umat=iunit         ! expected 15
  iunit=iunit+1
  open(umat,file=fmat_x,access='direct',recl=4) !  write per float

  uic=iunit          ! expected 16
  iunit=iunit+1
  open(uic,file=fmat_ic,access='direct',recl=4) ! write per integer

  uid=iunit          ! expected 17
  iunit=iunit+1
  open(uid,file=fmat_id,access='direct',recl=4)   ! write per integer

  udif=iunit         ! expected 18
  iunit=iunit+1
  open(udif,file='fdift',access='direct',recl=4) ! write per float  residues for inversion nd<nl<nt

  ures=iunit         ! expected 19
  iunit=iunit+1
  open(ures,file='fresi',access='direct',recl=4)     ! write per float  residues for each data
  do i=1,nt; write(ures,rec=i) -999.; enddo          ! write a fancy value for initialisation

     uwes=iunit         ! expected 20
     iunit=iunit+1
     open(uwes,file='fresu',access='direct',recl=4)     ! write per float  residues for each data with weight
     do i=1,nt; write(uwes,rec=i) -999.; enddo          ! write a fancy value for initialisation

        uidn=iunit         ! expected 21
        iunit=iunit+1
        open(uidn,file='fidnl',access='direct',recl=4)     ! write integer id_obs for accepted data

        ural=iunit         ! expected 22
        iunit=iunit+1
        open(ural,file='frays.xlen',access='direct',recl=4)     ! write per float   ray length

        umeca=iunit        ! expected 23
        iunit=iunit+1
        open(umeca,file='fmeca',access='direct',recl=10*4)     ! write idata,idsrc,idsta,-dtx,-dty,-dtz,azimut,dip,iup-down,iphase (1-P ou 2-S)
        write(umeca,rec=1) nt,ntp,nts,-999,-999,-999,-999,-999,-999,-999   ! first line: number of data
        do i=2,nt+1
           write(umeca,rec=i) -999,-999,-999,-999,-999,-999,-999,-999,-999,-999
        enddo          ! write a fancy value for initialisation


        write(flog,*) ' STORING DIRECTLY IN FILES SENSITIVITY MATRIX -------------'

        ! =====================   WE ASSUME (4-6) *rmax non-zero elements 
        !                         for the sensitivity matrix for each data/ray
        !nnz=40*rmax*nl+nc_src   ! an over-estimation of the non-zero values of the sensitivity matrix
        !JEAN     write(*,*) ' nnz ',nnz

        !allocate(mat(nnz))
        !allocate(ic(nnz))
        !allocate(id(nl+1))     ! +3 ???   back to +1 ???



        ! ##############################################################
        !------------------- COMPUTE DERIVATIVE MATRIX
        !
        !                      LOOPS OVER DATA
        ! ##############################################################

        rms=0.;rms_weight=0.;sum_weight=0.
        nvit=n1*n2*n3   ! bien sur la maille inverse (grossiere) (pas sous la sous-maille (grossiere))
        ! ----------------------- either P or S ...
        !@@@@@@@@@@@@@@@@@@@@@@@@@@ reduce the memory requirement
        !allocate(der(ni1,ni2,ni3))
        ! partial derivatives with respect to slowness or velocity (full storage) a verifier
        ! on the subgrid where parameters are inverted   (ixo:ni1, iyo:ni2, izo:ni3)
        !@@@@@@@@@@@@@@@@@@@@@@@@@@ on alloue der(ni1+1-ixo,ni2+1-iyo,ni3+1-izo)  SEULEMENT
        ni11=ni1+1-ixo;ni22=ni2+1-iyo;ni33=ni3+1-izo
        ! ----------------------- pour les derivees aux sources
        allocate(der(ni11,ni22,ni33))  ! should be erased at each ray tracing
        ! ----------------------- need the velocity structure
        allocate(vit(n1,n2,n3))

        nc=nc_src                ! start with unknows for quakes and blasts
        nc=nc+nvit               ! number of unknowns including P & quakes/blasts

        if(uvs /= 0) then
           nc=nc+nvit             ! number of unknowns including P & S & quakes/blasts
           nvit=nvit+nvit         ! P and S
        endif

        !allocate(der_line(nc))   ! allocate needed either for P OR S 
        !allocate(der_ig(nc))     ! global index of non-zero elements

        kptl=0                   ! SHOUlD BE SIMILAR TO nd at the end
        kptm=0                   ! increment on non-zero parameters

        allocate(coord_ray(3,rmax))  ! allocation for ray coordinates (using rmax)

        allocate(tcal(nt))         ! allocation for synthetic travel times   ALL DATA
        tcal(:)=-999.0

        noffset=0                  ! offset for P data
        open(uvp,file=fvp,access='direct',recl=4*n1*n2*n3)
        read(uvp,rec=1) vit
        close(uvp)

        nd=0        ! true data at this iteration
        ndp=0       ! for P
        nds=0       ! for S    nd=ndp+nds

        ndw=0       ! counter for data with synthetics but non-zero user- weights
        incr_t= 0   ! on incremente les parametres a modifier

        write(flog,*) ' expected number of unknowns at this iteration ',nc
        write(flog,*) ' unknowns for velocities ',nvit, ' and for location ',nc_src
        write(flog,*) ' expected number of rays  nl',nl
        write(flog,*) ' number of data  nt',nt
        write(flog,*) ' number of data in the inversion not yet estimated=0',nd

        ! ###########################################################
        !
        ! LOOP OVER DATA ..... P DATA FIRST
        !
        ! ###########################################################
        do idata=1,ntp
           if(mod(idata,10000).eq.0) then
              write(*,'(30x,a10,i12,a3,i12)') 'P data ',idata,' / ',ntp
              write(flog,'(30x,a10,i12,a3,i12)') 'P data ',idata,' / ',ntp
           endif
           der(:,:,:)=0.         ! erase la matrice de sensibilite en format plein pour une donnee
           idsrc=lut_src(idata)  ! id of the source     === for position and origin time derivatives 
           idsta=lut_sta(idata)  ! id of the station    === for station statics
           iray=lut_ray(idata)   ! iray associated ray for this data **** VERY IMPORTANT ****
           if(iray /= 0) then
              write(flog,*) ' P idata number and the ray number iray  : ',idata,iray
              write(flog,*) ' ID src and ID src_ray ',idsrc,ind_ray(1,iray)
              write(flog,*) ' ID sta and ID sta_ray ',idsta,ind_ray(2,iray)
           else
              write(flog,*) 'Source outside the domain for P wave iray=0, idsrc, idsta, iray ', idsrc,idsta, iray 
           endif
           !=============================================================
           ! scanning stations and sources for getting the correct one and corresponding values
           !=============================================================
           ! estimate the static delay time at the station through its ID
           iter=0
           do i=1,nsta
              if(idsta == id_sta(i)) then
                 dt_static=dtp_sta(i)
                 if(abs(dt_static+999.) < 1.) dt_static=0.     ! pas de statique connue  JEAN
                 iter=iter+1
              endif
           enddo
           if(iter /= 1) then
              write(flog,*) ' WARNING UNKNOWN STATION ID: missing or too many',idsta
              write(flog,*) ' In the derive_slowness.tomo estimation for P'
              dt_static=0.
           endif
           !     
           ! ............... finding the source related to the idsrc   
           !  for sr_to_cur,ki_pos_cur,ki_to_cur,plent_cur,slent_cur
           !
           iter=0
           do iloop=1,nsrc
              if(idsrc == ki_id(iloop)) then
                 sr_to_cur=sr_to(iloop)
                 plent_cur=plent(iloop)
                 ki_pos_cur=ki_pos(iloop)
                 ki_to_cur=ki_to(iloop)
                 ki_m_cur=ki_m(iloop)
                 iter=iter+1
              endif
           enddo
           if(iter /= 1) then
              write(flog,*) ' ERROR UNKNOWN SOURCE ID',idsrc
              write(flog,*) ' In the derive_slowness.tomo estimation for P'
              stop
           endif
           !============================================================= 
           ! only for data with an associated ray for P for sources inside the model
           !=============================================================
           if(iray /= 0 .and. ki_m_cur /= 0) then
              ! compute slowness derivatives in the 3D structure (as well as estimated travel-times)
              call submain(urai,der,vit,pt_ray,l_ray,tps_ray,coord_ray, &
                   xo,yo,zo,n1,n2,n3,h1,h2,h3,ixo,iyo,izo,ni1,ni2,ni3,incr,iray,flog)
              ! we may end with an error ... back-raytracing has failed
              !$$$$$$$$ JEAN   August 2012 $$$$$
              if(tps_ray <= 1.e-21) then    ! if tps_ray stays at zero ... the ray has only one fake node and it is a failure
                 write(*,*) ' MISSING RAY FOR THIS DATA: a failure should have been mentioned in raytracing ',iray
                 write(*,*) ' we may continue by ignoring this data '
                 tps_ray=-999.
                 !     store computed travel times even fake ones
                 tcal(idata)=tps_ray
                 !     store residual through incremental ndp (observed-time minus orig.-time) temps-sr_to ~ computed time
                 !    tres(ndp)=temps(idata)+dt_static-sr_to_cur-tps_ray
                 tdat=temps(idata)+dt_static-sr_to_cur-tps_ray
                 write(ures,rec=idata) tdat                  ! directly into a file at observation posision    nt
                 write(uwes,rec=idata) tdat*weight(idata)    ! directly into a file at observation position    nt
              else
                 ndp=ndp+1   ! success of derivatives: we may proceed for the inversion itself
                 ! data
                 tcal(idata)=tps_ray  ! save computed travel times
                 !    incremental ndp (observed-time minus orig.-time) temps-sr_to ~ computed time
                 !    tres(ndp)=temps(idata)+dt_static-sr_to_cur-tps_ray
                 tdat=temps(idata)+dt_static-sr_to_cur-tps_ray
                 !frechet derivative
                 incr_t=incr_t+incr
                 !JEAN   WRITE(87,*) ' iray, incr,incr_t ',iray,incr,incr_t
                 ! store matrix row in temporary 1D full vector   
                 call store_line(der,nu,ki_pos_cur,ki_to_cur,plent_cur,coord_ray,n1,n2,n3, &
                      ixo,iyo,izo,ni1,ni2,ni3,dtx,dty,dtz,               &
                      nsrc,neqks,nvit,noffset,iray,umat,uic,uid,kptl,kptm,flog)
                 write(flog,*) ' P:incremented number of parameters to invert nu ', nu
                 !     store computed travel times
                 tcal(idata)=tps_ray
                 !     store residual through incremental ndp (observed-time minus orig.-time) temps-sr_to ~ computed time
                 !    tres(ndp)=temps(idata)+dt_static-sr_to_cur-tps_ray
                 tdat=temps(idata)+dt_static-sr_to_cur-tps_ray
                 ! store row in a file DIRECTLY **************
                 !    call stli(der_line,der_ig,nu,umat,uic,uid,kptm,kptl,flog)
                 !     store ray length
                 !    ray_length(ndp)=l_ray            
                 write(ural,rec=ndp) l_ray                   ! directly into a file                            nd
                 write(udif,rec=ndp) tdat                    ! directly into a file                            nd
                 write(uidn,rec=ndp) id_dat(idata)           ! save the id of the data for tracking it         nd
                 write(ures,rec=idata) tdat                  ! directly into a file at observation position    nt
                 write(uwes,rec=idata) tdat*weight(idata)    ! directly into a file at observation position    nt



!!!JEAN
                 if(tdat > 40.) then  ! peut-etre une erreur d'arrondi de la minute
                    write(98,*) id_dat(idata),temps(idata),tps_ray,tdat,'sta ID',idsta,' 40s '
!                    write(*,*) 'P data id t_obs,t_syn,res > 40 ',id_dat(idata),temps(idata),tps_ray,tdat,'sta ID',idsta
                 endif
                 if(tdat < -40.) then  ! peut-etre une erreur d'arrondi de la minute
                    write(97,*) id_dat(idata),temps(idata),tps_ray,tdat,'sta ID',idsta,' 40s '
!                    write(*,*) 'P data id t_obs,t_syn,res < 40 ',id_dat(idata),temps(idata),tps_ray,tdat,'sta ID',idsta
                 endif
                 if(temps(idata) < 0.) then
                    write(96,*) id_dat(idata),temps(idata),tps_ray,tdat,'sta ID',idsta
!                    write(*,*) 'reminder P data id t_obs<0,t_syn,res ',id_dat(idata),temps(idata),tps_ray,tdat,'sta ID',idsta
                 endif
!!!JEAN







                 !================ JEAN AVRIL 2017  ... binary
                 call slowness2angle(-dtx,-dty,-dtz,azimuth,dip)
                 write(umeca,rec=idata+1) id_dat(idata),idsrc,idsta,-dtx,-dty,-dtz,azimuth,dip,0,1   ! always +1

                 
                 !----------------- make statistics   two rms here
                 rms=rms+tdat**2       ! only data connected to a true ray    NO WEIGHTING
                 if(weight(idata) > 0.) then  ! only count those with non-zero weight
                    ndw=ndw+1
                    sum_weight=sum_weight+weight(idata)
                    rms_weight=rms_weight+(tdat*weight(idata))**2  ! apply user weight   INITIAL
                 endif
              endif           ! end of if on tps_ray
           endif
        enddo

        deallocate(plent)
        deallocate(dtp_sta)

        !==============================================================
        !
        !   LOOP OVER DATA ..... S DATA SECOND ONLY IF uvs # 0
        !
        !============================================================== 
        if(uvs /= 0) then
           noffset=n1*n2*n3  ! offset for S data
           open(uvs,file=fvs,access='direct',recl=4*n1*n2*n3)
           read(uvs,rec=1) vit
           close(uvs)
           !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@à start the loop
           do idata=ntp+1,ntp+nts    ! ntp+nts=nt
              if(mod(idata,10000).eq.0) then
                 write(*,'(30x,a10,i12,a3,i12,a11,i12,a3,i12)') 'S data ',idata-ntp,' / ',nts,' : All data',idata,' / ',ntp+nts
                 write(flog,'(30x,a10,i12,a3,i12,a11,i12,a3,i12)') 'S data ',idata-ntp,' / ',nts,' : All data',idata,' / ',ntp+nts
              endif
              der(:,:,:)=0.          ! erase la matrice de sensibilite en format plein pour une donnee
              idsrc=lut_src(idata)   ! id of the source
              idsta=lut_sta(idata)   ! id of the station
              iray=lut_ray(idata)    ! iray associated ray for this data  **** VERY IMPORTANT ****
              if(iray /= 0) then
                 write(flog,*) ' S idata number and the ray number iray  : ',idata,iray
                 write(flog,*) ' ID src and ID src_ray ',idsrc,ind_ray(1,iray)
                 write(flog,*) ' ID sta and ID sta_ray ',idsta,ind_ray(2,iray)
              else
                 write(flog,*) 'Source outside the domain for S iray=0, idsrc, idsta, iray ', idsrc,idsta, iray 
              endif
              !=============================================================
              ! scanning stations and sources for getting the correct one and corresponding values
              !=============================================================
              ! estimate the static delay time at the station through its ID
              iter=0
              do i=1,nsta
                 if(idsta == id_sta(i)) then
                    dt_static=dts_sta(i)
                    if(abs(dt_static+999.) < 1.) dt_static=0.    ! correction jean
                    iter=iter+1
                 endif
              enddo
              if(iter /= 1) then
                 write(flog,*) ' WARNING UNKNOWN STATION ID: missing or too many',idsta
                 write(flog,*) ' In the derive_slowness.tomo estimation for S'
                 dt_static=0.
              endif
              !     
              ! ............... finding the source related to the idsrc   
              !  for sr_to_cur,ki_pos_cur,ki_to_cur,plent_cur,slent_cur
              !
              iter=0
              do iloop=1,nsrc
                 if(idsrc == ki_id(iloop)) then
                    sr_to_cur=sr_to(iloop)
                    slent_cur=slent(iloop)
                    ki_pos_cur=ki_pos(iloop)
                    ki_to_cur=ki_to(iloop)
                    ki_m_cur=ki_m(iloop)
                    iter=iter+1
                 endif
              enddo
              if(iter /= 1) then
                 write(flog,*) ' ERROR UNKNOWN SOURCE ID',idsrc
                 write(flog,*) ' IN THE derive_S.tomo estimation for S'
                 stop
              endif
              !============================================================= 
              ! only for data with an associated ray for S for sources inside the model
              !=============================================================
              if(iray /=0. .and. ki_m_cur /= 0) then
                 
                 !   compute slowness derivatives 3d vector
                 call submain(urai,der,vit,pt_ray,l_ray,tps_ray,coord_ray, &
                      xo,yo,zo,n1,n2,n3,h1,h2,h3,ixo,iyo,izo,ni1,ni2,ni3,incr,iray,flog)
                 
                 if(tps_ray <= 1.e-21) then
                    write(*,*) ' MISSING RAY FOR THIS DATA: a failure should have been mentioned in raytracing ',iray
                    write(*,*) ' we may continue by ignoring this data '
                    tps_ray=-999.
                    !     store computed travel times even fake ones                                                                                 
                    tcal(idata)=tps_ray
                    !     store residual through incremental ndp (observed-time minus orig.-time) temps-sr_to ~ computed time                        
                    !    tres(ndp)=temps(idata)+dt_static-sr_to_cur-tps_ray                                                                          
                    tdat=temps(idata)+dt_static-sr_to_cur-tps_ray
                    write(ures,rec=idata) tdat                  ! directly into a file at observation posision    nt                                 
                    write(uwes,rec=idata) tdat*weight(idata)    ! directly into a file at observation position    nt  
                 else
                    nds=nds+1   ! success of derivatives: we may proceed for including this data in the inversion procedure
                    !data
                    tcal(idata)=tps_ray   ! store the computed time
                    !   estimate residual through incremental nds (observed-time minus orig.-time)= temps-sr_to ~ computed time
                    !      tres(ndp+nds)=temps(idata)+dt_static-sr_to_cur-tps_ray    array not needed
                    tdat=temps(idata)+dt_static-sr_to_cur-tps_ray
                    !frechet derivative
                    incr_t=incr_t+incr    ! number of nodes where derivatives is non-zero  ... increment
                    !JEAN    WRITE(87,*) ' iray, incr,incr_t ',iray,incr,incr_t
                    ! store matrix data row in temporary 1D full vector  slent(:) a identifier via idsrc
                    call store_line(der,nu,ki_pos_cur,ki_to_cur,slent_cur,coord_ray,n1,n2,n3, &
                         ixo,iyo,izo,ni1,ni2,ni3,dtx,dty,dtz,               &
                         nsrc,neqks,nvit,noffset,iray,umat,uic,uid,kptl,kptm,flog)
                    write(flog,*) ' S:incremented number of parameters to be inverted nu',nu
                    ! store row in files DIRECTLY 
                    !      call stli(der_line,der_ig,nu,umat,uic,uid,kptm,kptl,flog)
                    !     store ray length
                    !      ray_length(ndp+nds)=l_ray
                    write(ural,rec=ndp+nds) l_ray               ! directly into a file                            nd
                    write(udif,rec=ndp+nds) tdat                ! store directly residuals incremental nds        nd
                    write(uidn,rec=ndp+nds) id_dat(idata)       ! save the id of the data for tracking it         nd
                    write(ures,rec=idata) tdat                  ! directly into a file at observation position    nt
                    write(uwes,rec=idata) tdat*weight(idata)    ! directly into a file at observation position    nt



!!!JEAN
                    if(tdat > 40.) then  ! peut-etre une erreur d'arrondi de la minute
                       write(98,*) id_dat(idata),temps(idata),tps_ray,tdat,'ID sta ',idsta, ' 40s '
!                       write(*,*) 'S data id t_obs,t_syn,res > 40 ',id_dat(idata),temps(idata),tps_ray,tdat,'ID sta ',idsta
                    endif
                    if(tdat < -40.) then  ! peut-etre une erreur d'arrondi de la minute
                       write(97,*) id_dat(idata),temps(idata),tps_ray,tdat,'ID sta ',idsta, ' -40s '
!                       write(*,*) 'S data id t_obs,t_syn,res < -40 ',id_dat(idata),temps(idata),tps_ray,tdat,'ID sta ',idsta
                    endif
                    if(temps(idata) < 0.) then
                       write(96,*) id_dat(idata),temps(idata),tps_ray,tdat,'ID sta ',idsta
!                       write(*,*) 'reminder S data id t_obs<0,t_syn,res ',id_dat(idata),temps(idata),tps_ray,tdat,'ID sta ',idsta
                    endif
!!!JEAN





                    !================ JEAN AVRIL 2017 ---> binary
                    call slowness2angle(-dtx,-dty,-dtz,azimuth,dip)
                    write(umeca,rec=idata+1) id_dat(idata),idsrc,idsta,-dtx,-dty,-dtz,azimuth,dip,0,2   ! always +1


                    !-------------- statistics   two rms here
                    rms=rms+tdat**2                             ! without NO WEIGHTING
                    if(weight(idata) > 0.) then                 ! only count non-zero weight
                       sum_weight=sum_weight+weight(idata)
                       rms_weight=rms_weight+(tdat*weight(idata))**2  ! apply user weight
                       ndw=ndw+1
                    endif
                 endif ! end of the if on tps_ray
              endif
           enddo
           deallocate(slent)   ! only allocated when uvs /= 0 and, therefore deallocated inside the if/endif
        endif
        deallocate(dts_sta) ! always allocated ... therefore we can deallocate it

        deallocate(der)     ! 3D
        !deallocate(der_line)! 1D
        !deallocate(der_ig)  ! 1D index

        deallocate(ind_ray)

        nd=ndp+nds

        write(flog,*) 'nbre of rays, nbre of related data and nbre of true data',nl,nd,nt
        write(flog,*) 'corresponding P values ray - associated data - all data',nlp,ndp,ntp
        write(flog,*) 'corresponding S values ray - associated data - all data',nls,nds,nts

        write(*,*) 'nbre of rays, nbre of related data and nbre of true data',nl,nd,nt
        write(*,*) 'corresponding P values ray - associated data - all data',nlp,ndp,ntp
        write(*,*) 'corresponding S values ray - associated data - all data',nls,nds,nts

        !JEAN write(*,*) ' counting maximum of ',8*incr_t,' parameters to be updated : merging is possible'

        if(nd < nl) then
           write(flog,*) ' nd (all true rays) < nl (all rays) < nt (all data) ',nd,nl,nt
           write(flog,*) ' nt true data '
           write(flog,*) ' nl data related to all rays '
           write(flog,*) ' nd data related to true rays in case of raytracing failure '
           write(flog,*) ' at least nd < nl when loosing data in derive_SLOWNESS.tomo '
        endif

        if (nl > nt) then
           write(flog,*) ' something is screwed-up nl (expected data related to a ray)  > nt (all data) ? ',nl,nt
           stop
        endif

        !id(nl+1)=kptm   ! store the number of inverted parameters in the final element of id
        write(uid,rec=nd+1) kptm   ! it is not stored in a vector 
        close(umat)     ! close files for sparse sensitivity matrix
        close(uic)
        close(uid)

        close(udif)    ! RHS .. residues
        close(ures)    ! residues for each observation   ! not a RHS of the linear system
        close(uidn)    ! id of the data
        close(ural)    ! ray length for selected data


        write(flog,*) ' true number of non-zero elements ',kptm
        write(*,*) ' =========================================='
        write(*,*) ' true number of non-zero elements ',kptm
        write(*,*) ' =========================================='

        nnz=kptm           ! update the maximum of non-zero elements

        !============== raw rms   for the nd data with synthetics
        rms=sqrt(rms/float(nd))
        open(49,file='rms.file',access='append')
        write(49,*) rms,' over data with synthetics ',nd
        close(49)
        !============== rms with user weight
        rms_weight=sqrt(rms_weight/sum_weight)
        open(49,file='rms.weig',access='append')
        write(49,*) rms_weight,' over data with non-zero weights ',ndw
        close(49)

        write(*,*) ' **************** '
        write(*,*) ' nd data with synthetics'
        write(*,*) ' ndw data with synthetics and non-zero weights '
        write(*,*) ' nd,ndw ',nd,ndw
        write(*,*) ' **************** '
        write(flog,*) ' **************** '
        write(flog,*) ' nd data with synthetics '
        write(flog,*) ' ndw data with synthetics and non-zero weights '
        write(flog,*) ' nd,ndw ',nd,ndw
        write(flog,*) ' **************** '
        !-------------------
        !
        !-------------------       STORE SPARSE MATRIX IN OUTPUTS FILES

        write(flog,*) ' write parameter file matrice_par ----------------'
        open(49,file='matrix.par')
        write(49,'(a)') fmat_x
        write(49,'(a)') fmat_ic
        write(49,'(a)') fmat_id
        write(49,*) nnz
        write(49,*) nd,ndp
        write(49,*) nc
        close(49)
        !
        write(flog,'(20x,a35,5x,a15)') 'elements non nuls :',fmat_x
        write(flog,'(20x,a35,5x,a15)') 'index des colonnes :',fmat_ic
        write(flog,'(20x,a35,5x,a15)') 'pointeur de lignes :',fmat_id
        write(flog,'(20x,a35,5x,i10)') 'nombre d''elements non nuls :',nnz
        write(flog,'(20x,a35,5x,i7)') 'nombre de lignes :',nd
        write(flog,'(20x,a35,5x,i7)') 'nombre de colonnes :',nc

        write(flog,*) ' write synthetic times -------------------------'
        open(iunit,file='fcal',access='direct',recl=6*4)    ! same format as fobs
        write(iunit,rec=1) nt,ntp,nts,i1,i2    ! write the total number of data (split into P & S data)
        do idata=1,nt 
           write(iunit,rec=1+idata) id_dat(idata),tcal(idata),dtemps(idata),lut_src(idata),lut_sta(idata),lut_ray(idata)
        enddo
        close(iunit)

        deallocate(id_dat)
        deallocate(lut_src)
        deallocate(lut_sta)
        deallocate(lut_ray)
        deallocate(tcal)
        deallocate(dtemps)
        deallocate(temps)
        deallocate(vit)
        deallocate(weight)

        !======================= NO MORE PERTINENT
        !allocate(tfake(nl))
        !write(flog,*) ' write identity matrix for weight for data with synthetics '
        !open(iunit,file='poids',access='direct',recl=4*nl)
        !tfake(:)=1.
        !write(iunit,rec=1) tfake
        !close(iunit)
        !deallocate(tfake)

        !======================= FAKE fnorm identity matrix
        allocate(tfake(nc))
        write(flog,*) ' write identity matrix for normalisation in model space  -------'
        open(iunit,file='fnorm',access='direct',recl=4*nc)
        tfake(:)=1.    ! we take one in case we do not use precond step
        write(iunit,rec=1) tfake
        close(iunit)
        deallocate(tfake)

        !open(iunit,file='derivee.par')    ! not very informative ... should be cancelled
        !write(iunit,'(a)') 'matrix.par'
        !write(iunit,'(a)') 'fdift'
        !write(iunit,*) nvit
        !close(iunit)

        !open(iunit,file='poids_stat',access='direct',recl=4*nl)
        !call sub_weight(tcal,histo,nclass,nl,poids_stat)
        !write(iunit,rec=1) poids_stat
        !close(iunit)

!!! add the saving of fdift before going to precond and inversion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! create fdift.sel ... where are only fdift used for RHS ... nd 2022 JEAN        
        open(udif,file='fdift',access='direct',recl=4)
        open(udif+1,file='fdift.sel',status='unknown')

        write(udif+1,*) ndp,nds,nd,'!',ntp,nts,nt

        do idata=1,nd
           read(udif,rec=idata) rms
           write(udif+1,*) rms
        enddo
        close(udif)
        close(udif+1)

        stop
      end program derivee
 
