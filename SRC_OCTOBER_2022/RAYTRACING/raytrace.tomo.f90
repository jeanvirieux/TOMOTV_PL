!###############################################################################|
!                                                                               |
!                     RAY TRACING                                               |
!                                                                               |
! computing rays through the time tables computed over a discrete grid          |
! using a finite difference solved of the eikonal equation (code of Podvin and  |
! Lecomte).                                                                     |
! travel times are not saved and considered in this module                      |
!                                                                               |
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@|
!   input (ascii files)                                                               |
!                timefd.par                                                     |
!   input (binary files)                                                        |
!                fsrc           (source description)                            | 
!                fsta           (stations description)                          |   
!                fobs           (arrival date)                                  |
!                model.fwd      (forward P velocity field)                      |
!               [ modelS.fwd ]  (forward S velocity field)                      |
!                                                                               |
!   output (binary files)                                                       |
!                frays.coor     (ray coordinates file)                          |
!                frays.tabl     (index scanning the ray file)                   |
!                                                                               |
!   output (ascii files)                                                        |
!                raytrace.par                                                   |
!                flog.raytrace                                                  |
!                                                                               |
!                                                                               |
!       ------------ FILE DESCRIPTION (input) ------------------                |
!                                                                               | 
!---timefd.par (ascii) parametres-----------------------------------------------|
!                                                                               |
!    chx          (integer)   : chx=1 => only P wave, chx=2 => S and P waves    | 
!    xo yo zo     (real)      : origin of the velocity model                    |
!    n1 n2 n3     (integer)   : number of nodes along the three directions      |
!    h            (real)      : spatial step                                    |
!    model.fwd    (character) : forward P velocity file                         |
!    [modelS.fwd] (character) : forward S velocity file                         |
!                                                                               |
!    this file has been created by the forward model construction module        |        
!    "module model.tomo"                                                        |
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   input (binary)                                                              |
!                                                                               |
!---fsrc (binaire rec=32) file describing sources      -------------------------|  
!                                                                               |
!       fsrc        (sources)    8*4 = 32 bytes for each record                 |
! BINARY FORMAT     nsrc,neqks,nshot,nblast,i4,i5,i6,i7                         |
!   record 1+irec   sr_x(irec),sr_y(irec),sr_z(irec),sr_to(irec),   &           |
!                   ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)             |
!                                                                               |
!             +---------------------------------------------------------------+ |
! rec=1     : |nsr   |neqks |nshot |nblast |neqks_out|   0   |   0    |  0    | |
! ...         +------+------+------+-------+---------+-------+--------+-------| |
! rec=i+1   : | x(i) | y(i) | z(i) | to(i) | kp(i)   | kt(i) | in(i)  |id(i)  | |
! ...         +------+------+------+-------+---------+-------+--------+-------| |
! rec=nsr+1 : |x(nsr)|y(nsr)|z(nsr)|to(nsr)| kp(nsr) |kt(nsr)| in(nsr)|id(nsr)| |
!             +------+------+------+-------+---------+-------+--------+-------| |
!                                                                               |
! nsr       integer(kind=4): number of events                                   |
! neqks     integer(kind=4): number of quakes (unknown position and date)       |
! nshot     integer(kind=4): number of shots (known position and date)          |
! nblast    integer(kind=4): number of blasts (known position; unknown data)    |
! neqks_out integer(kind=4): number of quakes outside the box                   |
!                                                                               |
!    x(j),y(j),z(j) real(kind=4): coordinates of the source (j<nsrc+1)          |
!    to(j)          real(kind=4): origin date of the source j                   | 
!    kp(j)       integer(kind=4): index over positions of quakes                |
!                                                                               | 
!           kp is a counter which is not zero only for quakes                   |
!           For these events, this counter is increasing & continuous           |        
!                                                                               |
!    kt(j)     integer(kind=4)  : index over origin date of sources             |
!                                                                               | 
!           kt is a counter which is not zero only for quakes and blasts        |
!           For these events, the counter is increasing & continuous            |
!                                                                               |
!    in(j)     integer(kind=4)  : if in(j)=1 =>  the source j is in the box     |
!                                 if in(j)=0 =>  the source j is outside the box|
!                                                                               |
!    id(j=)    integer(kind=4)  : ID of the source (UNIQUE)                     |
!===============================================================================| 
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!                                                                               |
!---fsta (binaire rec=28) file describing each station  ------------------------| 
!                                                                               |
!              +--------+-------+-------+-------+---------+---------+----------+|
! rec=1      : | nsta   |   0   |   0   |   0   |    0    |   0     |    0     ||
!     ...      +--------+-------+-------+-------+---------+---------+----------+|
! rec=i+1    : | id(i)  |  x(i) |  y(i) |  z(i) | dtp(i)  | dts(i)  | name(i)  || 
!     ...      +--------+-------+-------+-------+---------+---------+----------+| 
! rec=nsta+1 : |id(nsta)|x(nsta)|y(nsta)|z(nsta)|dtp(nsta)|dts(nsta)|name(nsta)||
!              +--------+-------+-------+-------+---------+---------+----------+|
!                                                                               | 
!     id_sta(i)    integer(kind=4): id of the station (UNIQUE)                  !
!     x(i),y(i),z(i) real(kind=4) : coordinates of the station j (j<nsta+1)     |
!     dtp_sta(i)     real(kind=4) : static delay time for P wave                |
!     dts_sta(i)     real(kind=4) : static delay time for S wave                | 
!     name(i) (character*4) name of the station (not used in the code)          |
!                                                                               |
!===============================================================================|
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
!                                                                               |
!     nt  integer(kind=4) : number of data (arrival date)                       |
!     ntp integer(kind=4) : number of P wave data                               |
!     nts integer(kind=4) : number of S wave data                               |
!                                                                               |
!     id(i) integer(kind=4): id of the data                                     |
!     t(i)  real(kind=4)   : arrival dates for P waves if i<ntp+1               |
!                                          or  S waves if nt>i>ntp              |
!     dt(i) real(kind=4)        : uncertainties in the picked time              |
!     lut_src(i) integer(kind=4): ID of the source                              |
!     lut_sta(i) integer(kind=4): ID of the station                             |
!     lut_ray(i) integer(kind=4): index of the ray                              |
!                                                                               |
!---model*.fwd (binaire rec=4*n1*n2*n3) modele de vitesse P ou S----------------|
!                                                                               |
!   model*.fwd contains an array v with a dimension n1*n2*n3                    |
!   the velocity at a given node (i,j,k) is given by v(n2*n1*(k-1)+n1*(j-1)+i). |
!   in other words, the following loop transforms 3 indexes into 1 index        |
!   indx=0                                                                      |
!   do k=1,n3                                                                   |
!       do j=1,n2                                                               |
!          do i=1,n1                                                            |
!              indx=indx+1                                                      |
!              v(indx)=vit(i,j,k)                                               |
!          end do                                                               |
!       end do                                                                  |
!   end do                                                                      |
!                                                                               |
!  In the real life, we store the array vit in the file model*.fwd :            |
!                                                                               | 
!  open(i,file='model*.fwd',access='direct',recl=4*n1*n2*n3)                    |
!  write(i,rec=1) vit                                                           |
!                                                                               |
!                                                                               |
!###############################################################################
program raytrace
  use s_read_par
  use s_hs
  use s_ray
  use s_read_write_obs
  use s_read_src
  use s_read_sta

  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     memory requirement roughly
  ! 
  !     2*(n1*n2*n3) + 2*(src*sta) + 7*src + 3*sta + 3*rmax + 2*nrays
  !
  !     rmax should be tuned by the used (see timefd.par input file)
  !
  !-------------------------------------------------------------------------------

  ! file unit

  integer(kind=4) :: iunit,usrc,uray,utab,uobs,usta,flog
  integer(kind=4) :: ndum1,ndum2,ndum3,ndum4,ndum5,ndum6     ! dummy variables

  ! P&L grid for FMM

  integer(kind=4) :: n1,n2,n3
  real(kind=4) :: xori,yori,zori,h,h1,h2,h3,rap
  real(kind=4) ,allocatable,dimension(:,:,:), target :: vit     ! either P or S
  real(kind=4) ,allocatable,dimension(:,:,:), target :: temps

  character(len=132) fvp,fvs

  ! stations
  integer(kind=4) :: nsta
  integer(kind=4), allocatable,dimension(:), target :: id_sta         ! id of the data 
  real(kind=4) ,allocatable,dimension(:), target :: x_sta,y_sta,z_sta ! coordinates of stations (x,y,z)
  ! we do not consider here the static delay times at each station 

  ! sources

  real(kind=4) ,allocatable,dimension(:), target :: x_src,y_src,z_src,to_src   ! coordinates of sources  (x,y,z,t0)
  integer(kind=4) ,allocatable,dimension(:), target :: ki_pos,ki_to,ki_m,ki_id ! various integers for sources 
  integer(kind=4) :: nsrc,neqks,nshot,nblast,neqks_out

  ! data

  integer(kind=4) :: nt,ntp,nts,noff,nt_cur
  integer(kind=4) ,allocatable,dimension(:), target :: id_dat,lut_src,lut_sta,lut_ray

  ! RAYS

  integer(kind=4) :: pri,nray,nrp,nrs,kpt,rmax,kmax,kmin   ! rays
  integer(kind=4) :: chx

  ! ADD BOX FOR CHECKING ... sources (and stations)

  character(len=1) :: carac
  integer(kind=4) :: n1inv,n2inv,n3inv
  real(kind=4) :: x1inv,x2inv,x3inv
  real(kind=4) :: h1inv,h2inv,h3inv
  real(kind=4) :: box_xmin,box_xmax
  real(kind=4) :: box_ymin,box_ymax
  real(kind=4) :: box_zmin,box_zmax
  
  ! 
  !-------------------------------------------------------------------------------
  iunit=10
  
  open(iunit,file='model.head',status='old')
  read(iunit,'(a)') carac
  read(iunit,*) chx
    ! ####################################  input
  read(iunit,'(a)') carac
  read(iunit,*) x1inv,x2inv,x3inv
  read(iunit,'(a)') carac
  read(iunit,*) n1inv,n2inv,n3inv
  read(iunit,'(a)') carac
  read(iunit,*) h1inv,h2inv,h3inv
  close(iunit)
  box_xmin=x1inv
  box_xmax=x1inv+(n1inv-1)*h1inv
  box_ymin=x2inv
  box_ymax=x2inv+(n2inv-1)*h2inv 
  box_zmin=x3inv
  box_zmax=x3inv+(n3inv-1)*h3inv
  
  !     fichier info
  flog=iunit
  iunit=iunit+1
  open(flog,file='flog.raytrace')

  write(flog,*) 'box xmin,xmax',box_xmin,box_xmax
  write(flog,*) 'box ymin,ymax',box_ymin,box_ymax
  write(flog,*) 'box zmin,zmax',box_zmin,box_zmax


  !     fichier stations : fsta      unite 11
  usta=iunit
  iunit=iunit+1

  !     fichier sources  : fsrc      unite 12
  usrc=iunit
  iunit=iunit+1

  !     fichier observation : fobs   unite 13
  uobs=iunit
  iunit=iunit+1

  ! deux fichiers rayons (l'un pour les coordonnees et l'autre pour les indices
  uray=iunit                       ! unite 14
  iunit=iunit+1
  open(uray,file='frays.coor',access='direct',recl=12)
  utab=iunit                       ! unite 15
  iunit=iunit+1
  open(utab,file='frays.tabl',access='direct',recl=12)

  ! -------------------   read parameters

  call read_par(iunit,fvp,fvs,xori,yori,zori,chx,n1,n2,n3,h,rmax,rap,flog) 
  h1=h;h2=h;h3=h

  ! -------------------   read other informations   stations, sources and data

  open(usta,file='fsta',access='direct',recl=28)
  read(usta,rec=1) nsta,ndum1,ndum2,ndum3,ndum4,ndum5,ndum6
  allocate(id_sta(nsta))
  allocate(x_sta(nsta))
  allocate(y_sta(nsta))
  allocate(z_sta(nsta))

  call read_fsta(usta,nsta,id_sta,x_sta,y_sta,z_sta, &
       box_xmin,box_xmax,box_ymin,box_ymax,box_zmin,box_zmax)
  close(usta)


  write(*,*) ' Reading station file OK '


  open(usrc,file='fsrc',access='direct',recl=32)
  read(usrc,rec=1) nsrc,neqks,nshot,nblast,neqks_out,ndum1,ndum2,ndum3
  allocate(x_src(nsrc))
  allocate(y_src(nsrc))
  allocate(z_src(nsrc))
  allocate(to_src(nsrc))
  allocate(ki_pos(nsrc))
  allocate(ki_to(nsrc))
  allocate(ki_m(nsrc))
  allocate(ki_id(nsrc))

  call read_fsrc(usrc,nsrc,neqks,nshot,nblast,neqks_out, &
       x_src,y_src,z_src,to_src,ki_pos,ki_to,ki_m,ki_id, &
       box_xmin,box_xmax,box_ymin,box_ymax,box_zmin,box_zmax)   !! fermture dans subroutine

  write(*,*) ' Reading source file OK '

!=========================

  open(uobs,file='fobs',access='direct',recl=24)
  read(uobs,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3

  write(*,*) ' number of total data ',nt

  allocate(id_dat(nt))
  allocate(lut_src(nt))
  allocate(lut_sta(nt))
  allocate(lut_ray(nt))

  call read_fobs(uobs,id_dat,lut_src,lut_sta,lut_ray,nt)
  close(uobs)

  write(*,*) ' OBSERVED TIMES TOTAL, P and S '
  write(*,*) 'nombre de points:',nt,ntp,nts

!=================================================
! checking consistency
!=================================================


  if(nsrc.gt.nsta) then
     pri=1  ! priority station for raytracing
  else
     pri=0  ! priority source for raytracing
  endif
  !pri=0    ! JEAN DEBUG
  !                        =========== 
  !======================= RAY TRACING
  !                        ===========
  kpt=0                     ! INDEX GLOBAL DES RAIS
  nray=0                    ! NOMBRE DE RAIS
  kmax=0                    ! NOMBRE MAX DE POINTS SUR UN RAYON
  kmin=100000               ! NOMBRE MIN DE POINTS SUR UN RAYON

  !
  !======================= MEMORY ALLOCATION FOR FMM and RAYS
  !

  allocate(vit(n1,n2,n3))
  allocate(temps(n1,n2,n3))

  !==============
  !   P WAVE RAY TRACING
  !=============
  write(*,*) '-----------------------------------'
  write(*,'(/30x,a/)') 'P WAVES RAY TRACING'
  write(*,*) '-----------------------------------'
  !   READING P VELOCITY FILE
  open(iunit,file=fvp,access='direct',recl=4*n1*n2*n3)
  read(iunit,rec=1) vit
  close(iunit)

  call sub_hs(vit,n1,n2,n3,h)    ! conversion en lenteur*distance   (1/sec)
  ! computing rays for P waves
  noff=0    ! considering P waves only
  nt_cur=ntp! in submain

  write(*,*) ' P RAY TRACING '

  open(38,file='ray_track',status='unknown')

  call submain(pri,x_src,y_src,z_src,ki_m,ki_id,nsrc,x_sta,y_sta,z_sta,id_sta,nsta, &
       xori,yori,zori,n1,n2,n3,h1,h2,h3,rap,            &
       lut_src,lut_sta,lut_ray,               &
       vit,temps,utab,uray,kpt,nray,rmax,kmax,kmin,noff,nt_cur)
  write(*,*) ' END OF P RAY TRACING '
  nrp=nray

  if(chx ==2) then
     !==============
     !   S WAVE RAY TRACING
     !=============
     write(*,*) '-----------------------------------'
     write(*,'(/30x,a/)') 'S WAVES RAY TRACING'
     write(*,*) '-----------------------------------'
     !   READING S VELOCITY FILE
     open(iunit,file=fvs,access='direct',recl=4*n1*n2*n3)
     read(iunit,rec=1) vit
     close(iunit)

     call sub_hs(vit,n1,n2,n3,h)  ! conversion en lenteur*distance (1/sec)
     ! computing rays for S waves
     noff=ntp
     nt_cur=nts

     write(*,*) ' S RAY TRACING '

     call submain(pri,x_src,y_src,z_src,ki_m,ki_id,nsrc,x_sta,y_sta,z_sta,id_sta,nsta, &
          xori,yori,zori,n1,n2,n3,h1,h2,h3,rap,            &
          lut_src,lut_sta,lut_ray,               &
          vit,temps,utab,uray,kpt,nray,rmax,kmax,kmin,noff,nt_cur)
     write(*,*) ' END OF S RAY TRACING '
  endif

  close(38)
  !
  nrs=nray-nrp
  write(utab,rec=1) nray,nrp,nrs
  write(utab,rec=nray+2) kpt+1,0,0  ! store the final number of points 
  ! for estimating number of points of the last ray 

  write(*,*) ' final number of rays ', nray
  write(*,*) ' number of rays for P waves ',nrp
  write(*,*) ' number of rays for S waves ',nrs
  write(*,*) ' total number of points for all rays ',kpt+1


  !    one should write down the lut_ray inside the fobs file
  !    keep the first record as it is
  open(uobs,file='fobs',access='direct',recl=24)
  call write_fobs(uobs,id_dat,lut_src,lut_sta,lut_ray,nt)
  close(uobs)

  !           outputs    is it interesting ?
  open(49,file='raytrace.par')
  write(49,*) nrp,nrs
  write(49,*) nsrc,nsta
  close(49)
  !
  !
  !
  write(*,*) ' **************************** '
  write(*,*) ' Maximum number of pts on a ray ',kmax
  write(*,*) ' Minimum number of pts on a ray ',kmin
  write(*,*) ' should be used for tuning rmax input parameter in timefd.par '
  write(*,*) ' should be used also to modify the input parameter in stdin '
  write(*,*) ' especially of kmin is small (< 10 for example) '
  write(*,*) ' **************************** '

  stop
end program raytrace

