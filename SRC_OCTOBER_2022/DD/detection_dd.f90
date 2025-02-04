!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! manipulation of files for creating a file detect_dd: each active data has the
!!!                    master event if its data is found
!!!
!!! does not depend on the structure of the system to be solved: search along rows in fact
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program detection_dd
  implicit none

  integer(kind=4) :: nnz,nl,nc,nlp

  real(kind=4), allocatable, dimension(:) :: temps,dtemps
  integer(kind=4), allocatable, dimension(:) :: id_dat,lut_src,lut_sta,lut_ray
  integer(kind=4), allocatable, dimension(:) :: idnl

  integer(kind=4) :: iunit
  integer(kind=4) :: nt,ntp,nts,i1,i2,i3
  
  real(kind=4), allocatable, dimension(:) :: xsrc,ysrc,zsrc
  integer(kind=4), allocatable, dimension(:) :: idsrc,icolor
  integer(kind=4) :: nsrc,isrc


  integer(kind=4), allocatable, dimension(:) :: idmaster
  integer(kind=4) :: ncolor

  integer(kind=4) :: flog,idump

  character(len=132) :: fmat,fic,fid

  iunit=10
  flog=iunit
  iunit=iunit+1

  open(flog,file='flog.detection_dd',status='unknown')

  write(flog,*) '---------------------------------------------------------------'
  write(flog,'(35x,a,a15)') ' reading matrix.par  '
  write(flog,*) '---------------------------------------------------------------'

  open(iunit,file='matrix.par')
  read(iunit,'(a)') fmat
  read(iunit,'(a)') fic
  read(iunit,'(a)') fid
  read(iunit,*) nnz
  read(iunit,*) nl,nlp                    ! nlp ne sert a rien ici
  read(iunit,*) nc                        ! column
  close(iunit)

  write(*,*) ' nbre of rows nl', nl
  write(*,*) ' nbre of P data nlp',nlp
  write(*,*) ' nbre of S data ',nl-nlp
  write(*,*) ' nbre of columns nc',nc
  write(*,*) ' nbre of non-zero elements nnz',nnz

  write(flog,*) ' nbre of rows nl', nl
  write(flog,*) ' nbre of P data nlp',nlp
  write(flog,*) ' nbre of S data ',nl-nlp
  write(flog,*) ' nbre of columns nc',nc
  write(flog,*) ' nbre of non-zero elements nnz',nnz

!!!====================================================
!!!   read fidnl (dimension nl) which provides the id_dat=idnl(*) of the related data
!!!====================================================
  allocate(idnl(nl))
  open(iunit,file='fidnl',access='direct',recl=4*nl)
  read(iunit,rec=1) idnl
  close(iunit)

!!!=========================================
!!!   read fobs for tracking used data       we need nt and id_dat
!!!=========================================
  open(iunit,file='fobs',access='direct',recl=6*4)
  read(iunit,rec=1) nt,ntp,nts,i1,i2,i3   ! read the total number of data (split into P & S data)
  allocate(id_dat(nt))    ! id of the data
  allocate(temps(nt))      ! both P and S waves   ntp first and then nts
  allocate(dtemps(nt))
  allocate(lut_src(nt))    ! id of the source
  allocate(lut_sta(nt))    ! id of the station
  allocate(lut_ray(nt))    ! lut for each data (=0 if no ray or iray if one ray)
  call read_fobs(iunit,id_dat,temps,dtemps,lut_src,lut_sta,lut_ray,nt,ntp,nts)
  close(iunit)
  deallocate(temps,dtemps,lut_ray)   ! not needed 

!!!=========================================
!!!   read fsrc_dd for cluster definition with color ID
!!!=========================================
  open(iunit,file='fsrc_dd',access='direct',recl=5*4)
  read(iunit,rec=1) nsrc,idump,idump,idump,idump
  allocate(idsrc(nsrc))
  allocate(xsrc(nsrc))
  allocate(ysrc(nsrc))
  allocate(zsrc(nsrc))
  allocate(icolor(nsrc)) 
  do isrc=1,nsrc
     read(iunit,rec=isrc+1) idsrc(isrc),xsrc(isrc),ysrc(isrc),zsrc(isrc),icolor(isrc)
  enddo
  close(iunit)

  ncolor=0
  do isrc=1,nsrc
     write(*,*) isrc,icolor(isrc)
     if(icolor(isrc) > ncolor) ncolor=icolor(isrc)
  enddo

  write(*,*) ' number of clusters ',ncolor
  write(flog,*) ' number of clusters ',ncolor

!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ define the master event  
  allocate(idmaster(ncolor))
  call bary_dd(nsrc,idsrc,xsrc,ysrc,zsrc,icolor,ncolor,idmaster)

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% build up the file   detect_dd.asc (ascii for the moment)

  call detect_dd(iunit,nl,nlp,nc,idnl, nt,id_dat,lut_src,lut_sta, nsrc,idsrc,icolor, ncolor,idmaster)

  deallocate(idmaster)
  deallocate(idsrc,icolor)
  deallocate(id_dat,lut_src,lut_sta)
  deallocate(idnl)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the master event for each cluster (barycenter strategy)
!!!           other strategies are possible (for example, number of data)
!!!           which should be elaborated later on  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine bary_dd(nsrc,idsrc,xsrc,ysrc,zsrc,icolor,ncolor,idmaster)
    implicit none
    real(kind=4),allocatable,dimension(:) :: xsrc,ysrc,zsrc
    integer(kind=4), dimension(:) :: idsrc,icolor
    integer(kind=4) :: nsrc,isrc

    integer(kind=4), dimension(:) :: idmaster
    integer(kind=4) :: ncolor

    integer(kind=4) ::  icol,ibary
    real(kind=4) :: xbary,ybary,zbary,dist,seuil

    do icol=1,ncolor

       xbary=0.; ybary=0.;zbary=0.
       do isrc=1,nsrc
          xbary=xbary+xsrc(isrc)
          ybary=ybary+ysrc(isrc)
          zbary=zbary+zsrc(isrc)
       enddo
       xbary=xbary/float(nsrc)
       ybary=ybary/float(nsrc)
       zbary=zbary/float(nsrc)
       seuil=1.e+29
       ibary=1
       do isrc=1,nsrc
          dist=sqrt((xsrc(isrc)-xbary)**2+(ysrc(isrc)-ybary)**2+(zsrc(isrc)-zbary)**2)
          if(dist < seuil) then
             ibary=isrc
             seuil=dist
          endif   
       enddo
       idmaster(icol)=idsrc(ibary)   ! get the ID of the master event for this cluster icol
       
    enddo   

  end subroutine bary_dd

  subroutine detect_dd(iunit,nl,nlp,nc,idnl, nt,id_dat,lut_src,lut_sta, nsrc,idsrc,icolor, ncolor,idmaster)

!!!§§§§§§§§§§§§§§§§§ mat,ic,id,nnz,nl,nc    system A
!!!::::::::::::::::: res,idnl               vector b
!!!&&&&&&&&&&&&&&&&& id_dat,lut_src,lut_sta  ID data connection
!!!<<<<<<<<<<<<<<<<  idsrc,icolor           cluster of events
!!!================  ncolor,idmaster        master event of each cluster

    implicit none

    integer(kind=4), dimension(:) :: idnl
    integer(kind=4) :: iunit,nl,nlp,nc

    integer(kind=4), dimension(:) :: id_dat,lut_src,lut_sta
    integer(kind=4) :: nt

    integer(kind=4), dimension(:) :: idsrc,icolor
    integer(kind=4) :: nsrc

    integer(kind=4), dimension(:) :: idmaster
    integer(kind=4) :: ncolor

    integer(kind=4), allocatable, dimension(:) :: lutsrc,lutsta,irefsrc

    integer(kind=4) :: irow,irow1,it,isrc,ioff,idouble,icluster,irec

    allocate(lutsrc(nl))                    !  ID of the event
    allocate(lutsta(nl))                    !  ID of the station
    allocate(irefsrc(nl))                   !  ID of the master event

    open(iunit,file='detect_dd',access='direct',status='unknown',recl=7*4)
    open(iunit+1,file='detect_dd.asc',status='unknown')
    irec=0
!!!**************************  P *****************************    
!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ search the master event for each data
!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ with or without an available data    
    do irow=1,nlp                                   ! looping over active data
       if(mod(irow,1000) == 0) write(*,*) ' P data search master ',irow
       do it=1,nt                                   ! looping over the database
          if(idnl(irow) == id_dat(it)) then         ! ok we find the data in the database
             lutsrc(irow)=lut_src(it)               ! this is the event of this data
             lutsta(irow)=lut_sta(it)               ! this is the station of this data
             do isrc=1,nsrc                         ! looping over sources
                if(lutsrc(irow) == idsrc(isrc)) then      ! we find the source in the source database
                   icluster=icolor(isrc)            ! get the cluster color
                   irefsrc(irow)=idmaster(icluster) ! we know the master event for this data
                   goto 1000
                endif
             enddo
          endif
       enddo
1000   continue
    enddo !irow

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find if the master event has a data
    do irow=1,nlp                            !looping again over active data
       idouble=0                            ! no detection yet
       if(mod(irow,1000) == 0) write(*,*) ' P data search master data',irow
       do irow1=1,nlp                   !find the row of the master event for the same station
!!! the irow1 comes from the master event (irefsrc)  with a picked data (lutsta)
          if(lutsrc(irow1) ==  irefsrc(irow) .and. lutsta(irow1) == lutsta(irow)) then ! yes
             idouble=1
             irec=irec+1
             write(iunit,rec=irec) irow,irow1,lutsrc(irow),lutsrc(irow1),lutsta(irow),idnl(irow),idnl(irow1)
             write(iunit+1,*) irow,irow1,lutsrc(irow),lutsrc(irow1),lutsta(irow),idnl(irow),idnl(irow1)
             goto 2000
          endif
       enddo  ! end on irow1

2000   continue
    enddo  ! end of irow

    ioff=nlp   !!! offset for S active data
    
    !!!**************************  S *****************************    
!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ search the master event for each data
!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ with or without an available data    
    do irow=ioff+1,nl                               ! looping over active data
       if(mod(irow,1000) == 0) write(*,*) ' S data search master ',irow
       do it=1,nt                                   ! looping over the database
          if(idnl(irow) == id_dat(it)) then         ! ok we find the data in the database
             lutsrc(irow)=lut_src(it)               ! this is the event of this data
             lutsta(irow)=lut_sta(it)               ! this is the station of this data
             do isrc=1,nsrc                         ! looping over sources
                if(lutsrc(irow) == idsrc(isrc)) then      ! we find the source in the source database
                   icluster=icolor(isrc)            ! get the cluster color
                   irefsrc(irow)=idmaster(icluster) ! we know the master event for this data
                   goto 1001
                endif
             enddo
          endif
       enddo
1001   continue
    enddo !irow

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find if the master event has a data
    do irow=ioff+1,nl                            !looping again over active data
       if(mod(irow,1000) == 0) write(*,*) ' S data search master data',irow
       idouble=0                                 ! no detection yet
       do irow1=ioff+1,nl                        !find the row of the master event for the same station
!!! the irow1 comes from the master event (irefsrc)  with a picked data (lutsta)
          if(lutsrc(irow1) ==  irefsrc(irow) .and. lutsta(irow1) == lutsta(irow)) then ! yes
             idouble=1
             irec=irec+1
             write(iunit,rec=irec) irow,irow1,lutsrc(irow),lutsrc(irow1),lutsta(irow),idnl(irow),idnl(irow1)
             write(iunit+1,*) irow,irow1,lutsrc(irow),lutsrc(irow1),lutsta(irow),idnl(irow),idnl(irow1)
             goto 2001
          endif
       enddo  ! end on irow1

2001   continue
    enddo  ! end of irow

    deallocate(lutsrc,lutsta,irefsrc)
    close(iunit); close(iunit+1)

  end subroutine detect_dd

  subroutine read_fobs(uobs,id_dat,temps,dtemps,lut_src,lut_sta,lut_ray,nt,ntp,nts)
    implicit none

    integer(kind=4) :: uobs                  ! id of file
    integer(kind=4) :: nt,ntp,nts,irec       ! total nbre of data

    real(kind=4), dimension(:) :: temps(:),dtemps(:)

    integer(kind=4), dimension(:) :: lut_src(:),lut_sta(:),lut_ray(:)
    integer(kind=4), dimension(:) :: id_dat(:)

    if(nt /= ntp+nts) then
       write(*,*) ' error in read_fobs number of total observed times different from the sum of P and S times'
       stop
    endif

    !==================== reading P wave travel times (ntp first values)
    !reading id of the data as well as the picked time (date in fact)
    !reading luts such that we can find station id, source id and ray id
    !                     this gives us unambiguous relations !!!
    !====================
    do irec=1,ntp ! lecture des temps P
       read(uobs,rec=irec+1) id_dat(irec),temps(irec),dtemps(irec),lut_src(irec),lut_sta(irec),lut_ray(irec)
    enddo
    !==================== reading S wave travel times (ntsp next values) if any
    !build up luts such that we can find data id, station id and source id
    !====================
    if(nts /= 0) then ! lecture des temps S
       do irec=1,nts
          read(uobs,rec=1+ntp+irec) id_dat(ntp+irec),temps(ntp+irec),dtemps(ntp+irec),lut_src(ntp+irec),&
                          lut_sta(ntp+irec),lut_ray(ntp+irec)
       enddo
    endif

    return
  end subroutine read_fobs


end program




















