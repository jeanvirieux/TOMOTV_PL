!###############################################################################|
! PROGRAM fsta2a     fsta >>> fsta_dd and fsrc_dd.asc
!
!       input (binary file)
!       fsta        (stations)   7*4 = 28 bytes for each record                 |
! BINARY FORMAT     nsta,i1,i2,i3,i4,i5,i6                                      |
!   record 1+irec   id_sta(ista),xdum,ydum,zdum,dtp_sta(ista),dts_sta(ista),name|
!                                                                               |
!
!     output
!       fsta_dd and fsrc_dd.asc
!
!
!###############################################################################|
program fsta2dd
  implicit none

  integer(kind=4) :: ista,nsta,id_sta,i1,i2,i3,i4,i5,i6
  real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta
  character(len=4) :: name

  real(kind=4) :: xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=4) :: dist,seuil
  integer(kind=4) :: jcolor,jsta

  real(kind=4), allocatable, dimension(:) :: xsta,ysta,zsta
  integer(kind=4), allocatable, dimension(:) :: icolor,idsta

  xmin=1.e29
  xmax=-1.e29
  ymin=1.e29
  ymax=-1.e29
  zmin=1.e29
  zmax=-1.e29

  write(*,*) ' Neighbooring analysis of the fsta file '
  write(*,*) ' each event must belong to a cluster '

  write(*,*) ' Enter the distance max (in meters usually) for clustering'
  read(*,*) seuil
  
  open(10,file='fsta',access='direct',recl=7*4) ! we open the stream on file fsta
  open(11,file='fsta_dd',access='direct',recl=5*4)
  open(12,file='fsta_dd.asc',status='unknown')
  

  read(10,rec=1) nsta,i1,i2,i3,i4,i5,i6
  write(11,rec=1) nsta,i1,i1,i1,i1
  write(*,*) ' number of stations ',nsta
  
  allocate(xsta(nsta))
  allocate(ysta(nsta))
  allocate(zsta(nsta))
  allocate(idsta(nsta))
  allocate(icolor(nsta))    ! colorisation of events

  do ista=1,nsta
     read(10,rec=1+ista) id_sta,x_sta,y_sta,z_sta,dtp_sta,dts_sta,name
     xsta(ista)=x_sta
     ysta(ista)=y_sta
     zsta(ista)=z_sta
     idsta(ista)=id_sta
     icolor(ista)=-999    ! starting with negative values
     if(xmin >= x_sta) xmin=x_sta
     if(ymin >= y_sta) ymin=y_sta
     if(zmin >= z_sta) zmin=z_sta
     if(xmax <= x_sta) xmax=x_sta
     if(ymax <= y_sta) ymax=y_sta
     if(zmax <= z_sta) zmax=z_sta
  enddo

  write(*,*) ' extension of the box for stations '
  write(*,*) ' xmin,xmax ',xmin,xmax
  write(*,*) ' ymin,ymax ',ymin,ymax
  write(*,*) ' zmin,zmax ',zmin,zmax
  write(*,*) ' ================================ '


    jcolor=0    ! initiate the counting of the color 
  do ista=1,nsta
     if(icolor(ista) < 0) then
        jcolor=jcolor+1
        icolor(ista)=jcolor   ! add a color if not yet set
     endif   
     do jsta=ista+1,nsta
        if(icolor(jsta) < 0) then    ! not yet in the cluster
           dist=sqrt((xsta(jsta)-xsta(ista))**2+(ysta(jsta)-ysta(ista))**2+(ysta(jsta)-ysta(ista))**2)
           if(dist < seuil) then
!!!              write(*,*) ista,jsta,dist,seuil         
              icolor(jsta)=icolor(ista)
           endif
        endif
     enddo
  enddo

  do ista=1,nsta     ! no negative value expected at the end
     if(icolor(ista) < 0) then
        write(*,*) ' fatal error at station ',ista
        stop
     endif
  enddo   

  write(*,*) ' number of colored clusters ',jcolor

  do ista=1,nsta
     write(11,rec=ista+1) idsta(ista),xsta(ista),ysta(ista),zsta(ista),icolor(ista)
     write(12,'(i8,3f25.6,i8)') idsta(ista),xsta(ista),ysta(ista),zsta(ista),icolor(ista)
  enddo
  
  close(10)
  close(11)
  close(12)

  stop
end program fsta2dd

