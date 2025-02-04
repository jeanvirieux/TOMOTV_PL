!###############################################################################|
! PROGRAM facqui2tv     FACQUI  to fsta,fsrc,fobs (format defined by S. Operto) |
! Common Receiver Gather (OBS)
!	  input facqui
!     output (binary file) fsrc,fsta,fobs
!
!  facqui is supposed to be in meters 
!###############################################################################|
program facqui2tv
  implicit none

  integer(kind=4) :: ista,nsta
  integer(kind=4) :: i1,i2,i3,i4,i5,i6
  real(kind=4) :: x_sta,y_sta,z_sta,dtp_sta,dts_sta

  integer(kind=4) :: nsrc,isrc

  character(len=4) :: name
  character(len=3) :: name3
  character(len=132) :: string
  integer(kind=4) :: iempty

  real(kind=4) :: x1,x2,t0,inc,eps,dist,xcode
  integer(kind=4) :: icode,iloop

  integer, parameter :: nsrc_max=100000
  real(kind=4), dimension(nsrc_max) :: x1_src,x2_src
  integer(kind=4), dimension(nsrc_max) :: id_src

  integer, parameter :: nsta_max=1000
  real(kind=4), dimension(nsta_max) :: x1_sta,x2_sta
  integer(kind=4), dimension(nsta_max) :: id_sta

  real(kind=4) :: time_lu,dt_lu
  integer(kind=4) :: iobs,ilu,id_lu,id_ray

  open(11,file='facqui',status='old')

  nsrc=0;i1=0;i2=0;i3=0;i4=0;i5=0;i6=0
  nsta=0
  ! starting with the first non blank line 

  ista=0
  isrc=0
  eps=0.01       ! km    precision 10 meters
  eps=eps*eps   ! take the square

2000 continue
  read(11,*,end=99) x1,x2,t0,inc,xcode

  icode=xcode+0.001

  if(icode == 0) then                  ! here is an OBS ... receiver
     ista=ista+1
     if(ista > nsta_max) then
        write(*,*) ' please increase nsta_max ',nsta_max
        stop
     endif
     write(*,*) ' first pass:  one OBS detected ',ista
     ! we do not check if it is a NEW OBS ... we assume that it is
     x1_sta(ista)=x1;x2_sta(ista)=x2; id_sta(ista)=ista  
     goto 2000 ! it is an OBS  position

  else !   if (icode == 1) then
     !
     ! check if new source
     !

     do iloop=1,isrc
        dist=(x1-x1_src(iloop))**2+(x2-x2_src(iloop))**2
        if(dist <= eps) goto 2000 ! already there
     enddo
     ! new sources
     isrc=isrc+1
     if(isrc > nsrc_max) then
        write(*,*) ' please increase nsrc_max ',nsrc_max
        stop
     endif
     x1_src(isrc)=x1; x2_src(isrc)=x2; id_src(isrc)=isrc   

  endif
  goto 2000

99 continue
  close(11)
  !===================== end of the source and receiver construction

  nsta=ista
  write(*,*) ' number of independent stations',nsta

  nsrc=isrc
  write(*,*) ' number of sources ',nsrc


  !================================================
  ! sources are receivers
  ! receivers are sources
  !================================================

  !=========================== receivers
  write(*,*) ' fsta building ... '
  open(10,file='fsta',access='direct',recl=7*4) ! we open the stream on file fsta
  write(10,rec=1) nsta,i1,i2,i3,i4,i5,i6   ! 


  dtp_sta=0.    ! no p static
  dts_sta=0.    ! no s static

  do ista=1,nsta
     write(name,'(I4.4)') id_sta(ista)
     !======================================= writing fsta file in meters
     write(10,rec=1+ista) id_sta(ista),x2_sta(ista),0.,x1_sta(ista),dtp_sta,dts_sta,name
  enddo
  close(10)

  !=========================== sources (only shots)

  write(*,*) ' fsrc building ... '
  open(10,file='fsrc',access='direct',recl=8*4)

  !====================== quakes
  !========================= shots
  !============================= blasts
  write(10,rec=1) nsrc,0,nsrc,0,0,0,0,0

  do isrc=1,nsrc

     !======================================================kp=0 except for quakes
     !========================================================kt=0 except for quakes and blasts
     !==========================================================km=1 if source inside the domain
     write(10,rec=1+isrc) x2_src(isrc),0.,x1_src(isrc),0.,0,0,1,id_src(isrc)

  enddo
  close(10)

  !============================ building the fobs file

  write(*,*) ' fobs building ... '

  open(11,file='facqui',status='old')

  open(12,file='fobs',access='direct',recl=6*4) ! we open the stream 12 on file fobs

  write(12,rec=1) 0,0,0,0,0,0 ! dummy line upto now for fobs
  id_ray=0   ! no information yet about the ray for this analysis

  id_lu=0
  dt_lu=0.
  ista=0

2001 continue

  read(11,*,end=98) x1,x2,t0,inc,xcode

  icode=xcode+0.001

  if(icode /= 0 ) then
     write(*,*) ' error in reading the facqui should read an OBS icode=0',xcode,icode,x1,x2
     stop
  endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ not very beautiful ...
2010 continue ! branchement spagetti ... 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ but easy ...

  if(icode == 0) then                  ! here is an OBS ... receiver  
    do iloop=1,nsta
        ista=iloop
        dist=(x1-x1_sta(iloop))**2+(x2-x2_sta(iloop))**2
        if(dist <= eps) goto 2002 ! already there
     enddo
     stop 'error in reading receivers '
2002 continue
     write(*,*) ' one OBS detected ',ista
     ilu=0
3001 continue
     ilu=ilu+1       ! local counter for the selected OBS
     id_lu=id_lu+1   ! global counter for the acquisition

     read(11,*,end=98) x1,x2,time_lu,inc,xcode
     icode=xcode+0.001

     if(icode == 0) then
       id_lu=id_lu-1   ! rewind one
       ilu=ilu-1
       goto 2010
     endif

     do iloop=1,nsrc
        isrc=iloop
        dist=(x1-x1_src(iloop))**2+(x2-x2_src(iloop))**2
        if(dist <= eps) goto 3002 ! already there
     enddo
     stop 'error in reading sources '
3002 continue
     write(*,*) ' couple obs, src ', ista,isrc
     write(12,rec=1+id_lu) id_lu,time_lu,dt_lu,id_src(isrc),id_sta(ista),id_ray
     goto 3001
97   continue   ! end of reading the OBS picking time
     close(9)
  endif
  goto 2001     ! next OBS
98 continue
  id_lu=id_lu-1
  write(12,rec=1) id_lu,id_lu,0,0,0,0 
  close(11)
  close(12)

  stop
end program facqui2tv

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

