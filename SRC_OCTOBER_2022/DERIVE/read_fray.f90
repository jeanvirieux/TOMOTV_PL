MODULE s_read_fray
contains

! ############################################################################
!
!
!     READ FILES : fray.coor and  fray.tab
!
!     input :  
!     iunit   : compteur unite logique fichiers  
!
!     output :
!     urai         : unite logigue du fichier frays.coor   open here
!     nl           : nombre total de rais
!     nlp          : nombre de rayons P
!     nls          : nombre de rayons S
!     rmax         : nombre max de points par rais!
!     pt_ray(nl+1) : pointeur sur le fichier frays.coor
!     ind_ray(2,nl): source_index and station_index for each ray   (IDS !!)
!
! ##############################################################################

subroutine read_fray(iunit,urai,pt_ray,ind_ray,rmax,nl,nlp,nls,flog)
implicit none
integer(kind=4) :: flog
integer(kind=4) :: k,iunit,rmax,nl,nlp,nls,urai
integer(kind=4), dimension(:,:) :: ind_ray
integer(kind=4), dimension(:) :: pt_ray

if(nl /= nlp+nls) then
  write(*,*) ' inconsistent number of rays in read_fray: could not survive '
  stop
endif

write(flog,'(//35x,a/)') 'read ray index and open ray coordinates file'

!--------- open the coordinate file but do nothing here
write(flog,'(10x,a,i5,a11)') 'open ',urai,'frays.coor'
open(urai,file='frays.coor',access='direct',recl=3*4)

pt_ray(:)=0
ind_ray(:,:)=0

call readt(iunit,pt_ray,ind_ray,nl,k,rmax)   ! read the pointer and the index
!     rmax gives the maximum number of points for the longest ray
close(iunit)   ! we have finished with indexes reading ... now in tables pt_ray & ind_ray

return
end subroutine read_fray

! ########################################################################
! ........................................................................
!     Lecture du fichier frays.tabl
!
!     input : 
!     iunit  : unite logique du fichier frays.tabl
!     nl : nombre total de rais
!
!     output :
!     pt_ray(nl+1)  : pointeur sur fichier frays.coor
!     ind_ray(2,nl) : indexes for source and station for each ray  (IDS?)
!     rmax :  nombre max de points par rais
!
! ........................................................................
! ########################################################################
subroutine readt(iunit,pt_ray,ind_ray,nl,k,rmax)
  implicit none
  integer(kind=4) :: idum1,idum2
  integer(kind=4) :: iunit,nl,k,rmax
  integer(kind=4), dimension(:) :: pt_ray(nl+1),ind_ray(2,nl)

  rmax=0   ! set to zero the max number of points in a ray

  !  k=1
  read(iunit,rec=2) pt_ray(1),ind_ray(1,1),ind_ray(2,1)

  open(96,file='flog.ray_pointer',status='unknown')

  do k=2,nl   ! loop over rays for reading indexes  (not coordinates)

     read(iunit,rec=k+1) pt_ray(k),ind_ray(1,k),ind_ray(2,k)
     if(mod(k,5000) == 0) then
        write(96,*) 'fray.tabl ',k
        write(96,*) ' for ray ',k, ' pointer to the position table',pt_ray(k)
        write(96,*) ' index source and index station (IDS maybe)',ind_ray(1,k),ind_ray(2,k)
     endif
     rmax=max(rmax,pt_ray(k)-pt_ray(k-1)+1)    ! number of points for ray (k-1)
  enddo   ! over k number of rays

  close(96)

  read(iunit,rec=nl+2) pt_ray(nl+1),idum1,idum2
  rmax=max(rmax,pt_ray(nl+1)-pt_ray(nl)+1)  ! number of points for ray (nl)

  return
end subroutine readt

END MODULE s_read_fray

