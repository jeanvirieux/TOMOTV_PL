!##############################################################################
!     Lecture du fichier donnees fobs
!     input : uobs,nt,ntp,nts
!     output: id_dat,temps,lut_src,lut_sta,lut_ray
!
!     nt  : nombre total de temps observes
!     ntp : nombre de temps P
!     nts : nombre de temps S   
!
!     temps(1:ntp) : temps observe P
!     temps(ntp+1:ntp+nts) : temps observe S
!
!
! true binary format : record=1 : nt,ntp,nts,idum,idum    ( 5 * 4 bytes)
!                      record=1+irec :   irec from 1 to nt
!                      id_dat(irec),temps(irec),lut_src(irec),lut_sta(irec),lut_ray(irec)
!
!     id_dat(nt)  : id of the data .... be sure that we do not confuse various data !
!
! LOOK-UP TABLES : tracking sources, stations and rays
!
!     lut_src(nt) : id of the source related to the selected data
!     lut_sta(nt) : id of the station related to the selected data
!     lut_ray(nt) : id of the ray related to the selected data
!                   if ==0 means no rays
!
!#############################################################################
MODULE s_read_obs
contains

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
read(uobs,rec=1+ntp+irec) id_dat(ntp+irec),temps(ntp+irec),dtemps(ntp+irec),lut_src(ntp+irec),lut_sta(ntp+irec),lut_ray(ntp+irec)
  enddo 
endif

return
end subroutine read_fobs

END MODULE s_read_obs









