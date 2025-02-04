!###############################################################################
!   Reading the file of station definition "fsta"
!   input  : usta file unit
!            nsta number of stations
!   output :id_sta(1:nsta) ! id of the station
!           x_sta(1:nsta),y_sta(1:nsta),z_sta(1:nsta) ! position of the station
!           name is ignored uptonow
!
!  NO CHECKING ABOUT THE STATION POSITION WITH RESPECT TO THE GRID
! ##############################################################################
MODULE s_read_sta
contains

  subroutine read_fsta(usta,nsta,id_sta,x_sta,y_sta,z_sta, &
       box_xmin,box_xmax,box_ymin,box_ymax,box_zmin,box_zmax)
    implicit none

    real(kind=4) :: x_sta(:),y_sta(:),z_sta(:)
    integer(kind=4) :: id_sta(:)
    integer(kind=4) :: nsta,usta,irec,irec1
    character(len=4) :: name

    real(kind=4) :: box_xmin,box_xmax
    real(kind=4) :: box_ymin,box_ymax
    real(kind=4) :: box_zmin,box_zmax

    do irec=1,nsta
       read(usta,rec=irec+1) id_sta(irec),x_sta(irec),y_sta(irec),z_sta(irec),name
    enddo
    !
    !   check unicity of IDs : this ID should be unique and as low as possible
    !

    do irec=1,nsta
       do irec1=irec+1,nsta
          if(id_sta(irec1) == id_sta(irec)) then
             write(*,*) ' FATAL ERROR IN fsta : same id for two stations irec,irec1 ',irec,irec1,id_sta(irec)
             stop
          endif
       enddo
       if(x_sta(irec) < box_xmin .or. x_sta(irec) > box_xmax .or. &
            y_sta(irec) < box_ymin .or. y_sta(irec) > box_ymax .or. &
            z_sta(irec) < box_zmin .or. z_sta(irec) > box_zmax) then
       write(*,*) ' FATAL ERROR IN fsta : the station is outside the box',irec,id_sta(irec)
       write(*,*) 'box xmin,xmax and station x',box_xmin,box_xmax,x_sta(irec)
       write(*,*) 'box ymin,ymax and station y',box_ymin,box_ymax,y_sta(irec)
       write(*,*) 'box zmin,zmax and station z',box_zmin,box_zmax,z_sta(irec)
       stop
    endif

 enddo

 return
end subroutine read_fsta

END MODULE s_read_sta









