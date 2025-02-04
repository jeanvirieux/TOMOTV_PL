!##############################################################################
!  Reading of the data file "fobs"
!
!     input : uobs file unit
!           : nt   number of data
!
!     output: id_dat IDs of data
!           : lut_src
!           : lut_sta
!           : lut_ray
!
!     id_dat(nt) : look-up table for each data : provide the data id
!     lut_src(nt) : look-up table for each data : provide the src id
!     lut_sta(nt) : look-up table for each data : provide the sta id
!     lut_ray(nt) : look-up table for each data : provide the ray number
!
!#############################################################################
!  Writing of the data file "fobs" essentially for storing the lut_ray
!  when performing the ray tracing
!#############################################################################
MODULE s_read_write_obs
contains

  subroutine read_fobs(uobs,id_dat,lut_src,lut_sta,lut_ray,nt)
    implicit none

    integer(kind=4) :: uobs
    integer(kind=4) :: nt,irec,irec1

    integer(kind=4) :: id_dat(nt),lut_src(nt),lut_sta(nt),lut_ray(nt)

    real(kind=4) :: t_lu,dt_lu

    do irec=1,nt ! lecture 
       read(uobs,rec=irec+1) id_dat(irec),t_lu,dt_lu,lut_src(irec),lut_sta(irec),lut_ray(irec)
    enddo

    !
    !   check unicity of IDs   quite expensive
    !
    if(nt < 10000) then
       do irec=1,nt
          do irec1=irec+1,nt
             if(id_dat(irec1) == id_dat(irec)) then
                write(*,*) ' FATAL ERROR IN FSRC : same id for two observations irec,irec1 ',irec,irec1,id_dat(irec)
                stop
             endif
          enddo
       enddo
    else
       write(*,*) ' checking the data redundancy will be too expensive for this amount of data ',nt
       write(*,*) ' should be done before '
    endif

    return
  end subroutine read_fobs

subroutine write_fobs(uobs,id_dat,lut_src,lut_sta,lut_ray,nt)
implicit none

integer(kind=4) :: uobs
integer(kind=4) :: nt,irec

integer(kind=4) :: id_dat(nt),lut_src(nt),lut_sta(nt),lut_ray(nt)

real(kind=4) :: t_lu,dt_lu
integer(kind=4) :: id_dat_cur,lut_src_cur,lut_sta_cur,lut_ray_cur

do irec=1,nt ! lecture 
   read(uobs,rec=irec+1) id_dat_cur,t_lu,dt_lu,lut_src_cur,lut_sta_cur,lut_ray_cur
! checking values
   if(id_dat(irec) /= id_dat_cur) then
     write(*,*) ' ERROR RAYTRACING>WRITE_FOBS (READ_FOBS) ID of the data different ',id_dat(irec),id_dat_cur
   endif
   if(lut_src(irec) /= lut_src_cur) then
     write(*,*) ' ERROR RAYTRACING>WRITE_FOBS (READ_FOBS) source LUT different ',lut_src(irec),lut_src_cur
   endif
   if(lut_sta(irec) /= lut_sta_cur) then
     write(*,*) ' ERROR RAYTRACING>WRITE_FOBS (READ_FOBS) ID of the data different ',lut_sta(irec),lut_sta_cur
   endif

  write(uobs,rec=irec+1) id_dat(irec),t_lu,dt_lu,lut_src(irec),lut_sta(irec),lut_ray(irec)
enddo                

return
end subroutine write_fobs

END MODULE s_read_write_obs









