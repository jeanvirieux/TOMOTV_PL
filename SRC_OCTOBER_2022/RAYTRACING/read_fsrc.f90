! ###############################################################################     
!    Readinf the source file "fsrc"
!
!     input : usrc file unit
!           : nsrc number of sources
!           : neqks number of quakes
!           : nshot number of shots (perfectly known in space and in time)
!           : nblast number of blasts (unknown origin time; known position)
!           : neqks_out quakes outside the box
!
!     output: x_src(;),y_src(:),z_src(:) positions of sources
!             to_src(:)                  origin times of sources
!             ki_pos(:)                  counter for inversion of positions (only quakes)
!             ki_to(:)                   counter for inversion of origin times (only quakes & blasts)
!             ki_m(:)                    when non zero means source in the box
!             ki_id(:)                   IDs of sources UNIQUE
!
! add the box for checking 2022
! ##############################################################################
MODULE s_read_src
contains

  subroutine read_fsrc(usrc,nsrc,neqks,nshot,nblast,neqks_out, &
       x_src,y_src,z_src,to_src,ki_pos,ki_to,ki_m,ki_id, &
       box_xmin,box_xmax,box_ymin,box_ymax,box_zmin,box_zmax)

    implicit none

    integer(kind=4) :: nsrc,irec,irec1,usrc
    integer(kind=4) :: ki_pos(:),ki_to(:),ki_m(:),ki_id(:)
    real(kind=4) :: x_src(:),y_src(:),z_src(:),to_src(:)

    integer(kind=4) :: neqks,nshot,nblast,neqks_out,kpos,kto
    integer(kind=4) :: neqk,neqk_out
    integer(kind=4) :: ndum1=0,ndum2=0,ndum3=0
    
    real(kind=4) :: box_xmin,box_xmax
    real(kind=4) :: box_ymin,box_ymax
    real(kind=4) :: box_zmin,box_zmax


    do irec=1,nsrc
       read(usrc,rec=irec+1) x_src(irec),y_src(irec),z_src(irec),to_src(irec), &
            ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)
    enddo

    !
    !   check unicity of IDs
    !

    do irec=1,nsrc
       do irec1=irec+1,nsrc
          if(ki_id(irec1) == ki_id(irec)) then
             write(*,*) ' FATAL ERROR IN fsrc : same id for two events irec,irec1 ',irec,irec1,ki_id(irec)
             stop
          endif
       enddo
    enddo

    !   we also check the coherence of this file (quite central in the inversion)

    neqk=0
    neqk_out=0
    kpos=0 !STEPHANIE pointeur dynamique kpos
    kto=0  !STEPHANIE pointeur dynamique kto
    do irec=1,neqks+neqks_out
       if(ki_m(irec) /= 0) then
          ! make sure that the event is inside the box
          if(x_src(irec) > box_xmin .and. x_src(irec) < box_xmax .and. &
               y_src(irec) > box_ymin .and. y_src(irec) < box_ymax .and. &
               z_src(irec) > box_zmin .and. z_src(irec) < box_zmax)  then
             neqk=neqk+1
             kpos=kpos+1 !STEPHANIE pointeur dynamique kpos
             kto=kto+1   !STEPHANIE pointeur dynamique kto
             ki_pos(irec)=kpos !STEPHANIE pointeur dynamique kpos
             ki_to(irec)=kto   !STEPHANIE pointeur dynamique kto
          else ! event outside the box and not declared ...
             neqk_out=neqk_out+1
             ki_m(irec)=0                              ! on force the event to be outside the box
             if(ki_pos(irec) /= 0) ki_pos(irec)=0      ! on force a zero
             if(ki_to(irec) /=0) ki_to(irec)=0         ! on force a zero   et on continue
             write(*,*) ' FORCING quake ID',ki_id(irec),' TO BE OUT OF THE BOX !'
          endif
       else ! the event is known to be outside
          write(*,*) ' WARNING quake ID',ki_id(irec),' IS OUT OF THE BOX !'
          neqk_out=neqk_out+1
          if(ki_pos(irec) /= 0) ki_pos(irec)=0      ! on force a zero
          if(ki_to(irec) /=0) ki_to(irec)=0         ! on force a zero   et on continue
       endif
    enddo

    write(*,*) ' VERIFICATION EQKS '
    write(*,*) ' Declared EQKS ',neqks,' found EQKS ',neqk
    write(*,*) ' Declared EQKS OUT ',neqks_out,' found EQKS OUT ',neqk_out

    do irec=neqks+neqks_out+1,neqks+neqks_out+nshot
       if(ki_m(irec) == 0) then
          write(*,*) ' WARNING shot ID',ki_id(irec),' IS OUT OF THE BOX !'
       endif
       if(ki_pos(irec) /= 0) ki_pos(irec)=0     ! shot should have same position
       if(ki_to(irec) /= 0) ki_to(irec)=0       ! shot should have same t0
    enddo

    do irec=neqks+neqks_out+nshot+1,neqks+neqks_out+nshot+nblast
       if(ki_m(irec) == 0) then
          write(*,*) ' WARNING blast ID',ki_id(irec),' IS OUT OF THE BOX !'
       endif
       if(ki_pos(irec) /= 0) ki_pos(irec)=0     ! blast should have same position
    enddo

!!! STEPHANIE  .... fsrc is modified !
!!! ecriture du fichier fsrc avec allocation dynamique de ki_pos et ki_to
    close(usrc)
    open(usrc,file='fsrc',access='direct',recl=32)
    do irec=1,nsrc
       write(usrc,rec=irec+1) x_src(irec),y_src(irec),z_src(irec),to_src(irec), &
            ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)
    enddo
    write(usrc,rec=1) nsrc,neqk,nshot,nblast,neqk_out,ndum1,ndum2,ndum3
    close(usrc)

    return
  end subroutine read_fsrc

END MODULE s_read_src




