!#######################################################################
!
!   reading input parameters from the file "timefd.par" built in
!                                 the previous module 
!
!   input :   iu   : file unit (set in the calling program)
!             flog : log file (defined in the calling program)
!
!   output : fvp, fvs   : file name 
!            xo,yo,zo   : P&L grid origin
!            n1,n2,n3   : number of nodes of the grid
!            h          : grid step (cubic grid)
!   output              : Vp file name
!            chx=2      : Vs file name
!            rap        : sampling along the ray
!       
!
!#######################################################################
MODULE s_read_par
contains

subroutine read_par(iu,fvp,fvs,xo,yo,zo,chx,n1,n2,n3,h,rmax,rap,flog)
implicit none

integer(kind=4) :: iu,n1,n2,n3,flog,chx,rmax
real(kind=4) xo,yo,zo,h,rap
character(len=*) fvp,fvs
!     read timefd.par
open(iu,file='timefd.par',status='old',err=1000)
read(iu,*) chx
read(iu,*) rmax
read(iu,*) xo,yo,zo
read(iu,*) n1,n2,n3
read(iu,*) h
read(iu,*) rap
read(iu,'(a)') fvp
if(chx == 2) then
  read(iu,'(a)') fvs
endif
close(iu)
!
write(flog,'(/20x,a)') 'read timefd.par'
write(flog,'(10x,a19,3f10.2)') 'Cartesian origin  :',xo,yo,zo
write(flog,'(10x,a19,3i10)')  'number of nodes   :',n1,n2,n3
write(flog,'(10x,a19,f10.2)') 'grid step         :',h
write(flog,'(10x,a19,f10.2)') 'ray sampling      :',rap
return              
1000 continue
write(flog,*) ' error missing file timefd.par'
write(*,*) ' error missing file timefd.par'
stop
end subroutine read_par

END MODULE s_read_par

