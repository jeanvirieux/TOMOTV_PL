!=============================================================
!
!=============================================================
program convert

implicit none
real(kind=4) :: vel
integer(kind=4) :: icount,ivel
character(len=1) :: carac
real(kind=4) :: xorg,yorg,zorg
real(kind=4) :: n1,n2,n3


open(7,file='modelP.dat',status='old',err=20)
open(8,file='modelP.ini',access='direct',recl=4)

!=============================================== we need to be more robust against input errors
open(9,file='modelP.head',status='old',err=30)
read(9,'(a)') carac
read(9,*) xorg,yorg,zorg
read(9,'(a)') carac
read(9,*) n1,n2,n3
close(9)
!===============================================


icount=0

2000 continue
!read(7,*,end=1000) idump,jdump,kdump,ivel
read(7,*,end=1000) ivel
vel=float(ivel)
if(vel < 1.) then
  write(*,*) ' unexpected small value for velocity ',vel
  write(*,*) ' check your input file modelP.dat '
stop
endif
if(vel > 10000.) then
  write(*,*) ' unexpected high value for velocity ',vel
  write(*,*) ' check your input file modelP.dat '
stop
endif
!=========================== ok decent value: we add it  
icount=icount+1
write(8,rec=icount) vel
goto 2000
1000 continue
!=========================== end of the reading
write(*,*) ' reading velocity file: number of values ',icount
!=====================
if(icount == 0) then
  write(*,*) ' error in reading the ascii velocity file '
  stop
endif
!=====================
if(icount < n1*n2*n3) then
  write(*,*) ' not enough data in the file modelP.dat '
  write(*,*) ' check your input file '
  stop
endif
if(icount > n1*n2*n3) then
  write(*,*) ' too many data in the file modelP.dat '
  write(*,*) ' check your input file '
  stop
endif
stop
!=====================
20 continue
write(*,*) ' missing the file modelP.dat'
stop
!=====================
30 continue
write(*,*) ' missing the file modelP.head '
stop
end program convert
