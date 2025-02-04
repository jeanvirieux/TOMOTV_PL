!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! test the conversion date to utime   (integer should give down to second
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
program utime2date
implicit none
integer(kind=4),dimension(6) :: idate
integer(kind=4) :: itime,id
real(kind=8) :: utime
open(7,file='fdate.txt')


100 continue
read(7,*,end=200) id,itime,utime
call unix2date(itime,idate)
write(*,*)  id,idate
goto 100
200 continue
close(7)
stop
end program utime2date
subroutine unix2date(utime, idate)
      implicit none
      integer utime, idate(6)
!utime  input  Unix system time, seconds since 1970.0
!idate  output Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
!-Author  Clive Page, Leicester University, UK.   1995-MAY-2
      integer mjday, nsecs
      real day
!Note the MJD algorithm only works from years 1901 to 2099.
      mjday    = int(utime/86400 + 40587)
      idate(1) = 1858 + int( (mjday + 321.51) / 365.25)
      day      = aint( mod(mjday + 262.25, 365.25) ) + 0.5
      idate(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
      idate(3) = 1 + int(mod(day,30.6))
      nsecs    = mod(utime, 86400)
      idate(6) = mod(nsecs, 60)
      nsecs    = nsecs / 60
      idate(5) = mod(nsecs, 60)
      idate(4) = nsecs / 60
return
end subroutine unix2date
