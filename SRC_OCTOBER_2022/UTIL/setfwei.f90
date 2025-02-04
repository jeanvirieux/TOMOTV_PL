!###############################################################################
! program setfwei
! we need a user weight and we set it to the value 1. in this program
! the number of data is all we need
!###############################################################################

program setfwei
implicit none
!------------------------------- variables
integer(kind=4) :: nt
real(kind=4), allocatable,dimension(:),target :: val

!-------------------------------- input
write(*,*) ' enter the number of observations '
read(*,*) nt
open(10,file='fwei',access='direct',recl=4*nt)
allocate(val(nt))
val(:)=1.
write(10,rec=1) val
close(10)
end program setfwei
