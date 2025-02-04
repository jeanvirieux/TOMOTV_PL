!=============================================================
!
!=============================================================
program homo

  implicit none

  integer(kind=4), parameter :: ik=4

  real(kind=ik) :: xorg,yorg,zorg
  integer(kind=ik) :: n1,n2,n3,iopt_P_S
  real(kind=ik) :: speed
  character(len=2) :: carac
  real(kind=ik),allocatable,dimension(:) :: vel



  !=============================================== we need to be more robust against input errors
  open(9,file='model.head',status='old')
  read(9,'(a)') carac
  read(9,*) iopt_P_S
  read(9,'(a)') carac
  read(9,*) xorg,yorg,zorg
  read(9,'(a)') carac
  read(9,*) n1,n2,n3
  close(9)
  !===============================================

  write(*,*) ' enter the velocity value '
  read(*,*) speed

  allocate(vel(n1*n2*n3))
  vel(:)=speed

  open(8,file='modelP.ini',access='direct',recl=n1*n2*n3*ik)
  write(8,rec=1) vel
  close(8)
  deallocate(vel)

end program homo

