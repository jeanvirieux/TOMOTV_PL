MODULE s_write_rms
contains
!
!###############################################################
!
subroutine write_rms(r,p,n)
implicit none
integer(kind=4) :: i,n
real(kind=4) :: som,rms
real(kind=4), dimension(:) :: r(:),p(:) ! r(n),p(n)   
rms=0.
som=0.
open(49,file='rms.pond',access='append')
do i=1,n
  rms=rms+r(i)**2
  som=som+p(i)
!write(75,*) r(i),p(i),som,rms !STEPHANIE
end do
!close(75) !STEPHANIE
write(*,*) 'END computation rms: rms,som:',rms,som
rms=sqrt(rms/som)
write(49,*) rms   
close(49)
return
end subroutine write_rms
END MODULE s_write_rms
