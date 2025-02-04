!==================================================
!  transformation of velocities (hs(n)) into slowness
!  used by the P&L eikonal solver (hs(n))
!
!  n1,n2,n3 : model size
!  hs(n) : velocity matrix (being modified)
!  h : P&L grid step (cubic geometry)
!==================================================
MODULE s_hs
contains

subroutine sub_hs(hs,n1,n2,n3,h)
implicit none
integer(kind=4) :: n1,n2,n3
integer(kind=4) :: i1,i2,i3
real(kind=4) :: h
real(kind=4) :: hs(n1,n2,n3)

Do i1=1,n1
Do i2=1,n2
Do i3=1,n3
if(hs(i1,i2,i3)== 0.)then
hs=0.
else
hs(i1,i2,i3)=h/hs(i1,i2,i3)
endif
Enddo
Enddo
Enddo
return
end subroutine sub_hs

END MODULE s_hs
