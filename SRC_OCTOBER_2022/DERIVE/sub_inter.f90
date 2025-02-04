MODULE s_interpol3ds
contains

! #################################################################
!
!
!
!
! #################################################################
subroutine interpol3ds(x1a,x1b,x2a,x2b,x3a,x3b,y1,y2,y3,y4,y5,y6,y7,y8,x1,x2,x3,y)
implicit none
real(kind=4) :: x1a,x1b,x2a,x2b,x3a,x3b,y1,y2,y3,y4,y5,y6,y7,y8
real(kind=4) :: x1,x2,x3,y,t,u,w
!
if(x1b == x1a) then    ! machine precision
  t=0.
else
  t=(x1-x1a)/(x1b-x1a)
endif
if(x2b == x2a) then
  u=0.
else
  u=(x2-x2a)/(x2b-x2a)
endif
if(x3b == x3a) then
  w=0.
else
  w=(x3-x3a)/(x3b-x3a)
endif
y= (1.-t) * (1.-u) * (1.-w)*  y1+   &
      t   * (1.-u) * (1.-w)*  y2+   &
      t   *    u   * (1.-w)*  y3+   &
   (1.-t) *    u   * (1.-w)*  y4+   &
   (1.-t) * (1.-u) *   w   *  y5+   &
      t   * (1.-u) *   w   *  y6+   &
      t   *   u    *   w   *  y7+   &
   (1.-t) *   u    *   w   *  y8  
!
return
end subroutine interpol3ds

END MODULE s_interpol3ds
