!***********************************************************************
! subinterpol : interpolation of the forward square velocity grid from the
!     model velocity structure which is a rectangular velocity grid
!     using a bilinear interpolation ..... of course we must copy last
!     right points of boundaries
!
! ..... v_modele                            parameter model
! ..... v_forward                           forward model
! ..... x1inv,x2inv,x3inv                   origin for the parameter model
! ..... h1inv,h2inv,h3inv,n1inv,n2inv,n3inv values for the rectangular grid
!...... x1orig,x2orig,x3orig                origin for the forward model
! ..... n1,n2,n3,h1,h2,h3                   h1=h2=h3 for a square grid
!***********************************************************************

MODULE s_interpol_mod
contains

subroutine subinterpol(v_modele,v_forward,x1inv,x2inv,x3inv, &
                       n1inv,n2inv,n3inv,h1inv,h2inv,h3inv, &
                       x1orig,x2orig,x3orig,n1,n2,n3,h1,h2,h3,vmin,vmax)
implicit none
integer(kind=4) :: n1inv,n2inv,n3inv,n1,n2,n3
real(kind=4) :: v_modele(n1inv,n2inv,n3inv)
real(kind=4) :: v_forward(n1,n2,n3)
real(kind=4) :: h1inv,h2inv,h3inv,h1,h2,h3
real(kind=4) :: x1inv,x2inv,x3inv,x1orig,x2orig,x3orig
!
integer(kind=4) :: i,j,k
real(kind=4) :: x1_f,x2_f,x3_f,vit
real(kind=4) :: vmin,vmax

vmin=1.e+29; vmax=-1.e29
! ------------------------ right interpolation
v_forward(:,:,:)=0.
do k=1,n3-1
  x3_f=(k-1)*h3+x3orig
  do j=1,n2-1
    x2_f=(j-1)*h2+x2orig
    do i=1,n1-1
      x1_f=(i-1)*h1+x1orig
      call vitest(v_modele,x1inv,x2inv,x3inv,    &
           n1inv,n2inv,n3inv,h1inv,h2inv,h3inv, &
           x1_f,x2_f,x3_f,vit) ! get value at this point from the original grid
      v_forward(i,j,k)=vit              ! put it in the new grid
      if(vit > vmax) vmax=vit
      if(vit < vmin) vmin=vit
    enddo
! ------------------- plane j,k
    v_forward(n1,j,k)=v_forward(n1-1,j,k)
  enddo
! ------------------- plane i,k
  do i=1,n1
    v_forward(i,n2,k)=v_forward(i,n2-1,k)
  enddo
enddo
! ------------------- plane i,j
do i=1,n1
  do j=1,n2
    v_forward(i,j,n3)=v_forward(i,j,n3-1)
  enddo
enddo

return
end subroutine subinterpol
! ***********************************************************************
! estimation of the velocity in the modele structure which
!     is the rectangular grid in this simple case
! ***********************************************************************
subroutine vitest(v_modele,x1inv,x2inv,x3inv,n1inv,n2inv,n3inv, &
                  h1inv,h2inv,h3inv,x1_f,x2_f,x3_f,vit)
implicit none
integer(kind=4) :: n1inv,n2inv,n3inv
real(kind=4) :: v_modele(n1inv,n2inv,n3inv)
real(kind=4) :: h1inv,h2inv,h3inv,x1inv,x2inv,x3inv
real(kind=4) :: x1_f,x2_f,x3_f,vit
!
integer(kind=4) :: icur,jcur,kcur
real(kind=4) x1,x2,x3
! ----------------------------- locate the x1_f,x2_f,x3_f point in the grid
icur=int((x1_f-x1inv)/h1inv+1.e-6)+1
jcur=int((x2_f-x2inv)/h2inv+1.e-6)+1
kcur=int((x3_f-x3inv)/h3inv+1.e-6)+1 ! should be equal to one for 2D geometry
! ----------------------------- last segments of grid are not allowed
if(icur >= n1inv) then
  write(*,*) ' error in vitest icur,n1inv ',icur,n1inv,x1_f
  stop
endif 
if(jcur >= n2inv) then
  write(*,*) ' error in vitest jcur,n2inv ',jcur,n2inv,x2_f
  stop
endif 
if(kcur >= n3inv) then
  write(*,*) ' error in vitest kcur,n3inv ',kcur,n3inv,x3_f
  stop
endif 

! ----------------------------- last segments of grid are not allowed
if(icur <= 0) then
  write(*,*) ' error in vitest icur,n1inv ',icur,n1inv,x1_f
  stop
endif 
if(jcur <= 0) then
  write(*,*) ' error in vitest jcur,n2inv ',jcur,n2inv,x2_f
  stop
endif 
if(kcur <= 0) then
  write(*,*) ' error in vitest kcur,n3inv ',kcur,n3inv,x3_f
  stop
endif

! ----------------------------- linear interpolation
! ----------------------------- reduced coordinates x1,x2,x3
x1=((x1_f-x1inv)-(icur-1)*h1inv)/h1inv
x2=((x2_f-x2inv)-(jcur-1)*h2inv)/h2inv
x3=((x3_f-x3inv)-(kcur-1)*h3inv)/h3inv
! ----------------------------- interpolation in slowness
if(n3inv /= 1) then
vit=1./v_modele(icur  ,jcur  ,kcur  )*(1.-x1)*(1.-x2)*(1.-x3) &
   +1./v_modele(icur+1,jcur  ,kcur  )*    x1 *(1.-x2)*(1.-x3) &
   +1./v_modele(icur  ,jcur+1,kcur  )*(1.-x1)*    x2 *(1.-x3) &
   +1./v_modele(icur+1,jcur+1,kcur  )*    x1 *    x2 *(1.-x3) &
   +1./v_modele(icur  ,jcur  ,kcur+1)*(1.-x1)*(1.-x2)*    x3  &
   +1./v_modele(icur+1,jcur  ,kcur+1)*    x1 *(1.-x2)*    x3  &
   +1./v_modele(icur  ,jcur+1,kcur+1)*(1.-x1)*    x2 *    x3  &
   +1./v_modele(icur+1,jcur+1,kcur+1)*    x1 *    x2 *    x3  
vit=1./vit
else
vit=1./v_modele(icur  ,jcur  ,kcur  )*(1.-x1)*(1.-x2)*(1.-x3) & !x3=0 in fact
   +1./v_modele(icur+1,jcur  ,kcur  )*    x1 *(1.-x2)*(1.-x3) &
   +1./v_modele(icur  ,jcur+1,kcur  )*(1.-x1)*    x2 *(1.-x3) &
   +1./v_modele(icur+1,jcur+1,kcur  )*    x1 *    x2 *(1.-x3)
vit=1./vit
endif
if(vit < 1.e-6 .or. vit > 1.e+6) then
  write(*,*) ' funny value for velocity ',vit,icur,jcur,kcur
  write(*,*) v_modele(icur  ,jcur  ,kcur)
  write(*,*) v_modele(icur+1,jcur  ,kcur)
  write(*,*) v_modele(icur+1,jcur+1,kcur)
  write(*,*) v_modele(icur  ,jcur+1,kcur)
  write(*,*) v_modele(icur  ,jcur  ,kcur+1)
  write(*,*) v_modele(icur+1,jcur  ,kcur+1)
  write(*,*) v_modele(icur+1,jcur+1,kcur+1)
  write(*,*) v_modele(icur  ,jcur+1,kcur+1)
  write(*,*) ' please check your indexes again '
  stop
endif
end subroutine vitest
END MODULE s_interpol_mod
