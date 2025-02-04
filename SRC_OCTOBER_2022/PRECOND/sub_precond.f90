MODULE s_precond
contains

!===================================================================
subroutine sub_D1(d1,r,poids,nl,r0,r1,long,l0,l1,coeff)
implicit none
integer(kind=4) :: j,nl
real(kind=4) :: r0,r1,l0,l1,coeff
real(kind=4), dimension(:) :: d1(:),poids(:),long(:),r(:)
!
! Ponderation en fonction des residus
!
do j=1,nl
  !write(76,*) j,poids(j),r(j),long(j) !STEPHANIE
  if(abs(r(j)) <= r0) then
    d1(j)=1.*poids(j)
  elseif((abs(r(j)) >  r0) .and. (abs(r(j)) < r1))then
    d1(j)=((r1-abs(r(j))) / (r1-r0)) * poids(j)
  elseif(abs(r(j)) >= r1) then
    d1(j)=0.
  endif
enddo
!
! Ponderation en fonction de la longueur du rai
!
do j=1,nl
  if((long(j) < l1) .and. (long(j) > l0))then
    d1(j)=d1(j)*((l1-long(j)) / (l1-l0))
  elseif(long(j)  >= l1) then
    d1(j)=0.
  endif
enddo
! Ajout pour ponderation avec les DD ... coeff=1. dans le cas standard  (en attente Jean)
!do j=1,nl  !! stephanie
!  d1(j)=d1(j)*coeff
!enddo
coeff=1.
return
end subroutine sub_D1
! =============================================================================
!  Building the matrix D2 for getting adimensional gradients for each parameter class
!   d is the D2 matrix and has for dimension d(nc)
!          nc column
!          nl row   (ligne)
!   vip,vis,sei,to are integer offsets for the different classes of parameters
!   mat,ic,id,nnz,nl are standard inputs for the sparse matrix A
!   cp,cs,cpo,cto are a priori scaling required by the user through experimental
!                     numerical testing
! =================================================================
subroutine sub_D2(d,nc,vip,vis,sei,to,mat,ic,id,nnz,nl,cp,cs,cpo,cto)
  implicit none
  integer(kind=4) :: nc,sei,to,nnz,nl,vip,vis
  integer(kind=4), dimension(:) :: ic(:),id(:)
  real(kind=4), dimension(:) :: d(:),mat(:)
  real(kind=4) :: cp,cs,cpo,cto
  !
  integer(kind=4) :: i,id0,nl0
  real(kind=4) :: nrp,nrs,nre,nrt

  id0=id(1)                  ! just for avoiding warning from compilers
  nl0=nl                     ! idem 
  !
  ! We must take care of the different groups of parameters
  !
  nrp=-1.e+29
  nrs=-1.e+29
  nre=-1.e+29
  nrt=-1.e+29
  if(vip < 0) vip=0
  !write(*,*) ' JEAN ',nnz,vip,vis,sei,to
  do i=1,nnz
     !JEAN
     !write(77,*) ic(i),mat(i)
     d(ic(i))=d(ic(i))+mat(i)*mat(i)  ! gradient**2 is computed for each parameter
     if(ic(i) <= vip) then            ! Up or Vp      SECOND/METER OR  METER/SECOND
        nrp=max(nrp,d(ic(i)))
     elseif(ic(i) > vip .and. ic(i) <= vip+vis) then ! Us or Vs  SECOND/METER OR METER/SECOND
        nrs=max(nrs,d(ic(i)))
     elseif(ic(i) > vip+vis .and. ic(i) <= vip+vis+sei) then   ! xo,yo,zo      METER
        nre=max(nre,d(ic(i)))
     elseif(ic(i) >  vip+vis+sei .and. ic(i) <= vip+vis+sei+to) then   ! to    SECOND
        nrt=max(nrt,d(ic(i)))
     else 
        write(*,'(a)') ' WARNING1: pb sub_precond D2'
        write(*,'(a11,4i10)') 'i,ic(i),vip,nnz',i,ic(i),vip,nnz
        write(*,*) 'WARNING1: pb sub_precond D2'
     endif
  end do
  !JEAN
  !close(77)
  write(*,*) ' maximum value for gradient square per category when setup '
  if(nrp > -1.e+29) write(*,*) '  Vp sum of gradient**2 and dimension of subspace ',nrp,vip
  if(nrs > -1.e+29) write(*,*) '  Vs sum of gradient**2 and dimension of subspace ',nrs,vis
  if(nre > -1.e+29) write(*,*) ' pos sum of gradient**2 and dimension of subspace ',nre,sei
  if(nrt > -1.e+29) write(*,*) '  to sum of gradient**2 and dimension of subspace ',nrt,to

  !
  !  scaling by user values and normalizing
  !
  do i=1,nc
     if(i <= vip) then           ! Up ou Vp     on normalise les gradients par classe
        d(i)=cp/sqrt(nrp)         !              ET PAS PAR COMPOSANTE
     elseif(i > vip .and. i <= vip+vis) then   ! Us ou Vp
        d(i)=cs/sqrt(nrs)
     elseif(i > vip+vis .and. i <= vip+vis+sei) then  ! xo,yo,zo
        d(i)=cpo/sqrt(nre)
     elseif(i > vip+vis+sei .and. i <= vip+vis+sei+to) then ! to
        d(i)=cto/sqrt(nrt)
     else 
        write(*,'(a)') ' WARNING2: pb sub_precond D2 ... '
        write(*,'(a)') ' unexpected index for normalization over parameters'
        write(*,*) 'WARNING2: pb sub_precond D2',i,vip,vip+vis+sei+to,nc
     endif
  enddo
  return
end subroutine sub_D2

!-------------------------------------------------------------------------------
!     D*M
!     multiplication a gauche d'une matrice sparse par une matrice diagonale
!-------------------------------------------------------------------------------
subroutine mut_d_m(d,mat,id,nnz,nl)
implicit none
integer(kind=4) :: i,k,nnz,nl,nnz0
integer(kind=4), dimension(:) :: id(:)      ! id(nl+1)
real(kind=4), dimension(:) :: d(:),mat(:)   ! d(nl),mat(nnz)
nnz0=nnz                     ! just to avoid warning from compilers 
do i=1,nl
  do k=id(i),id(i+1)-1
    mat(k)=d(i)*mat(k)
  enddo
enddo
return
end subroutine mut_d_m

!-------------------------------------------------------------------------------
!     M*D
!     multiplication a droite d'une matrice sparse par une matrice diagonale
!-------------------------------------------------------------------------------
subroutine mut_m_d(d,mat,ic,nnz,nc)
implicit none
integer i,nnz,nc,nc0
integer(kind=4), dimension(:) :: ic(:)       ! ic(nnz)
real(kind=4), dimension(:) :: d(:),mat(:)    !d(nc),mat(nnz)
nc0=nc                      ! just to avoid warning from compilers
do i=1,nnz
  mat(i)=mat(i)*d(ic(i))
end do
return
end subroutine mut_m_d

!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine mut_d_v(d,v,n)
implicit none 
integer(kind=4) :: i,n
real(kind=4), dimension(:) :: d(:),v(:)
do i=1,n
  v(i)=v(i)*d(i)
enddo
return
end subroutine mut_d_v


!===============================================================================
! combining poids and Rpoids   coming from different analysis (data and statistics
!
!===============================================================================
subroutine sub_m_p(poids,Rpoids,n)
implicit none
integer(kind=4) :: i,n
real(kind=4), dimension(:) :: poids(:),Rpoids(:)
do i=1,n
  poids(i)=poids(i)*Rpoids(i)
enddo
return
end subroutine sub_m_p

END MODULE s_precond
