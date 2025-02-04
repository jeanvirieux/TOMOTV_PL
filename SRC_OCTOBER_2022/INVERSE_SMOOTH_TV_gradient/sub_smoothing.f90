!=================================================================
!apply triangle filter to input xx: yy=triangle*xx
! from Pengliang ... and from Fomel
!=================================================================
subroutine triangle(nbox,nd,xx,yy)
  implicit none

  integer, intent(in)::nbox,nd
  integer::i,np,nq
  real,dimension(nd),intent(in)::xx
  real,dimension(nd),intent(out)::yy
  real,dimension(:),allocatable::pp,qq

  allocate(pp(nd+nbox-1))
  allocate(qq(nd+2*nbox-2))
  call boxconv(nbox,nd,xx,pp)
  np=nbox+nd-1
  call boxconv(nbox,np,pp,qq)
  nq=nbox+np-1
  do i=1,nd
     yy(i)=qq(i+nbox-1)
  enddo
  do i=1,nbox-1
     yy(i)=yy(i)+qq(nbox-i) !fold back near end
  enddo
  do i=1,nbox-1
     yy(nd-i+1)=yy(nd-i+1)+qq(nd+(nbox-1)+i) !fold back far end
  enddo
  deallocate(pp)
  deallocate(qq)
end subroutine triangle
!==================================================================
subroutine boxconv(nbox, nx, xx, yy)
  implicit none
  integer,intent(in):: nx,nbox
  integer::i,ny
  real,dimension(nx),intent(in)::xx
  real,dimension(nx+nbox-1),intent(out)::yy
  real,dimension(:),allocatable::bb

  allocate(bb(nx+nbox))
  if (nbox < 1 .or. nbox > nx) then
     write(*,*) "boxconv: error in the length of input!",nbox,nx
  endif
  ny = nx+nbox-1
  do i= 1, ny
     bb(i) = 0.
  enddo
  bb(1) = xx(1)
  do i= 2, nx
     bb(i) = bb(i-1) + xx(i)  ! make B(Z) = X(Z)/(1-Z)
  enddo
  do i= nx+1, ny
     bb(i) = bb(i-1)
  enddo
  do i= 1, nbox
     yy(i) = bb(i)
  enddo
  do i= nbox+1, ny
     yy(i) = bb(i) - bb(i-nbox) ! make Y(Z) = B(Z)*(1-Z**nbox)
  enddo
  do i= 1, ny
     yy(i) = yy(i) / nbox
  enddo
  deallocate(bb)
end subroutine boxconv

