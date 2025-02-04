program smooth2d
  implicit none

  integer::n1,n2,r1,r2, repeat
  real,dimension(:,:),allocatable::mod0,mod1
  integer::i1,i2,i
  character(len=80)::inputfile,outputfile

  read(*,*) inputfile,outputfile
  read(*,*) n1,n2,r1,r2,repeat !size of data, smoothing radius,repeating times
!  print *,inputfile,outputfile,n1,n2,r1,r2,repeat


  allocate(mod0(n1,n2))
  allocate(mod1(n1,n2))

  open(10,file=inputfile,access='direct',recl=4*n1*n2)
  read(10,rec=1) mod0
  close(10)
  do i=1,repeat
     do i2=1,n2
        call triangle(r1,n1,mod0(:,i2),mod1(:,i2))
     enddo
     do i1=1,n1
        call triangle(r2,n2,mod1(i1,:),mod0(i1,:))
     enddo
  enddo

  open(10,file=outputfile,access='direct',recl=4*n1*n2,status='replace')
  write(10,rec=1) mod0
  close(10)


  deallocate(mod0)
  deallocate(mod1)

  call exit(0)
end program smooth2d

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
     write(0,*) "boxconv: error in the length of input!"
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

!=================================================================
!apply triangle filter to input xx: yy=triangle*xx
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
