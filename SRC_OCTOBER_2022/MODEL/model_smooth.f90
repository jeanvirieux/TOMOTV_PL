!=======================================================================
!
!***********************************************************************
!
!***********************************************************************
program model_smooth
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modsmooth
  !======================================================= smoothing
  real(kind=4), dimension(:), allocatable :: input,output
  
  !----------------------------- variables
  real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  integer(kind=4) n1inv,n2inv,n3inv,chx   ! dimensions for inversion grid 
  integer(kind=4) flog,finput             ! flux entree et sortie

  integer(kind=4)                    :: i1,i2,i3
  character(len=1)                   :: carac

  integer(kind=4) :: i,j,k,ir
  integer(kind=4) :: nrepeat,ir1,ir2,ir3

  logical :: file_exists

  flog=51
  open(flog,file='flog.model_smooth')
  finput=flog+1
  ! #################################### choix des ondes
  open(finput,file='model.head',status='old')
  read(finput,'(a)') carac
  read(finput,*) chx
  write(flog,*) ' flag 1 pour ondes P, 2 pour ondes P et S',chx
  ! ####################################  input
  read(finput,'(a)') carac
  read(finput,*) x1inv,x2inv,x3inv
  write(flog,*) ' (x1,x2,x3) origin of the inversion grid ',x1inv,x2inv,x3inv
  read(finput,'(a)') carac
  read(finput,*) n1inv,n2inv,n3inv
  write(flog,*) ' (n1,n2,n3) inversion grid ',n1inv,n2inv,n3inv
  read(finput,'(a)') carac
  read(finput,*) h1inv,h2inv,h3inv
  write(flog,*) ' (h1,h2,h3) grid sampling ', h1inv,h2inv,h3inv
  close(finput)


  ALLOCATE(modsmooth(n1inv,n2inv,n3inv))
  ! #####################################  array P
  ALLOCATE(modinvp(n1inv,n2inv,n3inv))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then                  !###  array S
     ALLOCATE(modinvs(n1inv,n2inv,n3inv))
  endif

  ! #####################################  io
  open(7,file='modelP',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  !------------ if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     open(9,file='modelS',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
  endif

  ! #####################################  reading the inverse model
  read(7,rec=1) modinvp             ! true grid in P
  modinvp(:,:,:)=1./modinvp(:,:,:)
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     read(9,rec=1) modinvs             ! true grid in S
     modinvs(:,:,:)=1./modinvs(:,:,:) 
  endif
  close(7)
  close(9)

  INQUIRE(FILE="model_smooth.head", EXIST=file_exists)

  if(file_exists) then
     open(71,file='model_smooth.head',status='old')
     read(71,*) nrepeat,ir1,ir2,ir3
     close(71)
  else
     nrepeat=2; ir1=2; ir2=2; ir3=2
  endif

  do ir=1,nrepeat
     ! smoothing along x
     allocate(input(n1inv)); allocate(output(n1inv))
     do k=1,n3inv               ! loop over invertable nodes
        do j=1,n2inv

           do i=1,n1inv
              input(i)=modinvp(i,j,k)
           enddo
           call triangle(ir1,n1inv,input,output)
           do i=1,n1inv
              modinvp(i,j,k)=output(i)
           enddo
        enddo   ! j=1,n2inv
     enddo   ! k=1,n3inv
     deallocate(input); deallocate(output)

     ! smoothing along y
     allocate(input(n2inv)); allocate(output(n2inv))
     do k=1,n3inv               ! loop over invertable nodes
        do i=1,n1inv
           do j=1,n2inv
              input(j)=modinvp(i,j,k) 
           enddo
           call triangle(ir2,n2inv,input,output)
           do j=1,n2inv
              modinvp(i,j,k)=output(j)
           enddo
        enddo   ! i=1,n1inv
     enddo   ! k=1,n3inv
     deallocate(input); deallocate(output)

     ! smoothing along z
     allocate(input(n3inv)); allocate(output(n3inv))
     do j=1,n2inv               ! loop over invertable nodes
        do i=1,n1inv

           do k=1,n3inv
              input(k)=modinvp(i,j,k)  
           enddo
           call triangle(ir3,n3inv,input,output)
           do k=1,n3inv
              modinvp(i,j,k)=output(k)
           enddo

        enddo   ! i=1,n1inv
     enddo   ! j=1,n2inv
     deallocate(input); deallocate(output)
  enddo ! ir
  
  modinvp(:,:,:)=1./modinvp(:,:,:)
  write(*,*) ' min,max P values',MINVAL(modinvp),MAXVAL(modinvp)
  open(7,file='modelP.smooth',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
  write(7,rec=1) modinvp
  close(7)

  !=============================================
  ! end of smoothing model P
  !=============================================

  if(chx.eq.2) then
     do ir=1,nrepeat
        ! smoothing along x
        allocate(input(n1inv)); allocate(output(n1inv))
        do k=1,n3inv               ! loop over invertable nodes
           do j=1,n2inv

              do i=1,n1inv
                 input(i)=modinvs(i,j,k)
              enddo
              call triangle(ir1,n1inv,input,output)
              do i=1,n1inv
                 modinvs(i,j,k)=output(i)
              enddo
           enddo   ! j=1,n2inv
        enddo   ! k=1,n3inv
        deallocate(input); deallocate(output)

        ! smoothing along y
        allocate(input(n2inv)); allocate(output(n2inv))
        do k=1,n3inv               ! loop over invertable nodes
           do i=1,n1inv
              do j=1,n2inv
                 input(j)=modinvs(i,j,k) 
              enddo
              call triangle(ir2,n2inv,input,output)
              do j=1,n2inv
                 modinvs(i,j,k)=output(j)
              enddo
           enddo   ! i=1,n1inv
        enddo   ! k=1,n3inv
        deallocate(input); deallocate(output)

        ! smoothing along z
        allocate(input(n3inv)); allocate(output(n3inv))
        do j=1,n2inv               ! loop over invertable nodes
           do i=1,n1inv

              do k=1,n3inv
                 input(k)=modinvs(i,j,k)  
              enddo
              call triangle(ir3,n3inv,input,output)
              do k=1,n3inv
                 modinvs(i,j,k)=output(k)
              enddo

           enddo   ! i=1,n1inv
        enddo   ! j=1,n2inv
        deallocate(input); deallocate(output)
     enddo ! ir

     !=============================================
     ! end of smoothing model S
     !=============================================
     modinvs(:,:,:)=1./modinvs(:,:,:)
     write(*,*) ' min,max S values',MINVAL(modinvs),MAXVAL(modinvs)
     open(7,file='modelS.smooth',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(7,rec=1) modinvs
     close(7)

  endif

contains

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
       write(*,*) "boxconv: error in the length of input!"
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

end program model_smooth






