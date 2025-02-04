  ! ----------------------------------------------------------------------------------------------------
  ! PROGRAMM SMOOTHGAUSS3DFAST: smooth a 3D grid with a Gaussian filter; smoothing is applied
  ! below topography defined in a binary file of dim(n2,n3)
  ! From Stephane Operto ... 
  ! should not be distributed without his permission
  ! ----------------------------------------------------------------------------------------------------
  
  IMPLICIT NONE

  INTEGER :: n1,n2,n3,itypsmo,i1,i2,i3,il3,il2,il1,k3,k2,k1,j1,j2,j3,jj1,jj2,jj3,it
  REAL :: h,tau1,tau2,tau3,xl3,xl2,xl1,tau12,tau22,tau32,d,betatot
  CHARACTER(LEN=80) :: name_in,name_out,name_topo
  REAL,ALLOCATABLE :: x(:,:,:),x1(:,:,:),topo(:,:),beta3(:),beta2(:),beta1(:)
  INTEGER,ALLOCATABLE :: itopo(:,:)

  ! -----------------------------------------------------------------------------------

  write(*,*) "Input / Output file names (name_in name_out)"
  read(*,*) name_in,name_out

  write(*,*) 'Grid dimensions (n1 n2 n3)'
  read(*,*) n1,n2,n3

  write(*,*) 'Grid interval (h)'
  read(*,*) h

  write(*,*) 'Correlation lengths of Gaussian filter (tau1 tau2)'
  read(*,*) tau1,tau2,tau3

  write(*,*) 'Smooth parameter (0) or its inverse (1)'
  read(*,*) itypsmo

  write(*,*) 'it(0/1) - Use (1) or not (0) topography,Topography file name (name_topo)'
  read(*,*) it,name_topo

  ! -----------------------------------------------------------------------------------

  open(1,file=name_in,access='direct',recl=n1*n2*n3*4)
  open(2,file=name_out,access='direct',recl=n1*n2*n3*4)
  ! ----------------------------------------------------------------
  ! DEBUG
  open(3,file='fdebug',access='direct',recl=n1*n2*4)
  ! ----------------------------------------------------------------
  if (it.eq.1) open(10,file=name_topo,access='direct',recl=n2*n3*4)

  ! -----------------------------------------------------------------------------------

  ALLOCATE (x(n1,n2,n3))
  ALLOCATE (x1(n1,n2,n3))
  ALLOCATE (topo(n2,n3))
  ALLOCATE (itopo(n2,n3))

  write(*,*) 'READ INPUT FD GRID'
  read(1,rec=1) (((x(i1,i2,i3),i1=1,n1),i2=1,n2),i3=1,n3)
  if (it.eq.1) read(10,rec=1) ((topo(i2,i3),i2=1,n2),i3=1,n3)

  x1(:,:,:)=0.

  if (it.eq.1) then
     do i3=1,n3
        do i2=1,n2
           itopo(i2,i3)=nint(topo(i2,i3)/h)+1
           write(*,*) i2,i3,itopo(i2,i3)
        end do
     end do
  else
     itopo(:,:)=1
  end if

  if (tau3.eq.0.and.tau2.eq.0..and.tau1.eq.0.) then

     x1(:,:,:)=x(:,:,:)

  else

     ! velocity to slowness  conversion

     if (itypsmo.eq.1) then
        write(*,*) 'CONVERT TO SLOWNESS'
        do i3=1,n3
           do i2=1,n2
              do i1=1,n1
                 x(i1,i2,i3)=1./x(i1,i2,i3)
              end do
           end do
        end do
     end if

     x1(:,:,:)=x(:,:,:)

     xl3=3.*tau3
     xl2=3.*tau2
     xl1=3.*tau1
     il3=int(xl3/h)+1
     il2=int(xl2/h)+1
     il1=int(xl1/h)+1

     write(*,*) "il3 il2 il1 = ",il3,il2,il1

     ALLOCATE (beta3(2*il3+1))
     ALLOCATE (beta2(2*il2+1))
     ALLOCATE (beta1(2*il1+1))

     tau12=tau1*tau1
     tau22=tau2*tau2
     tau32=tau3*tau3

     k3=0
     do i3=-il3,il3
        k3=k3+1
        d=float(i3)*h
        beta3(k3)=exp(-d**2/tau32)
     end do

     k2=0
     do i2=-il2,il2
        k2=k2+1
        d=float(i2)*h
        beta2(k2)=exp(-d**2/tau22)
     end do

     k1=0
     do i1=-il1,il1
        k1=k1+1
        d=float(i1)*h
        beta1(k1)=exp(-d**2/tau12)
     end do

     write(*,*) '---------------------------------------------------'
     write(*,*) 'SMOOTH ALONG DIMENSION 3'
     write(*,*) '---------------------------------------------------'

     do i3=1,n3
        do i2=1,n2
           do i1=itopo(i2,i3),n1
              x1(i1,i2,i3)=0.
              betatot=0.
              jj3=0
              do j3=-il3,il3
                 jj3=jj3+1
                 k3=j3+i3
                 if (k3.lt.1.or.k3.gt.n3) go to 3
                 x1(i1,i2,i3)=x1(i1,i2,i3)+beta3(jj3)*x(i1,i2,k3)
                 betatot=betatot+beta3(jj3)
3             end do

              x1(i1,i2,i3)=x1(i1,i2,i3)/betatot
           end do
        end do
     end do

     write(*,*) '---------------------------------------------------'
     write(*,*) 'SMOOTH ALONG DIMENSION 2'
     write(*,*) '---------------------------------------------------'


     do i3=1,n3
        do i2=1,n2
           do i1=itopo(i2,i3),n1
              x(i1,i2,i3)=0.
              betatot=0.
              jj2=0
              do j2=-il2,il2
                 jj2=jj2+1
                 k2=j2+i2
                 if (k2.lt.1.or.k2.gt.n2) go to 31
                 x(i1,i2,i3)=x(i1,i2,i3)+beta2(jj2)*x1(i1,k2,i3)
                 betatot=betatot+beta2(jj2)
31            end do

              x(i1,i2,i3)=x(i1,i2,i3)/betatot
           end do
        end do
     end do

     ! -------------------------------------------------------------------------
     ! DEBUG
     !         if (itypsmo.eq.1) then
     !        do i2=1,n2
     !        do i1=itopo(i2),n1
     !        x1(i1,i2)=1./x1(i1,i2)
     !        end do
     !        end do
     !        end if
     !
     !         write(3,rec=1) ((x1(i1,i2),i1=1,n1),i2=1,n2)
     ! ------------------------------------------------------------------------

     write(*,*) '---------------------------------------------------'
     write(*,*) 'SMOOTH ALONG DIMENSION 1'
     write(*,*) '---------------------------------------------------'

     do i3=1,n3
        do i2=1,n2
           do i1=itopo(i2,i3),n1
              x1(i1,i2,i3)=0.
              betatot=0.
              jj1=0
              do j1=-il1,il1
                 jj1=jj1+1
                 k1=j1+i1
                 if (k1.lt.itopo(i2,i3).or.k1.gt.n1) go to 4
                 x1(i1,i2,i3)=x1(i1,i2,i3)+beta1(jj1)*x(k1,i2,i3)
                 betatot=betatot+beta1(jj1)
4             end do

              x1(i1,i2,i3)=x1(i1,i2,i3)/betatot
           end do
        end do
     end do

     if (itypsmo.eq.1) then
        write(*,*) 'CONVERT FROM SLOWNESS TO VELOCITY'
        do i3=1,n3
           do i2=1,n2
              do i1=1,n1
                 x1(i1,i2,i3)=1./x1(i1,i2,i3)
              end do
           end do
        end do
     end if

  end if

  write(*,*) 'WRITE IN OUTPUT FILE'
  write(2,rec=1) (((x1(i1,i2,i3),i1=1,n1),i2=1,n2),i3=1,n3)

  close(2)
  close(1)
  if (it.eq.1) close(10)

  DEALLOCATE (x)
  DEALLOCATE (x1)
  DEALLOCATE (topo)
  DEALLOCATE (itopo)
  DEALLOCATE (beta3)
  DEALLOCATE (beta2)
  DEALLOCATE (beta1)
  write(*,*) 'END'

  stop
end program
