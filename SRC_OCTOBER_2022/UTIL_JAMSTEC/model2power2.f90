!=======================================================================
! make an interpolation for getting power 2 sampling ...
!
!
!***********************************************************************
program model
  use s_interpol_mod
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: vel_fmg
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: vel_inp

  integer(kind=4) :: n1inp,n2inp,n3inp
  integer(kind=4) :: n1fmg,n2fmg,n3fmg

  real(kind=4) :: x1inp,x2inp,x3inp
  real(kind=4) :: x1fmg,x2fmg,x3fmg

  real(kind=4) :: h1inp,h2inp,h3inp
  real(kind=4) :: h1fmg,h2fmg,h3fmg

  real(kind=4) :: x1tot,x2tot,x3tot

  real(kind=4) :: vpmin,vpmax

  integer(kind=4) :: i1,i3

write(*,*) ' enter n1inp,n2inp,n3inp '
read(*,*) n1inp,n2inp,n3inp
write(*,*) ' enter x1inp,x2inp,x3inp '
read(*,*) x1inp,x2inp,x3inp
write(*,*) ' enter h1inp,h2inp,h3inp '
read(*,*) h1inp,h2inp,h3inp

x1tot=(n1inp-1)*h1inp
x2tot=(n2inp-1)*h2inp
x3tot=(n3inp-1)*h3inp

!setting up the new model
!same origin
x1fmg=x1inp
x2fmg=x2inp
x3fmg=x3inp
!sampling power of two
call power2(n1inp,n1fmg)
call power2(n3inp,n3fmg)
n1fmg=n1fmg+1   !add one to get the next power
n3fmg=n3fmg+1
n2fmg=n2inp   ! pas de changement selon y
! stepping 
h1fmg=x1tot/float(n1fmg-1)
h2fmg=h2inp
h3fmg=x3tot/float(n3fmg-1)
!======================= rounding
h1fmg=int(h1fmg)! +2.00
h3fmg=int(h3fmg)! +2.00
x1tot=h1fmg*float(n1fmg-1)
x3tot=h3fmg*float(n3fmg-1)
!
!
!
write(*,*) ' size of the previous model ',x1tot,x2tot,x3tot
write(*,*) ' previous sampling ',n1inp,n2inp,n3inp
write(*,*) ' new stepping ',h1inp,h2inp,h3inp
!
!
!
x1tot=(n1fmg-1)*h1fmg
x2tot=(n2fmg-1)*h2fmg
x3tot=(n3fmg-1)*h3fmg
write(*,*) ' size of the new model ',x1tot,x2tot,x3tot
write(*,*) ' new sampling ',n1fmg,n2fmg,n3fmg
write(*,*) ' new stepping ',h1fmg,h2fmg,h3fmg

allocate(vel_inp(n1inp,n2inp,n3inp))
allocate(vel_fmg(n1fmg,n2fmg,n3fmg))

open(8,file='velocity_jamstec',access='direct',recl=4*n1inp*n2inp*n3inp)
read(8,rec=1) vel_inp
close(8)

open(8,file='velocity_jamstec_2D',access='direct',recl=4*n1inp*n3inp)
write(8,rec=1) vel_inp(:,2,:)
close(8)

 call subinterpol(vel_inp,vel_fmg,   &
                  x1inp,x2inp,x3inp,n1inp,n2inp,n3inp,h1inp,h2inp,h3inp,   &
                  x1fmg,x2fmg,x3fmg,n1fmg,n2fmg,n3fmg,h1fmg,h2fmg,h3fmg,vpmin,vpmax)

write(*,*) ' min,max velocities ',vpmin,vpmax

open(8,file='velocity_jamstec_power2',access='direct',recl=4*n1fmg*n2fmg*n3fmg)
write(8,rec=1) vel_fmg
close(8)

open(8,file='velocity_jamstec_power2_2D',access='direct',recl=4*n1fmg*n3fmg)
write(8,rec=1) ((vel_fmg(i1,2,i3),i1=1,n1fmg),i3=1,n3fmg)
!write(8,rec=1) vel_fmg(:,2,:)
close(8)

stop
end program model


SUBROUTINE power2(np,npout)
 INTEGER ::  np2(14)
 DATA (np2(i),i=1,14) /16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072/
 DO i=1,14
 IF (np.le.np2(i)) THEN
 IF (i.eq.1) THEN
 IF (np.lt.np2(i)) THEN
 WRITE(*,*) "input error; nsampl < 8"
 go to 999
 ELSE
 npout=8
 END IF
 ELSE
 npout=np2(i)
 go to 10
 END IF
 END IF
 END DO
10      continue
 if (np.ne.npout) write(*,*) "warning: npout > np"
 RETURN
999     STOP
end subroutine power2




