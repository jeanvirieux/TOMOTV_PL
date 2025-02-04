MODULE s_aprod
contains
! ###################################################################
! mode =1       dt = dt + A * dma
!
! mode =2       dma = dma + At *dt
!
! ###################################################################
!
!  A est stocke sous forme sparse via rw(nnz),iw(ictot),jw(idtot)
!
!
!               -fmat contient les elements non nuls de la matrice A. 
!               -fic contient les numeros de colonne de chaque elements de fmat:
!                    l'elemement mat(j) est dans la colonne iw(j) de A.
!               -fid est un pointeur servant a extraire une ligne i de la matrice A:
!                    les elements non nuls de la ligne i sont:
!
!                   [mat(jw(i)),mat(jw(i)+1),...,mat(jw(i)+k),...,mat(jw(i+1)-1)]
!
!                Ces trois fichiers suffisent pour reconstituer la matrice:
!                    l'element mat(jw(i)+k) se situe sur la ligne i de A et
!                    se situe sur la colonne iw(jw(i)+k) de A.
!                    
!                              A(i,iw(jw(i)+k))=mat(jw(i)+k).
!
!
!
subroutine aprod(mode,nrow,npar,x,y,leniw,lenjw,lenrw,iw,jw,rw)
!  official     (mode,m   ,n   ,x,y,leniw,lenjw,lenrw,iw,jw,rw)
!                                         =====          ==    added
!
integer(kind=4) :: mode
real(kind=4), dimension(:) ::  rw(:)
integer(kind=4), dimension(:) :: iw(:)     ! alias ic   parameter column
integer(kind=4), dimension(:) :: jw(:)     ! alias id   data row
real(kind=4), dimension(:) :: x(:)
real(kind=4), dimension(:) :: y(:)
integer(kind=4) :: npar,nrow     !alias n et m   (m = data et n = parameters)
integer(kind=4) :: leniw,lenjw,lenrw       ! alias lenic et lenid
integer(kind=4) :: irow,icol

integer(kind=4) :: leniw_0,lenjw_0,lenrw_0,npar_0
leniw_0=leniw;lenjw_0=lenjw;lenrw_0=lenrw;npar_0=npar      ! just for avoiding warnings from compilers

!
if(mode == 1) then
!   write(*,*) 'mode 1'      y = y + A x
  do irow=1,nrow                    ! toujours via les lignes stockage sparse naturel dans ce sens
!    y(irow)=0.
    do itot=jw(irow),jw(irow+1)-1   ! elements non nuls
      icol=iw(itot)                 ! position de la colonne
      y(irow)=y(irow)+ rw(itot)*x(icol)
!JEAN      write(71,*) ' irow,icol, rw ',irow,icol,rw(itot)
    enddo
  enddo
else if (mode == 2) then
!   write(*,*) 'mode 2'      x = x + At y
!  x(:)=0.
  do irow=1,nrow                    ! toujours via les lignes stockage sparse ... x(:) =0 si jamais rien
    do itot=jw(irow),jw(irow+1)-1   ! elements non nuls
      icol=iw(itot)                 ! position de la colonne
      x(icol)=x(icol)+rw(itot)*y(irow)
!JEAN      write(72,*) ' irow,icol, rw ',irow,icol,rw(itot)
    enddo
  enddo
else
   write(*,*) 'pb dans lsqr avec le mode de apro ni un ni deux',mode
   stop
end if

return
end subroutine aprod

end module s_aprod
