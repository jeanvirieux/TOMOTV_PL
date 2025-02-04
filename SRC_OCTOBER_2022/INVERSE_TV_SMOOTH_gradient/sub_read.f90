MODULE s_read
contains

!############################## READ INPUT #############################
!     
!     input : flog : unite logique fichier log
!                i : compteur unites logique
!
!     output : xo,yo,zo : origine du modele
!              n1,n2,n3 : grille du modele
!              h1,h2,h3 : pas du modele

!              ixo,iyo,izo  : origine volume d'inversion
!              ni1,ni2,ni3  : noeuds du volume d'inversion 
!              nvit : nbre total de noeuds a inverser
!              uvp  : unite logique de modele P ou 1./U
!              fvp  : name of the P velocity file
!              uvs  : unite logique de modele S ou 1./V  
!              fvs  : name of the S velocity file
!              damp : damping lsqr
!
!
subroutine read_par(iunit,xo,yo,zo,n1,n2,n3,h1,h2,h3,ixo,iyo,izo,ni1,ni2,ni3, &
                    nvit,uvp,uvs,fvp,fvs,damp, &
                    dvp_max,dvs_max,dhori_max,dvert_max,dto_max,itnlim,flog)
implicit none
integer(kind=4) :: iunit,ni1,ni2,ni3,nvit,ixo,iyo,izo,n1,n2,n3,uvp,uvs,flog
real(kind=4) :: xo,yo,zo,h1,h2,h3,damp
REAL(kind=4) :: dvp_max,dvs_max,dhori_max,dvert_max,dto_max 
integer(kind=4) :: itnlim                                   
integer(kind=4) :: chx,iloop
character(len=132) :: fvp,fvs
character(len=1) :: carac

open(49,file='inversion.par',status='old',err=1000)
read(49,*) chx
read(49,*) xo,yo,zo
read(49,*) n1,n2,n3
read(49,*) h1,h2,h3
!--------------------- read P velocity file name
read(49,'(a)') fvp
uvp=iunit
iunit=iunit+1

read(49,*) ixo,iyo,izo
read(49,*) ni1,ni2,ni3

if(chx.eq.2) then
  uvs=iunit 
  iunit=iunit+1
  read(49,'(a)') fvs
else
  uvs=0
endif

open(49,file='inversion.head',status='old',err=2000)
do iloop=1,21
read(49,'(a)') carac
enddo
read(49,*) damp
write(flog,*) ' damping for lsqr ',damp
! ##############################################################
read(49,'(a)') carac
read(49,*) itnlim
write(flog,*) ' maximum number of iterations for lsqr ',itnlim
read(49,'(a)') carac
read(49,*)  dvp_max
write(flog,*) ' Maximum perturbation for P wave velocity ',dvp_max
read(49,'(a)') carac
read(49,*)  dvs_max
write(flog,*) ' Maximum perturbation for S wave velocity ',dvs_max
read(49,'(a)') carac
read(49,*)  dhori_max
write(flog,*) ' Maximum perturbation for horizontal shift ', dhori_max
read(49,'(a)') carac
read(49,*)  dvert_max
write(flog,*) ' Maximum perturbation for vertical shift ',dvert_max
read(49,'(a)') carac
read(49,*)  dto_max
write(flog,*) ' Maximum perturbation for origin time ',dto_max
close(49)
! ###############################################################
!
if(ni1.gt.0.and.ni2.gt.0.and.ni3.gt.0) then
  nvit=ni1*ni2*ni3
else
  nvit=0
endif
return
!================================ error for the file inversion.par
1000 continue
write(flog,*) ' error in reading the file inversion.par '
write(*,*) ' error in reading the file inversion.par '
stop
2000 continue
write(flog,*) ' error in reading the file inversion.head '
write(*,*) ' error in reading the file inversion.head '
stop
end subroutine read_par

! ###############################################################
!--
!     read residuals
!
!     input : nt : nombre de lignes
!             iu : unite logique du fichier fdift
!
!     output : res(nl) : residus
!
subroutine read_fdift(iu,res,nt)
implicit none

integer(kind=4) :: iu,nt
real(kind=4), dimension(:) :: res(:)

open(iu,file='fdift',access='direct',recl=4*nt,status='old',err=1000)
read(iu,rec=1) res
close(iu)
return
1000 continue
write(*,*) ' error missing residues file fdift '
stop

end subroutine read_fdift

! ###############################################################
!--
!     read normalisation matrix
!
!     input : n  : nombre de colonnes
!             iu : unite logique du fichier fnorm
!
!     output : renorm (n alias nc) : vecteur de normalisation.
!
subroutine read_renorm(iu,renorm,n)
implicit none

integer(kind=4) :: iu,n
real(kind=4) :: renorm(n)

open(iu,file='fnorm',access='direct',recl=4*n)
read(iu,rec=1) renorm
close(iu)

return
end subroutine read_renorm

!#####################################  UPDATE NEW MODEL ####################
!
! renormalisation du resultat de l'inversion
! x(n) resultat de l'inversion
! y(n) coefficients de renormalisation
!
!  on devait resoudre Ax=b,
!  on preconditionnait ACy=b, (C diagonale)
!  LSQR calcule y, et pour trouver x,
!  on renormalise y, x=Cy.
!  ( C(i,i)=1/y(i) ,C(i,j)=0 si i<>j)
!
subroutine subrenorm(x,y,n)
implicit none

integer(kind=4) :: j,n
real(kind=4), dimension(:) :: x(n),y(n)
do j=1,n
!JEAN DEBUG  write(78,*) ' j,sol, scale ',j,x(j),y(j)
  x(j)=x(j)*y(j)
enddo

return
end subroutine subrenorm

END MODULE s_read
