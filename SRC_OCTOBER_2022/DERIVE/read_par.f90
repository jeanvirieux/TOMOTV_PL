! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIOUS SUBROUTINES FOR INPUT AND OUTPUT
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE s_read_par
contains

!############################################################################
!
!     READ FILE : INVERSION.PAR
!
!
!     input :
!     flog : unite logique du fichier log
!     i: compteur unite logique de fichier
! 
!
!     output :
!     xo,yo,zo    : origine du modele
!     n1,n2,n3    : grille du modele de vitesse
!     h1,h2,h3    : pas 
!     ixo,iyo,izo : origine de la grille d'inversion
!     ni1,ni2,ni3 : grille d'inversion
!     fvp, fvs    : nom des fichiers pour modele P et modele S
!     flog : unite logique du fichier log
!
! #############################################################################
subroutine read_par(iu, chx, fvp, fvs, &
           xo, yo, zo, n1i, n2i, n3i, h1i, h2i, h3i, &
           ixo,iyo,izo,ni1,ni2,ni3,flog)
implicit none
integer(kind=4) :: iu,n1,n2,n3,flog,chx
integer(kind=4) ::  n1i,n2i,n3i,ixo,iyo,izo,ni1,ni2,ni3
real(kind=4) :: xo,yo,zo,h1,h2,h3
real(kind=4) :: xoi,yoi,zoi,h1i,h2i,h3i
integer(kind=4) :: rmax_dum
character(len=132) :: fvp,fvs
!real(kind=4), allocatable, dimension(:,:,:) :: vit_P,vit_P_inv

!
! lecture des parametres             INTERET DE CETTE LECTURE ?
!
write(flog,'(/35x,a)') 'READ PARAMETERS'
open(iu,file='timefd.par',status='old')
read(iu,*) chx          ! P=1 or P&S=2
read(iu,*) rmax_dum     ! now we know the number of points along the longest ray
read(iu,*) xo,yo,zo     ! origin coordinates
read(iu,*) n1,n2,n3     ! number of nodes in each direction 
read(iu,*) h1,h2,h3     ! space stepping in each direction (could be equal)
read(iu,'(a)') fvp      ! name of the velocity file  (for P waves) 
if(chx == 2) then
read(iu,'(a)') fvs      ! name of the velocity file  (for S waves)
endif
close(iu)

!
! lecture des options d'inversion
!
open(iu,file='inversion.par',status='old')
read(iu,*) chx          ! P=1 or P&S=2   should be consistent !
read(iu,*) xoi,yoi,zoi  !
read(iu,*) n1i,n2i,n3i
read(iu,*) h1i,h2i,h3i
read(iu,'(a)') fvp
read(iu,*) ixo,iyo,izo
read(iu,*) ni1,ni2,ni3
if(chx == 2) then
read(iu,'(a)') fvs
endif
close(iu)
!
! on should check for coherence between timefd and inversion files
!
write(flog,'(/20x,a)') 'read timefd.par'
write(flog,'(10x,a9,3f10.2)') 'origin  :',xo,yo,zo
write(flog,'( 10x,a9,3i10)')  'nodes   :',n1,n2,n3
write(flog,'( 10x,a9,f10.2)') 'step    :',h1,h2,h3

return
end subroutine read_par

END MODULE s_read_par

