MODULE s_lap
contains
!###############################################################################
!
!           construction de la matrice de lissage par laplacien
!
!---parametres:-----------------------------------------------------------------
!
!   nnz1 : nombre de parametres de la matrice de lissage (variable)             
!   nl1  : nombre de lignes de la matrice de lissage     (variable)
!   lap  : type de lissage (1,2 ou 3)                    (fixe)
!   lmat(nnz1) : matrice de lissage sous forme sparse    (variable)
!   lic(nnz1)  : index colonnnes de la matrice de lissage(variable)  
!   lid(nl1+1) : index lignes de la matrice de lissage   (variable)
!   lres(nl1)  : second membre du systme (=0. ici)       (variable)
!   ondesS   : =.true. si on utilise les ondes S         (fixe)
!   lx,ly,lz : coefficients de lissage                   (fixe)
!   nix,niy,niz :taille du modele                        (fixe) 
!   nnz : taille de la matrice sparse des derivees       (fixe)   
!           
!   nl  : nombre de lignes de la matrice des derivees
!
!-----------------------------------------------------------------------------
!
!---objet---------------------------------------------------------------------
!
! construction de la matrice sparse de lissage: (lmat,lic,lid,lres)
! modification de nnz1 et nl1 
! 
!  
!###############################################################################
subroutine sub_lap(lap,lmat,lic,lid,lres,nnz1,nl1,   &
                   lx,ly,lz,nnz,nl,nix,niy,niz,ondesS)
implicit none
integer(kind=4) :: lap,nnz1,nl1,nnz,nl,nix,niy,niz,ligne,kpt
real(kind=4) :: lx,ly,lz
real(kind=4), dimension(:) :: lmat(:),lres(:)  ! lmat(nnz1),lres(nl1)
integer(kind=4), dimension(:) :: lic(:),lid(:) ! lic(nnz1),lid(nl1+1)
logical ondesS
!
ligne=0 ! ligne courante de la matrice de lissage
kpt=0   ! parametre courant de la matrice de lissage
!
!------------- pour debug -------------------------------
!     calcul de la matrice de lissage

if(lap == 1) then ! laplacien(x,y,z) =0 
!
  call lissage_1(lmat,lic,lid,lres,nnz1,nl1,   &
                  lx,ly,lz,nix,niy,niz,nnz,nl,0,ligne,kpt)
  if(ondesS) then
    call lissage_1(lmat,lic,lid,lres,nnz1,nl1,           &
                   lx,ly,lz,nix,niy,niz,nnz,nl,nix*niy*niz,ligne,kpt)
  endif
!
elseif(lap.eq.2) then ! laplacien (x,y) =0, laplacien (z) =0
  call lissage_2(lmat,lic,lid,lres,nnz1,nl1,        &
                         lx,ly,lz,nix,niy,niz,nnz,nl,0,ligne,kpt)
  if(ondesS) then
    call lissage_2(lmat,lic,lid,lres,nnz1,nl1,     &
                           lx,ly,lz,nix,niy,niz,nnz,nl,nix*niy*niz,ligne,kpt)
  endif
!
elseif(lap.eq.3) then  ! laplacien (x) =0, laplacien (y) =0, laplacien (z) =0
  call lissage_3(lmat,lic,lid,lres,nnz1,nl1,        &
                           lx,ly,lz,nix,niy,niz,nnz,nl,0,ligne,kpt)
  if(ondesS) then
    call lissage_3(lmat,lic,lid,lres,nnz1,nl1,     &
                           lx,ly,lz,nix,niy,niz,nnz,nl,nix*niy*niz,ligne,kpt) 
  endif
!
end if
!          
nl1=ligne  ! nombre de lignes de la matrice de lissage
nnz1=kpt   ! nombre d'elements de la matrice de lissage 
lid(nl1+1)=nnz+nnz1+1 ! dernier index des lignes de la matrice de lissage

return
end subroutine sub_lap 

!
!###################################################################################
!  construction de la matrice de lissage de type:
!
!                 lx*lap(x) + ly*lap(y) + lz*lap(z) = 0
!
!
! variables fixes :
! 
!         nix,niy,niz  : taille du modele
!         nnz          : taille de la matrice sparse des derivees 
!         nl           : nb de lignes de la matrice sparse des derivees (pour dedugage) 
!         lx, ly, lz   : ponderation du laplacien
!         vito         : =0 pour les ondes P, =nix*niy*niz pour les ondes S. 
!         nnz1         : taille de la martice de lissage
!         nl1          : nombre de lignes de la matrice de lissage
!
! variables non fixes :
!
!  l : ligne courante de la  matrice
!  k : element courant de la matrice
!  lmat(nnz1) : matrice sparse de lissage 
!  lic(nnz1)  : index colonnes de la matrice de lissage
!  lid(nl1+1) : index ligne de la matrice de lissage
!  lres(nl1)  : second membre (=0.)
!
subroutine lissage_1(lmat,lic,lid,lres,nnz1,nl1,lx,ly,lz,nix,niy,niz,nnz,nl,vito,l,k)
!
integer(kind=4) :: nnz1,nl1,ii,jj,kk,nnz,nl,nix,niy,niz,k,l,vito
real(kind=4) :: lx,ly,lz
real(kind=4), dimension(:) :: lmat(:),lres(:)    !  lmat(nnz1),lres(nl1)
integer(kind=4), dimension(:) :: lic(:),lid(:)   !  lic(nnz1),lid(nl1+1)

integer(kind=4) :: nl_0,nl1_0,nnz1_0
nl_0=nl; nl1_0=nl1; nnz1_0=nnz1  ! just for avoiding warnings from compilers

!     pour l'instant aucune contrainte sur la frontiere du domaine
!     boucle du le domaine de calcul :
do kk=2,niz-1
  do jj=2,niy-1
    do ii=2,nix-1
!     laplacien
      l=l+1          ! ligne suivante de la matrice
      lres(l)=0.     ! second membre de la ligne
      lid(l)=k+nnz+1 ! index de la ligne
!               write(*,*) l+nl,l,lid(l)
!     CONSTRUCTION  DE LA LIGNE :
!     1 .node (ii,jj,kk)
      k=k+1                                    ! element suivant de la matrice
      lmat(k)=-2.*lx-2.*ly-2.*lz               ! valeur de lissage
      lic(k)=ii+nix*(jj-1)+nix*niy*(kk-1)+vito ! index colonnne
!     2. node (ii+1,jj,kk)
      if(ii < nix) then
        k=k+1
        lmat(k)=lx
        lic(k)=ii+1+nix*(jj-1)+nix*niy*(kk-1)+vito
      endif
!     3. node (ii-1,jj,kk)
      if(ii > 1) then
        k=k+1
        lmat(k)=lx
        lic(k)=ii-1+nix*(jj-1)+nix*niy*(kk-1)+vito
      endif
!     4. node (ii,jj+1,kk)
      if(jj < niy) then
        k=k+1
        lmat(k)=ly
        lic(k)=ii+nix*(jj)+nix*niy*(kk-1)+vito
      endif
!     5. node (ii,jj-1,kk)
      if(jj > 1) then
        k=k+1
        lmat(k)=ly
        lic(k)=ii+nix*(jj-2)+nix*niy*(kk-1)+vito
      endif
!     6. node (ii,jj,kk+1)
      if(kk < niz) then
        k=k+1
        lmat(k)=lz
        lic(k)=ii+nix*(jj-1)+nix*niy*kk+vito
      endif
!     7. node (ii,jj,kk-1)
      if(kk > 1) then
        k=k+1
        lmat(k)=lz
        lic(k)=ii+nix*(jj-1)+nix*niy*(kk-2)+vito
      endif
!
!     FIN DE CONSTRUCTION DE LA LIGNE
!
    enddo
  enddo
enddo

return
end subroutine lissage_1
!
!##################################################################################
!
!###################################################################################
!  construction de la matrice de lissage de type:
!
!                 lx*lap(x) + ly*lap(y) = 0 
!                 lz*lap(z) = 0
!
!
! variables fixes :
! 
!         nix,niy,niz  : taille du modele
!         nnz          : taille de la matrice sparse des derivees 
!         nl           : nb de lignes de la matrice sparse des derivees (pour dedugage) 
!         lx, ly, lz   : ponderation du laplacien
!         vito         : =0 pour les ondes P, =nix*niy*niz pour les ondes S. 
!         nnz1         : taille de la martice de lissage
!         nl1          : nombre de lignes de la matrice de lissage
!
! variables non fixes :
!
!  l : ligne courante de la  matrice
!  k : element courant de la matrice
!  lmat(nnz1) : matrice sparse de lissage 
!  lic(nnz1)  : index colonnes de la matrice de lissage
!  lid(nl1+1) : index ligne de la matrice de lissage
!  lres(nl1)  : second membre (=0.)
!
subroutine lissage_2(lmat,lic,lid,lres,nnz1,nl1,lx,ly,lz,nix,niy,niz,nnz,nl,vito,l,k)
!
integer(kind=4) :: nnz1,nl1,ii,jj,kk,nnz,nl,nix,niy,niz,k,l,vito
real(kind=4) :: lx,ly,lz
real(kind=4), dimension(:) :: lmat(:),lres(:)    !   lmat(nnz1),lres(nl1)
integer(kind=4), dimension(:) :: lic(:),lid(:)   !   lic(nnz1),lid(nl1+1)

integer(kind=4) :: nl_0,nl1_0,nnz1_0
nl_0=nl; nl1_0=nl1; nnz1_0=nnz1  ! just for avoiding warnings from compilers

!     pour l'instant aucunes contraintes sur la frontiere 
!     du domaine
do kk=2,niz-1
  do jj=2,niy-1
    do ii=2,nix-1
!
!     ligne laplacien horizontal
      l=l+1
      lres(l)=0.
      lid(l)=k+nnz+1
!     1 .node (ii,jj,kk)
      k=k+1
      lmat(k)=-2.*lx-2.*ly
      lic(k)=ii+nix*(jj-1)+nix*niy*(kk-1)+vito
!     2. node (ii+1,jj,kk)
      if(ii < nix) then
        k=k+1
        lmat(k)=lx
        lic(k)=ii+1+nix*(jj-1)+nix*niy*(kk-1)+vito
      end if
!     3. node (ii-1,jj,kk)
      if(ii > 1) then
        k=k+1
        lmat(k)=lx
        lic(k)=ii-1+nix*(jj-1)+nix*niy*(kk-1)+vito
      end if
!     4. node (ii,jj+1,kk)
      if(jj < niy) then
        k=k+1
        lmat(k)=ly
        lic(k)=ii+nix*(jj)+nix*niy*(kk-1)+vito
      endif
!     5. node (ii,jj-1,kk)
      if(jj > 1) then
        k=k+1
        lmat(k)=ly
        lic(k)=ii+nix*(jj-2)+nix*niy*(kk-1)+vito
      endif
!
!     ligne laplacien vertical
      l=l+1
      lres(l)=0.
      lid(l)=k+nnz+1
!     1 .node (ii,jj,kk)
      k=k+1
      lmat(k)=-2.*lz
      lic(k)=ii+nix*(jj-1)+nix*niy*(kk-1)+vito
!     6. node (ii,jj,kk+1)
      if(kk < niz) then
        k=k+1
        lmat(k)=lz
        lic(k)=ii+nix*(jj-1)+nix*niy*kk+vito
      endif
!     7. node (ii,jj,kk-1)
      if(kk > 1) then
        k=k+1
        lmat(k)=lz
        lic(k)=ii+nix*(jj-1)+nix*niy*(kk-2)
      endif
!
    enddo
  enddo
enddo

return
end subroutine lissage_2
!
!##################################################################################
!###################################################################################
!  construction de la matrice de lissage de type:
!
!                 lx*lap(x) = 0 
!                 ly*lap(y) = 0 
!                 lz*lap(z) = 0
!
!
! variables fixes :
! 
!         nix,niy,niz  : taille du modele
!         nnz          : taille de la matrice sparse des derivees 
!         nl           : nb de lignes de la matrice sparse des derivees (pour dedugage) 
!         lx, ly, lz   : ponderation du laplacien
!         vito         : =0 pour les ondes P, =nix*niy*niz pour les ondes S. 
!         nnz1         : taille de la martice de lissage
!         nl1          : nombre de lignes de la matrice de lissage
!
! variables non fixes :
!
!  l : ligne courante de la  matrice
!  k : element courant de la matrice
!  lmat(nnz1) : matrice sparse de lissage 
!  lic(nnz1)  : index colonnes de la matrice de lissage
!  lid(nl1+1) : index ligne de la matrice de lissage
!  lres(nl1)  : second membre (=0.)
!
subroutine lissage_3(lmat,lic,lid,lres,nnz1,nl1,lx,ly,lz,nix,niy,niz,nnz,nl,vito,l,k)
!
      integer*4 nnz1,nl1,ii,jj,kk,nnz,nl,nix,niy,niz,k,l,vito
      real*4 lx,ly,lz
real(kind=4), dimension(:) :: lmat(:),lres(:)    ! lmat(nnz1),lres(nl1)
integer(kind=4), dimension(:) :: lic(:),lid(:)   ! lic(nnz1),lid(nl1+1)

integer(kind=4) :: nl_0,nl1_0,nnz1_0
nl_0=nl; nl1_0=nl1; nnz1_0=nnz1  ! just for avoiding warnings from compilers

!     pour l'instant aucunes contraintes sur la frontiere 
!     du domaine
do kk=2,niz-1
  do jj=2,niy-1
    do ii=2,nix-1
!
!     ligne laplacien x
      l=l+1
      lres(l)=0.
      lid(l)=k+nnz+1
!               write(*,*) l+nl,l,lid(l)
!     1 .node (ii,jj,kk)
      k=k+1
      lmat(k)=-2.*lx
      lic(k)=ii+nix*(jj-1)+nix*niy*(kk-1)+vito
!     2. node (ii+1,jj,kk)
      if(ii < nix) then
        k=k+1
        lmat(k)=lx
        lic(k)=ii+1+nix*(jj-1)+nix*niy*(kk-1)+vito
      endif
!     3. node (ii-1,jj,kk)
      if(ii >  1) then
        k=k+1
        lmat(k)=lx
        lic(k)=ii-1+nix*(jj-1)+nix*niy*(kk-1)+vito
      endif
!
!     ligne laplacien y
      l=l+1
      lres(l)=0.
      lid(l)=k+nnz+1
!     1 .node (ii,jj,kk)
      k=k+1
      lmat(k)=-2.*ly
      lic(k)=ii+nix*(jj-1)+nix*niy*(kk-1)+vito
!     4. node (ii,jj+1,kk)
      if(jj < niy) then
        k=k+1
        lmat(k)=ly
        lic(k)=ii+nix*(jj)+nix*niy*(kk-1)+vito
      endif
!     5. node (ii,jj-1,kk)
      if(jj.gt.1) then
        k=k+1
        lmat(k)=ly
        lic(k)=ii+nix*(jj-2)+nix*niy*(kk-1)+vito
      endif
!
!     ligne laplacien en z
      l=l+1
      lres(l)=0.
      lid(l)=k+nnz+1
!     1 .node (ii,jj,kk)
      k=k+1
      lmat(k)=-2.*lz
      lic(k)=ii+nix*(jj-1)+nix*niy*(kk-1)+vito
!     6. node (ii,jj,kk+1)
      if(kk < niz) then
        k=k+1
        lmat(k)=lz
        lic(k)=ii+nix*(jj-1)+nix*niy*kk
      endif
!     7. node (ii,jj,kk-1)
      if(kk > 1) then
        k=k+1
        lmat(k)=lz
        lic(k)=ii+nix*(jj-1)+nix*niy*(kk-2)+vito
      endif
!
    enddo
  enddo
enddo
return
end subroutine lissage_3

END MODULE s_lap


