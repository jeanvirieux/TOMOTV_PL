MODULE s_store_line
contains

! ############################################################################
!     Store one derivative matrix line DIRECTLY TO FILES 
!
!     input :
!     nsrc  : nombre total de sources
!     neqks : nombre de sources inversibles (seismes)
!     nvit  : nombre de noeuds inversibles
!     ioff  : =0 pour rai P, =nvit pour rai S
!     rai   : coordonnees du rai courant (ray connected)
!
!     iray  : current ray   just for tracing values if needed (not used)
!
!     slow_cur  : slowness vector at the source
!     ki_pos_cur: position index sources = index of invertible sources (eqks)
!     ki_to_cur : time index sources = index of invertible sources (eqks/blasts)
!
!     nx,ny,nz : grille d'inversion
!     der      : slowness/velocity derivatives
!
!     output :   WRITE DIRECTLY TO FILES
!     nu       : if nu == 0, we must cancel the data because no derivatives
!
!
! if the line has no partial derivatives, cancel the increment kptl (NEW TEST)
!
! ############################################################################
subroutine store_line(der,nu,ki_pos_cur,ki_to_cur,slow_cur,rai,nx,ny,nz, &
                  ixo,iyo,izo,nix,niy,niz,dtx,dty,dtz,               &
                  nsrc,neqks,nvit,ioff,iray,umat,uic,uid,kptl,kptm,flog)
implicit none
integer(kind=4) :: iray,nx,ny,nz,iray0,nx0,ny0,nz0
integer(kind=4) :: ixo,iyo,izo,nix,niy,niz
integer(kind=4) :: nvit,ioff,indx
integer(kind=4) :: nsrc,neqks,nsrc0
integer(kind=4) :: umat,uic,uid,flog
integer(kind=4) :: kptl,kptm   ! increment of the unknowns

integer(kind=4) :: ki_pos_cur,ki_to_cur

real(kind=4) :: slow_cur

!real(kind=4), dimension(:) :: u(:)      ! alias der_line
!integer(kind=4), dimension(:) :: iu(:)  ! alias der_ig
integer(kind=4) :: nu

real(kind=4), dimension(:,:,:) :: der(:,:,:)

real(kind=4), dimension(:,:) :: rai(:,:)

integer(kind=4) :: i,j,k,ii
real(kind=4) :: dtx,dty,dtz,xnorme,valeur

real(kind=4) :: WATER_LEVEL  ! (in meter)

iray0=iray              ! just to avoid annoying warning from compiler
nsrc0=nsrc              ! idem
nx0=nx;ny0=ny;nz0=nz    ! idem ... used for debugging


WATER_LEVEL=0.

dtx=-999.;dty=-999.;dtz=-999.    ! setting slowness at the source ... correct values only for quakes

ii=0   ! index des elements non nuls de la ligne
indx=0 ! index noeuds 

! indexes @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
kptl=kptl+1     ! increment the index of the line (data)
!               id(kptl)=kptm+1 ! store in the id the starting of non-zero elements
write(uid,rec=kptl) kptm+1   ! it is id(kptl) @@@@@@@@@@@@@@@@@

if(nix > 0 .and. niy > 0 .and. niz > 0) then
!=======================   ON PARCOURT LA SOUS-GRILLE  IXO:NI1; IYO:NI2; IZO:NI3
  do k=1,niz-izo+1
    do j=1,niy-iyo+1
      do i=1,nix-ixo+1
!   ig=i+ixo-1  jg=j+iyo-1  kg=k+izo-1
!   indx=ig + nx*(jg-1)+nx*ny*(kg-1) 
        indx=i+ixo-1+nx*(j+iyo-2)+nx*ny*(k+izo-2)
        valeur=der(i,j,k)
        if(valeur > WATER_LEVEL) then
          ii=ii+1
!JEAN          write(*,*) iray,ii,indx,ioff,valeur
!JEAN          write(77,*) ' iray: ii,indx,ioff,valeur', iray,ii,indx,ioff,valeur
!          u(ii)=valeur
!          iu(ii)=indx+ioff ! index of the selected node where derivative is stored in u
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ write directly to the file
          kptm=kptm+1
!    write in a file now   mat(kptm)
          write(umat,rec=kptm) valeur
!    write in a file now  directly the global index in the full presentation of the sensitivity matrix
          write(uic,rec=kptm) indx+ioff    ! it is ic(kptm) @@@@@@@@@@@@@@@@@@
        endif
      enddo
    enddo
  enddo
endif
!     
!     derivees par rapport a la source
!

if(ki_pos_cur /= 0) then   ! inverted eqks    les eqks_out n'arrivent pas ici a cause ki_m=0
  dtx=rai(1,2)-rai(1,1)
  dty=rai(2,2)-rai(2,1)
  dtz=rai(3,2)-rai(3,1)
  xnorme=sngl(dsqrt(dble(dtx**2)+dble(dty**2)+dble(dtz**2)))
!     - vecteur lenteur au foyer
  dtx=-slow_cur*dtx/xnorme
  dty=-slow_cur*dty/xnorme
  dtz=-slow_cur*dtz/xnorme

!     dt/dx
  ii=ii+1
!  u(ii)=dtx
!  iu(ii)=nvit+3*(ki_pos_cur-1)+1  ! ki_pos_cur IS the index of the eqks 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ write directly to the file
  kptm=kptm+1
!    write in a file now   mat(kptm)
  write(umat,rec=kptm) dtx
!    write in a file now  directly the global index in the full presentation of the sensitivity matrix
  write(uic,rec=kptm) nvit+3*(ki_pos_cur-1)+1    ! it is ic(kptm) @@@@@@@@@@@@@@@@@@
!     dt/dy                       ! et donc il positionne les parametres a la bonne colonne
  ii=ii+1                         ! ki_pos doit donc etre bien incremente !!!
!  u(ii)=dty
!  iu(ii)=nvit+3*(ki_pos_cur-1)+2
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ write directly to the file
  kptm=kptm+1
!    write in a file now   mat(kptm)
  write(umat,rec=kptm) dty
!    write in a file now directly the global index in the full presentation of the sensitivity matrix
  write(uic,rec=kptm) nvit+3*(ki_pos_cur-1)+2    ! it is ic(kptm) @@@@@@@@@@@@@@@@@@

!     dt/dz
  ii=ii+1
!  u(ii)=dtz
!  iu(ii)=nvit+3*(ki_pos_cur-1)+3
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ write directly to the file
  kptm=kptm+1
!    write in a file now   mat(kptm)
  write(umat,rec=kptm) dtz
!    write in a file now directly the global index in the full presentation of the sensitivity matrix
  write(uic,rec=kptm) nvit+3*(ki_pos_cur-1)+3    ! it is ic(kptm) @@@@@@@@@@@@@@@@@@

!JEAN  write(*,*) 'pos',iu(ii-2),iu(ii-1),iu(ii)
!
endif
!
if(ki_to_cur /= 0) then   ! inverted eqks/blasts
!     dt/dto
  ii=ii+1
!  u(ii)=1.                        ! on doit aussi avoir le bon increment aussi
!  iu(ii)=nvit+3*neqks+ki_to_cur   ! ki_to_cur IS the index of the eqks/blasts 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ write directly to the file
  kptm=kptm+1
!    write in a file now   mat(kptm)
  write(umat,rec=kptm) 1.
!    write in a file now directly the global index in the full presentation of the sensitivity matrix
  write(uic,rec=kptm) nvit+3*neqks+ki_to_cur   ! it is ic(kptm) @@@@@@@@@@@@@@@@@@
!JEAN  write(*,*) 'to',iu(ii)
!
endif
!
nu=ii
!     cas ligne nulle
if(ii.eq.0) then
!  nu=1
!  u(nu)=0.
!  iu(nu)=1
  kptl=kptl-1      ! we cancel this data ! STRANGE
  write(flog,*) ' this data has no partial derivatives ! '
endif
!
!
return
end subroutine store_line

END MODULE s_store_line

