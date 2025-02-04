MODULE s_submain
contains

! ######################################################################
!
!     Calcul:
!
!        -du temps de parcours tps_ray 
!        -des derivees par rapport aux lenteurs (en fait les vitesses) (der), 
!        -de la longeur du rai (long_ray),
!     Input:
!         urai : file unit where to read ray coordinates
!         der  : dimension(nix,niy,niz) derivatives
!         vit  : dimension(nx,ny,nz) velocity
!       pt_ray : global index for ray file
!     long_ray : total length of the ray (through simple computation)
!                will be compared with cumulative length through dichotomy l_ray
!     x_orig,* : origin of the total grid
!     nx,*     : nbre of nodes in each direction
!     hx,*     : grid step in each direction
!     ixo,*    : index of the node for the inversion
!     nix,*    : nbre of nodes for the inversion (starting at ixo,*)
!     iray     : index of the ray in relation with the data we consider
!
! ######################################################################
subroutine submain(urai,der,vit,pt_ray,long_ray,tps_ray,coord_ray, &
                   x_orig,y_orig,z_orig,nx,ny,nz,hx,hy,hz, &
                   ixo,iyo,izo,nix,niy,niz,incr, &
                   iray,flog)


implicit none

! grids model and inversion 
integer(kind=4) :: nx,ny,nz,nix,niy,niz,ixo,iyo,izo
real(kind=4) :: x_orig,y_orig,z_orig,hx,hy,hz
!real(kind=4), dimension(:,:,:) :: vit(:,:,:)     ! velocity
real(kind=4), dimension(:,:,:) :: vit(nx,ny,nz)     ! velocity
real(kind=4), dimension(:,:,:) :: der(:,:,:)  ! derivatives (in a reduced grid)
!                        en fait  der(nix+1-ixo,niy+1-iyo,niz+1-izo)
integer(kind=4) :: urai,iray
real(kind=4), dimension(:,:) :: coord_ray(:,:)
integer(kind=4), dimension(:) :: pt_ray(:)
integer(kind=4) :: npt,itab
real(kind=4) :: long_ray,tps_ray,l_ray
real(kind=4) :: l_tol


real(kind=8) :: XRO,YRO,ZRO,xo,yo,zo    ! c'est la lettre o ou O
real(kind=8) :: XRE,YRE,ZRE,xe,ye,ze

integer(kind=4) :: i,flog,INCR
logical :: test

l_tol=1.   ! tolerance in meters
l_ray=0.
tps_ray=0.
long_ray=0.
INCR=0

!
!     lecture du rai et calcul de la longeur du rai.
!
itab = pt_ray(iray)                 ! global index for the first element of the ray record  
npt = pt_ray(iray+1)-pt_ray(iray)   ! number of points for this ray

! JEAN write(*,*) 'npt,itab ',npt,itab    ! debug

call subreadrays(urai,coord_ray,npt,itab,long_ray) ! read it and store it into coord_ray

if(npt == 1) then                         ! debug JEAN
write(*,*) npt,itab,coord_ray(1,1),coord_ray(2,1),coord_ray(3,1) 
endif

!
!     Boucle sur les segments [X(i),X(i+1)] du rai
!

do i=1,npt-1
!        calcul des derivees et du temps de parcours en utilisant une dichotomie.
  XRO=dble(coord_ray(1,i))
  YRO=dble(coord_ray(2,i))
  ZRO=dble(coord_ray(3,i))
  XRE=dble(coord_ray(1,i+1))
  YRE=dble(coord_ray(2,i+1))
  ZRE=dble(coord_ray(3,i+1))

  call sub_test_seg(XRO,YRO,ZRO,XRE,YRE,ZRE,   &
                    x_orig,y_orig,z_orig,hx,hy,hz,test)
!
  DO WHILE (TEST)
!
!           le segment de rai traverse la frontiere de la cellule:
!                           |
!                X(i)       |           X(i+1)
!                +----------|-----------+ 
!                           |
!           cellule 1       |   hors cellule 1
!                           |       
!
    call sub_dicho(XRO,YRO,ZRO,XRE,YRE,ZRE,xo,yo,zo,xe,ye,ze,   &
                   x_orig,y_orig,z_orig,hx,hy,hz)
!
!
!              apres la dichotomie (sub_dicho) on se trouve dans la situation:
!     
!              |                  |                |                |
!              |    X(i)      Xo  | Xe             |       X(i+1)   |
!              |     +---------+--|-+--------------+-------+        |     
!              |                  |                |                |
!              |   cellule 1      |   cellule 2    | eventuelement  |
!              |                  |                | hors cellule 2 |
!
!
!              on integre sur [X(i),Xo] a l'interieur de la cellule 1
    call sub_int_l(tps_ray,vit,XRO,YRO,ZRO,xo,yo,zo, &
                   x_orig,y_orig,z_orig,hx,hy,hz,nx,ny,nz,l_ray)
!           calcul des derivee le long du segment [X(i),Xo]
    if(nix.gt.0.and.niy.gt.0.and.niz.gt.0) then
      call sub_int_dl(der,XRO,YRO,ZRO,xo,yo,zo,INCR,  &
                  x_orig,y_orig,z_orig,ixo,iyo,izo,  &
                  hx,hy,hz,nix,niy,niz)
    end if
    XRO=XE
    YRO=YE
    ZRO=ZE
    CALL SUB_TEST_SEG(XRO,YRO,ZRO,XRE,YRE,ZRE,    &
                  x_orig,y_orig,z_orig,hx,hy,hz,test)
       
  ENDDO ! on while(test)
!
!            Le segment de rai ne traverse pas la frontiere de la cellule:
!
!           |                                      |
!           |     X(i)                 X(i+1)      |
!           |      +---------------------+         |
!           |                                      |
!           |            cellule 1                 |
!           |                                      |
!     
!           Calcul du temps de parcours de X(i) a X(i+1).
  call sub_int_l(tps_ray,vit,xro,yro,zro,xre,yre,zre,    &
                     x_orig,y_orig,z_orig,hx,hy,hz,nx,ny,nz,l_ray)
!           Calcul des derivees par rapport a la lenteur.
  if(nix.gt.0.and.niy.gt.0.and.niz.gt.0) then
    call sub_int_dl(der,XRO,YRO,ZRO,XRE,YRE,ZRE,INCR,    &
                  x_orig,y_orig,z_orig,ixo,iyo,izo,  &
                  hx,hy,hz,nix,niy,niz)
  endif
end do    ! on the ray index i

!JEANwrite(*,*) ' number of parametres for this ray iray,INCR', iray,INCR

!
!     Changement de variables lenteurs -> vitesse    NOT HERE 
!
!if(nix.gt.0.and.niy.gt.0.and.niz.gt.0) then
!  call sub_multiplie(der,vit,nx,ny,nz,nix,niy,niz,ixo,iyo,izo)
!endif
!
if(abs(l_ray-long_ray).gt.l_tol) then
  write(flog,*) 'Dans submain, la dichotomie elimine une longueur de rai '
  write(flog,*) 'superieure a ',l_tol,' metres'
  write(flog,*) 'augmenter niter svp'
  write(flog,*) 'difference actuel entre les deux longueurs',l_ray-long_ray
  write(*,*) ' iray ',iray, 'missing length',l_ray-long_ray
endif
return
end subroutine submain

! #######################################################################
!    from slowness derivative to velocity derivative
! #######################################################################
subroutine sub_multiplie(der,vit,nx,ny,nz,nix,niy,niz,ixo,iyo,izo)
implicit none

integer(kind=4) :: nx,ny,nz,nix,niy,niz,ixo,iyo,izo
real(kind=4), dimension(:,:,:) :: der(nix,niy,niz)
real(kind=4), dimension(:,:,:) :: vit(nx,ny,nz)
integer(kind=4) :: i,j,k

! loops over the inverted grid 
do k=1,niz
  do j=1,niy
    do i=1,nix
      if(der(i,j,k).ne.0.) then
        der(i,j,k)=-der(i,j,k)/(vit(ixo+i-1,iyo+j-1,izo+k-1)**2)   ! derivative in velocity 
      endif
    enddo    ! i
  enddo      ! j
enddo        ! k
return
end subroutine sub_multiplie

! ######################################################################
!
!c     routine qui renvoie la variable logique test:
!c
!c     test=vrai si 
!c                   Xo          |          Xe
!c                   +-----------|---------+
!c                     cellule 1 | cellule 2
!c
!c     
!c     test=faux si
!c                   Xo                   Xe |
!c                   +--------------------+  |
!c                       cellule 1           |
!c
!c
!c     ou Xo=(xo,yo,zo) et Xe=(xe,ye,ze).
! ######################################################################
subroutine sub_test_seg(xo,yo,zo,xe,ye,ze,    &
           x_orig,y_orig,z_orig,hx,hy,hz,test)
implicit none
integer(kind=4) :: d_xo,d_yo,d_zo
integer(kind=4) :: d_xe,d_ye,d_ze
real(kind=4) :: hx,hy,hz
real(kind=4) :: x_orig,y_orig,z_orig
real(kind=8) :: xo,yo,zo
real(kind=8) :: xe,ye,ze
logical :: test
!
d_xo=int( ( ( xo - dble(x_orig) ) / dble(hx) )  )  +  1
d_yo=int( ( ( yo - dble(y_orig) ) / dble(hy) )  )  +  1
d_zo=int( ( ( zo - dble(z_orig) ) / dble(hz) )  )  +  1
d_xe=int( ( ( xe - dble(x_orig) ) / dble(hx) )  )  +  1
d_ye=int( ( ( ye - dble(y_orig) ) / dble(hy) )  )  +  1
d_ze=int( ( ( ze - dble(z_orig) ) / dble(hz) )  )  +  1
!     test pour savoir si le segment traverse la frontiere d'une cellule
! same cell ?
if((d_xo-d_xe).eq.0.and.(d_yo-d_ye).eq.0.and.(d_zo-d_ze).eq.0) then
  test=.false.
else
  test=.true.
endif
return
end subroutine sub_test_seg

! ########################################################################
! reading one ray in the big file :    itab is the offset
!             we compute the length at the same time 
! ########################################################################
subroutine subreadrays(urai,rai,npt,itab,long_rai)
implicit none
integer(kind=4) :: itab,npt,i,urai
real(kind=4), dimension(:,:)  :: rai(:,:)
real(kind=4) :: long_rai,dl 

! reading first point
read(urai,rec=itab+1-1) rai(1,1),rai(2,1),rai(3,1)

do i=2,npt
  read(urai,rec=itab+i-1) rai(1,i),rai(2,i),rai(3,i)
  dl=sngl(dsqrt(dble(rai(1,i)-rai(1,i-1))**2+dble(rai(2,i)-rai(2,i-1))**2+dble(rai(3,i)-rai(3,i-1))**2))
  long_rai=long_rai+dl
enddo

return
end subroutine subreadrays

! #######################################################################
!     Le segment de rai coupe la frontiere de la cellule, on veut 
!     encadrer ce point avec une precision de dl/2**niter(dl longeur
!     du segment de rai). La methode suivie est une simple dichotomie.
! #######################################################################
subroutine sub_dicho(xro,yro,zro,xre,yre,zre,xo,yo,zo,xe,ye,ze, &
                     x_orig,y_orig,z_orig,hx,hy,hz)
implicit none
integer(kind=4) :: i,niter
real(kind=8) :: ratio
real(kind=8) :: dl,dls
real(kind=8) :: xo,yo,zo
real(kind=8) :: xm,ym,zm
real(kind=8) :: xe,ye,ze
real(kind=8) :: xro,yro,zro
real(kind=8) :: xre,yre,zre
real(kind=4) :: x_orig,y_orig,z_orig
real(kind=4) :: hx,hy,hz
logical test
!
!     Segment de rai : [Xo(xro,yro,zro),Xe(xre,yre,zre)]
!     Point I a encadrer.
!
!             |-------x---------------|
!            Xo       I               Xe
!
niter=20
xo=xro
yo=yro
zo=zro
xe=xre
ye=yre
ze=zre
ratio=1.
dl=dsqrt((xe-xo)**2+(ye-yo)**2+(ze-zo)**2)
do i=1,niter
  xm=(xo+xe)*0.5d0
  ym=(yo+ye)*0.5d0
  zm=(zo+ze)*0.5d0
  call sub_test_seg(xo,yo,zo,xm,ym,zm,x_orig,y_orig,z_orig,hx,hy,hz,test)
  if(test) then
!                    I
!             |------x----*-----------|
!            Xo          Xm           Xe
    xe=xm
    ye=ym
    ze=zm
  else
!                               I
!             |-----------*-----x-----|
!            Xo          Xm           Xe
    xo=xm
    yo=ym
    zo=zm
  endif
! L'intervalle encadrant I a ete reduit, on continue ou on arrete 
! si on estime avoir I avec une precision suffisante.
  dls=dsqrt((xe-xo)**2+(ye-yo)**2+(ze-zo)**2)
  ratio=dl/dls
enddo
return
end subroutine sub_dicho

! #######################################################################
!
!     Calcul du temps de parcours le long du rai
!     pas les derivees
!
! #######################################################################
subroutine sub_int_l(tps,vit,xo,yo,zo,xe,ye,ze,      &
                     x_orig,y_orig,z_orig,hx,hy,hz,nx,ny,nz,l_rai)
use s_interpol3ds

implicit none
integer(kind=4) :: nx,ny,nz
real(kind=4), dimension(:,:) :: vit(nx,ny,nz)
real(kind=4) :: tps
integer(kind=4) :: d_xo,d_yo,d_zo
real(kind=4) :: hx,hy,hz
real(kind=4) :: x_orig,y_orig,z_orig
real(kind=8) :: xo,yo,zo,xxo,yyo,zzo
real(kind=8) :: xxm,yym,zzm
real(kind=8) :: xe,ye,ze,xxe,yye,zze
real(kind=8) :: dl
real(kind=4) :: l_rai
real(kind=4) :: lenteur1,lenteur2,lenteur3
real(kind=4) :: l1,l2,l3,l4,l5,l6,l7,l8
real(kind=8) :: xcur,ycur,zcur
! indices in the global grid     
      d_xo=int(( xo - dble(x_orig) ) / dble(hx))  +  1
      d_yo=int(( yo - dble(y_orig) ) / dble(hy))  +  1
      d_zo=int(( zo - dble(z_orig) ) / dble(hz))  +  1
! corner of the grid in the local frame (without the origine)
      xcur=dble((d_xo-1))*dble(hx)
      ycur=dble((d_yo-1))*dble(hy)
      zcur=dble((d_zo-1))*dble(hz)
! coordinates of both points     
      xxo=( (xo-dble(x_orig)) )
      yyo=( (yo-dble(y_orig)) )
      zzo=( (zo-dble(z_orig)) )
      xxe=( (xe-dble(x_orig)) )
      yye=( (ye-dble(y_orig)) )
      zze=( (ze-dble(z_orig)) ) 
!     Lenteur au 8 noeuds de la cellule 
      if(d_xo < 1 .or. d_xo > nx-1) then
        write(*,*) ' ray is going outside the box in x direction', d_xo,' above 1 and below ',nx-1
        write(*,*) ' FAILURE : please check the box dimension '
        stop
      endif
      if(d_yo < 1 .or. d_yo > ny-1) then
        write(*,*) ' ray is going outside the box in y direction', d_yo,' above 1 and below ',ny-1
        write(*,*) ' FAILURE : please check the box dimension '
        stop
      endif
      if(d_zo < 1 .or. d_zo > nz-1) then
        write(*,*) ' ray is going outside the box in z direction', d_zo,' above 1 and below ',nz-1
        write(*,*) ' FAILURE : please check the box dimension '
        stop
      endif
      l1=1./vit(d_xo,  d_yo,  d_zo)
      l2=1./vit(d_xo+1,d_yo,  d_zo)
      l3=1./vit(d_xo+1,d_yo+1,d_zo)
      l4=1./vit(d_xo,  d_yo+1,d_zo)
      l5=1./vit(d_xo,  d_yo,  d_zo+1)
      l6=1./vit(d_xo+1,d_yo,  d_zo+1)
      l7=1./vit(d_xo+1,d_yo+1,d_zo+1)
      l8=1./vit(d_xo,  d_yo+1,d_zo+1)
!      write(*,*) 'lenteurs',l1,l2,l3,l4,l5,l6,l7,l8
!     Interpolation 
      call interpol3ds(sngl(xcur),sngl(xcur)+hx,sngl(ycur),       &
           sngl(ycur)+hy,sngl(zcur),sngl(zcur)+hz,l1,l2,l3,      &
           l4,l5,l6,l7,l8,sngl(xxo),sngl(yyo),sngl(zzo),lenteur1)
      xxm=0.5d0*(xxo+xxe)
      yym=0.5d0*(yyo+yye)
      zzm=0.5d0*(zzo+zze)
      call interpol3ds(sngl(xcur),sngl(xcur)+hx,sngl(ycur),       &
           sngl(ycur)+hy,sngl(zcur),sngl(zcur)+hz,l1,l2,l3,      &
           l4,l5,l6,l7,l8,sngl(xxm),sngl(yym),sngl(zzm),lenteur2)
      call interpol3ds(sngl(xcur),sngl(xcur)+hx,sngl(ycur),       &
           sngl(ycur)+hy,sngl(zcur),sngl(zcur)+hz,l1,l2,l3,      &
           l4,l5,l6,l7,l8,sngl(xxe),sngl(yye),sngl(zze),lenteur3)
!      write(*,*) 'lent',lenteur1,lenteur2,lenteur3
!      write(*,'(a6,3f15.6,a7,f17.4)') 'seg',xo,yo,zo,'vit',1./lenteur1
!      write(*,*) 'vit',1./lenteur1,1./lenteur2,1./lenteur3
      dl=dsqrt((xe-xo)**2+(ye-yo)**2+(ze-zo)**2)
!----------------------------------------------------------------------
!     INTEGRATION
!     Formule de Simpson dl*[ (1/3)*f(x1) + (4/3)*f(x2) + (1/3)*f(x3) ]
!     avec:
!     >|-----dl----|<   
!      +-----------+-----------+
!     x1         x2           x3
!
!----------------------------------------------------------------------
tps=tps+sngl(dl*0.5d0*(lenteur1*(1.d0/3.d0)+lenteur2*(4.d0/3.d0)+lenteur3*(1.d0/3.d0)))
!    longeur du rai
l_rai=l_rai+sngl(dl)
return 
end subroutine sub_int_l

! ######################################################################
!       derivee  
!
! ######################################################################
subroutine sub_int_dl(der,xo,yo,zo,xe,ye,ze,incr,    &
                      x_orig,y_orig,z_orig,ixo,iyo,izo,hx,hy,hz,nix,niy,niz)
implicit none
integer(kind=4) :: nix,niy,niz,ixo,iyo,izo
real(kind=4), dimension(:,:,:) :: der(nix,niy,niz)
integer(kind=4) :: d_xo,d_yo,d_zo
real(kind=4) :: hx,hy,hz
real(kind=4) :: x_orig,y_orig,z_orig
real(kind=8) :: xo,yo,zo,xxo,yyo,zzo
real(kind=8) :: xe,ye,ze,xxe,yye,zze
real(kind=8) :: xxm,yym,zzm
real(kind=8) :: dl
real(kind=8) :: dlenteur1,dlenteur2,dlenteur3
!
integer(kind=4) :: d_xo_s,d_yo_s,d_zo_s,incr
common/inversion/d_xo_s,d_yo_s,d_zo_s
external indexe            ! block data for having initial values of indexes

!
!            ^  CELLULE
!          Z |
!            | 
!            5-------8      
!            |\      |\
!            | \     | \
!            |  \    |  \
!            |   6---+---7
!            1---|---4---+---> Y
!             \  |    \  |
!           X  \ |     \ |
!               \|      \|
!                2-------3
!                 \
!                  \
!                   v
!     dl longeur du segment de rai
!      write(*,*) 'CALCUL LE LONG DU RAI'
dl=dsqrt((xe-xo)**2+(ye-yo)**2+(ze-zo)**2)
!     coordonnees du premier noeud dans le repere de la grille d'inversion
!      write(*,*) nix,niy,niz
d_xo=int(( xo - dble(x_orig) ) / dble(hx)) + 1 - ixo +  1
d_yo=int(( yo - dble(y_orig) ) / dble(hy)) + 1 - iyo +  1
d_zo=int(( zo - dble(z_orig) ) / dble(hz)) + 1 - izo +  1
if(d_xo /= d_xo_s .or. d_yo /= d_yo_s .or. d_zo /= d_zo_s) then
  incr=incr+1
  d_xo_s=d_xo
  d_yo_s=d_yo
  d_zo_s=d_zo     ! for avoiding increments
endif
!JEAN write(77,*) d_xo,d_yo,d_zo
!JEAN write(78,*) xo,yo,zo,dl,d_xo,d_yo,d_zo
!      write(*,*) (ixo-1)*hx,(iyo-1)*hy,(izo-1)*hz
!      write(*,*) int((xo-(ixo-1)*hx)/hx)+1,d_xo
!      write(*,*) int((yo-(iyo-1)*hy)/hy)+1,d_yo
!      write(*,*) int((zo-(izo-1)*hz)/hz)+1,d_zo
!      write(*,*) 'fin'
!      pause
!     coordonnees locales  dans la grille d'inversion
xxo=( (xo-(dble(x_orig)+dble(ixo-1)*dble(hx)))-dble(d_xo-1)*dble(hx) )/dble(hx)
!      write(*,*) 'xo',xo
!      write(*,*) 'x_orig',x_orig
!      write(*,*) 'ixo',ixo
!      write(*,*) 'hx',hx
!      write(*,*) 'd_xo',d_xo
!      write(*,*) 'xxo',xxo
!      write(*,*) '-----------------------'
yyo=( (yo-(dble(y_orig)+dble(iyo-1)*dble(hy)))-dble(d_yo-1)*dble(hy) )/dble(hy)
zzo=( (zo-(dble(z_orig)+dble(izo-1)*dble(hz)))-dble(d_zo-1)*dble(hz) )/dble(hz)
xxe=( (xe-(dble(x_orig)+dble(ixo-1)*dble(hx)))-dble(d_xo-1)*dble(hx) )/dble(hx)
yye=( (ye-(dble(y_orig)+dble(iyo-1)*dble(hy)))-dble(d_yo-1)*dble(hy) )/dble(hy)
zze=( (ze-(dble(z_orig)+dble(izo-1)*dble(hz)))-dble(d_zo-1)*dble(hz) )/dble(hz)

xxm=(xxo+xxe)*0.5d0
yym=(yyo+yye)*0.5d0
zzm=(zzo+zze)*0.5d0
!
!     CALCUL DES DERIVEES AUX 8 NOEUDS
!
!     Integration utilisant la formule de Simpson:
!     dl*[ (1/3)*f(x1) + (4/3)*f(x2) + (1/3)*f(x3) ]
!     avec:
!     >|-----dl----|<   
!      +-----------+-----------+
!     x1         x2           x3
!
!     noeud 1
dlenteur1=(1.d0-xxo)*(1.d0-yyo)*(1.d0-zzo)
dlenteur2=(1.d0-xxm)*(1.d0-yym)*(1.d0-zzm)
dlenteur3=(1.d0-xxe)*(1.d0-yye)*(1.d0-zze)
if(d_xo.le.nix.and.d_xo.ge.1.and. &
   d_yo.le.niy.and.d_yo.ge.1.and. &
   d_zo.le.niz.and.d_zo.ge.1) then
  der(d_xo,d_yo,d_zo)=der(d_xo,d_yo,d_zo)+ &
          sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0)  &
         +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
endif
!     noeud 2
if(d_xo.lt.nix.and.d_xo.ge.1.and.     &
   d_yo.le.niy.and.d_yo.ge.1.and.     &
   d_zo.le.niz.and.d_zo.ge.1) then
  dlenteur1=(xxo)*(1.d0-yyo)*(1.d0-zzo)
  dlenteur2=(xxm)*(1.d0-yym)*(1.d0-zzm)
  dlenteur3=(xxe)*(1.d0-yye)*(1.d0-zze)
  der(d_xo+1,d_yo,d_zo)=der(d_xo+1,d_yo,d_zo)+ &
            sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0)  &
           +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
endif
!     noeud 3
if(d_xo.lt.nix.and.d_xo.ge.1.and.     &
   d_yo.lt.niy.and.d_yo.ge.1.and.     &
   d_zo.le.niz.and.d_zo.ge.1) then
  dlenteur1=(xxo)*(yyo)*(1.d0-zzo)
  dlenteur2=(xxm)*(yym)*(1.d0-zzm)
  dlenteur3=(xxe)*(yye)*(1.d0-zze)
  der(d_xo+1,d_yo+1,d_zo)=der(d_xo+1,d_yo+1,d_zo)+   &
            sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0)     &
           +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
endif
!     noeud 4
if(d_xo.le.nix.and.d_xo.ge.1.and.     &
   d_yo.lt.niy.and.d_yo.ge.1.and.     &
   d_zo.le.niz.and.d_zo.ge.1) then
  dlenteur1=(1.d0-xxo)*(yyo)*(1.d0-zzo)
  dlenteur2=(1.d0-xxm)*(yym)*(1.d0-zzm)
  dlenteur3=(1.d0-xxe)*(yye)*(1.d0-zze)
  der(d_xo,d_yo+1,d_zo)=der(d_xo,d_yo+1,d_zo)+  &
           sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0) &
          +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
      end if
!     noeud 5
if(d_xo.le.nix.and.d_xo.ge.1.and.     &
   d_yo.le.niy.and.d_yo.ge.1.and.     &
   d_zo.lt.niz.and.d_zo.ge.1) then
  dlenteur1=(1.d0-xxo)*(1.d0-yyo)*(zzo)
  dlenteur2=(1.d0-xxm)*(1.d0-yym)*(zzm)
  dlenteur3=(1.d0-xxe)*(1.d0-yye)*(zze)
  der(d_xo,d_yo,d_zo+1)=der(d_xo,d_yo,d_zo+1)+   &
          sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0)   &
         +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
endif
!     noeud 6
if(d_xo.lt.nix.and.d_xo.ge.1.and.     &
   d_yo.le.niy.and.d_yo.ge.1.and.     &
   d_zo.lt.niz.and.d_zo.ge.1) then
  dlenteur1=(xxo)*(1.d0-yyo)*(zzo)
  dlenteur2=(xxm)*(1.d0-yym)*(zzm)
  dlenteur3=(xxe)*(1.d0-yye)*(zze)
  der(d_xo+1,d_yo,d_zo+1)=der(d_xo+1,d_yo,d_zo+1)+    &
            sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0)      &
           +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
endif
!     noeud 7
if(d_xo.lt.nix.and.d_xo.ge.1.and.        &
   d_yo.lt.niy.and.d_yo.ge.1.and.        &
   d_zo.lt.niz.and.d_zo.ge.1) then
  dlenteur1=(xxo)*(yyo)*(zzo)
  dlenteur2=(xxm)*(yym)*(zzm)
  dlenteur3=(xxe)*(yye)*(zze)
  der(d_xo+1,d_yo+1,d_zo+1)=der(d_xo+1,d_yo+1,d_zo+1)+ &
            sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0)       &
           +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
endif
!     noeud 8
if(d_xo.le.nix.and.d_xo.ge.1.and.        &
   d_yo.lt.niy.and.d_yo.ge.1.and.        &
   d_zo.lt.niz.and.d_zo.ge.1) then
  dlenteur1=(1.d0-xxo)*(yyo)*(zzo)
  dlenteur2=(1.d0-xxm)*(yym)*(zzm)
  dlenteur3=(1.d0-xxe)*(yye)*(zze)
  der(d_xo,d_yo+1,d_zo+1)=der(d_xo,d_yo+1,d_zo+1)+    &
          sngl(dl*0.5d0*(dlenteur1*(1.d0/3.d0)        &
         +dlenteur2*(4.d0/3.d0)+dlenteur3*(1.d0/3.d0)))
endif
!
return 
end subroutine sub_int_dl

END MODULE s_submain







