! ###########################################################################
!
! modification JANUARY 2015 .... for ORTENSIA ... ray may be bent strongly
! making the back raytracing moving away from the end point during a while
! how long depends on the heterogeneities ... no obvious rules ...
!
! ###########################################################################
!     Trace de rai: Calcul des trajets sources-stations et 
!     stockage des rais sur le disque.
!
!     input :  
!              x_src(nsrc) : coordonnees des sources
!              ki_m(nsrc)  : index des sources 
!              ki_id(nsrc) : ID of sources
!              nsrc        : nombre de sources  
!              x_sta(nsta) : coordonnees des stations
!              xg,yg,zg    : point origine du modele
!              n1,n2,n3    : taille du modele de vitesse = grille P&L
!              h1,h2,h3    : stepping of the grid (h1=h2=h3)    and rap is the ratio sampling for rays
!              nsta  : nombre de stations
!              hs(n1,n2,n3): modele P&L 
!              uray : fichier contenant les coordonnees des rais
!              utab : fichier contenant le pointeur sur le fichier des rais
!              nray : nombre de rais
!              utab : fichier contenant le pointeur sur le fichier des rais
!              uray : fichier contenant les coordonnees des rais
!              rmax : nombre de points max pour trace de rais 
!              kmax : estimated number of points for raytracing
!              noff : offset in the data (for P is zero; for S is ntp)
!              nt_cur : number of data (for P is ntp; for S is nts)
!
!     workspace : rai(3,rmax)  : coordonnnees des rais 
!                 t(n1,n2,n3)  : temps de parcours (calcule par P&L)
!
!     output/input : nray : nombre de rais    set to zero in the main program at the beginning
!                    kray_tot (alias kpt) : global index for the storage of ray coordinates
!
!     output : files of ray positions as well as pointer for a single file.
!              rays are identified by two indices (isrc,ista) always in this order
!
!     no management of ier at this level
!
!
!================================================ Janvier 2016
!
!================================================ april 2016
! delete lut_dat as not used ... does not harm anyhow
!
! take into account the fact that a source can be outside while there is still a time data
! for this source ... the related id_src of the data has a zero index ...
!
! ########################################################################
MODULE s_ray
contains

  subroutine submain(pri,x_src,y_src,z_src,ki_m,ki_id,nsrc,x_sta,y_sta,z_sta,id_sta,nsta, &
       xg,yg,zg,n1,n2,n3,h1,h2,h3,rap,            &
       lut_src,lut_sta,lut_ray,       &
       hs,t,utab,uray,kray_tot,nray,rmax,kmax,kmin,noff,nt_cur)
    implicit none

    integer(kind=4) :: pri,rmax,kmax,kmin,kray_tot,kr
    ! velocity and time 3D grids
    integer(kind=4) :: noff,nt_cur      ! P or S waves
    real(kind=4) :: hs(n1,n2,n3),t(n1,n2,n3)
    ! true source
    integer(kind=4), dimension(:) :: ki_m(:),ki_id(:)
    real(kind=4), dimension(:) :: x_src(:),y_src(:),z_src(:)
    integer(kind=4) :: nsrc,isrc,idsrc,idsrc_min,idsrc_max,ndsrc
    ! true station
    real(kind=4), dimension(:) :: x_sta(:),y_sta(:),z_sta(:)
    integer(kind=4) :: nsta,jsta,idsta,idsta_min,idsta_max,ndsta
    integer(kind=4), dimension(:) :: id_sta(:)
    ! data
    integer(kind=4), dimension(:) :: lut_src(:),lut_sta(:),lut_ray(:)
    integer(kind=4) :: idata
    ! ray coordinates
    integer(kind=4) :: nray,utab,uray
    real(kind=4), allocatable,dimension(:,:) :: rai
    ! grid
    integer(kind=4) :: n1,n2,n3
    real(kind=4) :: xg,yg,zg,h1,h2,h3,rap
    ! virtual source
    real(kind=4) :: xs,ys,zs

    ! flag for making the inverse relation between ID and index either for stations or sources
    integer(kind=4), allocatable,dimension(:) :: lut_inv

    integer(kind=4) :: itemps,time_3d
    integer(kind=4) :: msg,ier,iray_missing
    real(kind=4) :: eps
    !

    integer(kind=4) :: debug    ! debug = 0 (no output) 1 () 2 (output)

    !===================================== debugging if you need
    debug=0

    msg=0
    eps=0.001

    allocate(rai(3,rmax))
    kr=0
    iray_missing=0

    ! used for raytracing (saving on file will correct for that)
    !if(nsrc > nsta) then
    !  pri=1  ! priority station for raytracing    nsta < nsrc
    !  write(*,*) ' assuming eikonal solver from stations as virtual sources '
    !else
    !  pri=0  ! priority source for raytracing     nsrc < nsta
    !  write(*,*) ' assuming eikonal solver from sources as virtual sources '
    !endif
    ! DONE IN THE MAIN PROGRAM
    ! 
    !
    !  less number of stations than number of sources   double boucle sur les stations et les donnees
    !
    if(pri == 1) then
       !--------------------------
       idsrc_min=ki_id(1)        ! we must have the inverse relation    isrc=lut_inv(idsrc) au lieu de idsrc=lut_src(isrc)
       idsrc_max=ki_id(1)
       do isrc=2,nsrc
          idsrc_min=MIN(idsrc_min,ki_id(isrc))
          idsrc_max=MAX(idsrc_max,ki_id(isrc))
       enddo
       ndsrc=idsrc_max-idsrc_min+1  ! range of IDs     EVITER D'AVOIR N'IMPORTE QUOI COMME IDs
       write(*,*) ' PRIORITY ON STATIONS : range of IDs source ',ndsrc,idsrc_min,nsrc
       !=================== build the inverse table
       allocate(lut_inv(ndsrc))
       lut_inv(:)=0     ! si un ID manque
       do isrc=1,nsrc
          idsrc=ki_id(isrc)
          if(idsrc-idsrc_min+1 > ndsrc) then
             write(*,*) ' submain error in loop over lut_inv source ',idsrc,isrc,ndsrc,idsrc_min
          endif
          lut_inv(idsrc-idsrc_min+1)=isrc    ! la relation inverse pour trouver rapidement le compteur via l'ID
       enddo
       !=================== done
       !--------------------------   STATIONS CAR NSTA < NSRC 
       do jsta=1,nsta             ! boucle sur les stations
          idsta=id_sta(jsta)       ! ID OF THE STATION
          !           xs,ys,zs : source position in the P&L grid
          call ccoo(xg,yg,zg,n1,n2,n3,h1,h2,h3,x_sta(jsta),y_sta(jsta),z_sta(jsta),xs,ys,zs,ier)
          if(ier == 0) then ! ok the station is inside the grid
             if(mod(jsta,10).eq.0) then
                write(*,'(30x,a10,i12,a3,i12)') 'station ',jsta,' / ',nsta
             endif
             !           calcul de la table des temps (P&L)
             if(debug > 0) write(*,*) ' ray tracing submain P&L ',kr,rmax
             itemps=time_3d(hs,t,n1,n2,n3,xs-1,ys-1,zs-1,eps,msg)    ! be sure of 'minus 1' for P&L source 
             !           trace de rais (boucle sur les donnees maintenant)
             do idata=noff+1,noff+nt_cur
                if(lut_sta(idata) == idsta) then     ! on est bien a la station idsta pour cette donnee
                   idsrc=lut_src(idata)               ! on prend l'ID de la source given by the data
                   if(idsrc-idsrc_min+1 > ndsrc) then
                      write(*,*) ' submain error in loop over lut_inv ',idsrc,isrc,ndsrc
                   endif
                   isrc=lut_inv(idsrc-idsrc_min+1)    ! on prend le compteur actuel pour l'ID de la source
                   !JEAN write(*,*) ' TEST JEAN',idsrc,idsrc_min
                   if(isrc == 0) then
                      !write(*,*) ' error in RAYTRACING>SUBMAIN in LUT_INV zero index ',isrc,' for the id source ',idsrc
                      !write(*,*) ' FATAL ERROR : COULD NOT SURVIVE '
                      !                      stop
                      write(*,*) ' for this data ',idata,' the source is out ',idsrc
                      lut_ray(idata)=0
                      goto 888   !!! JEAN PATCH
                   endif
                   ! is the source inside the domain (check on the ki_m flag)
                   if(ki_m(isrc) /= 0) then
                      nray=nray+1
                      lut_ray(idata)=nray              ! on stocke le rayon pour la donnee idata
                      !==================================================== trace du rai nray nbre of points kr
                      call sub_rai(t,x_sta(jsta),y_sta(jsta),z_sta(jsta),   &
                           x_src(isrc),y_src(isrc),z_src(isrc),   &
                           rai,rmax,xg,yg,zg,n1,n2,n3,h1,h2,h3,rap,kr,ier)
                      !================================================================================ DEBUGGING ... related to rmax
                      !JEAN write(*,*) ' OK SUCCES ',x_sta(jsta),y_sta(jsta),z_sta(jsta),idsta
                      !JEAN write(*,*) ' OK SUCCES ',x_src(isrc),y_src(isrc),z_src(isrc),idsrc
                      !                      if(kr <= 10) then
                      !                      write(*,*) ' OK SUCCES ', kr
                      !                      write(*,*) x_sta(jsta),y_sta(jsta),z_sta(jsta),idsta
                      !                      write(*,*) x_src(isrc),y_src(isrc),z_src(isrc),idsrc
                      !                      write(*,*) ' ONE MORE RAY '
                      !                      endif
                      !==============================================================================================================
                      kmax=max(kmax,kr)    ! recherche du maximum
                      if(ier .ne. -2) then
                         kmin=min(kmin,kr)    ! recherche du minimum 
                      else
                         iray_missing=iray_missing+1
                      endif
                      !===================================  ecriture du rai sur le disque
                      !===================================  output is always from TRUE source to REAL station
                      call st_rai(pri,uray,utab,rai,rmax,kray_tot,kr,idsrc,idsta,nray)! always from source to station    
                   else 
                      lut_ray(idata)=0                   !STEPHANIE: on met le numero du rayon  zero quand la source est OUT.
                   endif   ! ki_m(isrc)
                   
888                continue !!!!! JEAN PATCH
                   
                endif   ! lut_sta   (ID de la station)
             enddo  ! idata sur noff+1,noff+nt_cur
          endif    ! on ier  for station inside the grid
       enddo      ! jsta on 1,nsta
       deallocate(lut_inv)

    else

       !
       ! less number of sources than number of stations     double boucle sur les sources et les donnees 
       !

       idsta_min=id_sta(1)
       idsta_max=id_sta(1)
       do jsta=2,nsta
          idsta_min=MIN(idsta_min,id_sta(jsta))
          idsta_max=MAX(idsta_max,id_sta(jsta))
       enddo
       ndsta=idsta_max-idsta_min+1  ! range of IDs     EVITER D'AVOIR N'IMPORTE QUOI COMME IDs
       write(*,*) ' for PRIORITY OF SOURCES : range of IDs for stations ',ndsta
       allocate(lut_inv(ndsta))
       lut_inv(:)=0
       do jsta=1,nsta
          idsta=id_sta(jsta)
          lut_inv(idsta-idsta_min+1)=jsta    ! la relation inverse pour trouver rapidement le compteur via l'ID
          if(idsta-idsta_min+1 > ndsta) then
             write(*,*) ' submain error in loop over lut_inv station ',idsta,jsta,ndsta
          endif
          !JEAN    write(*,*) jsta,idsta_min,idsta,idsta-idsta_min+1
       enddo
       do isrc=1,nsrc
          idsrc=ki_id(isrc)
          !JEAN     write(*,*) ' isrc,idsrc ',isrc,idsrc
          ! is the source inside the domain (check on the ki_m flag)
          if(ki_m(isrc) /= 0) then      ! only when the source is inside the global grid
             !         xs,ys,zs : source position in the P&L grid
             call ccoo(xg,yg,zg,n1,n2,n3,h1,h2,h3,x_src(isrc),y_src(isrc),z_src(isrc),xs,ys,zs,ier)  
             if(ier /= 0) then 
                write(*,*) ' RAYTRACING>SUBMAIN contradiction : the true source should be inside'
                write(*,*) ' and has been found outside in the ccoo subroutine STOP '
                stop
             endif
             if(mod(isrc,10).eq.0) then
                write(*,'(30x,a10,i12,a3,i12)') 'source ',isrc,' / ',nsrc
             endif
             !          calcul de la table des temps (P&L)
             if(debug > 0) write(*,*) ' ray tracing submain P&L ',kr,rmax
             itemps=time_3d(hs,t,n1,n2,n3,xs-1,ys-1,zs-1,eps,msg)    ! be sure of 'minus 1' for P&L source 
             !           trace de rais (boucle sur les donnees maintenant)
             !##################################################### we look over data to get those related to the source identified by idsrc
             !##################################################### might not be the best way ... 
             do idata=noff+1,noff+nt_cur
                !JEAN write(*,*) ' idata,lut_src,lut_sta,lut_ray ', idata,lut_src(idata),lut_sta(idata),lut_ray(idata)
                !JEAN write(*,*) ' idsrc, ki_id, isrc ',idsrc,ki_id(isrc),isrc
                kr=1
                if(lut_src(idata) == idsrc) then     ! on est bien a la source idsrc pour cette donnee
                   idsta=lut_sta(idata)               ! on prend l'ID de la station de cette donnee
                   if(idsta-idsta_min+1 > ndsta .OR. idsta-idsta_min+1 <= 0) then
                      write(*,*) ' submain error in loop over lut_inv ',idsta,jsta,ndsta,idsta_min
                   endif
                   jsta=lut_inv(idsta-idsta_min+1)    ! on prend le compteur actuel pour l'ID de la station
                   !JEAN write(*,*) ' TEST JEAN jsta,idsta ',jsta,idsta,idsta_min
                   if(jsta == 0) then
                      write(*,*) ' error in RAYTRACING>SUBMAIN in LUT_INV zero index for the id station ',jsta,idsta
                      write(*,*) ' FATAL ERROR : could not survive '
                      stop
                   endif
                   call ccoo(xg,yg,zg,n1,n2,n3,h1,h2,h3,x_sta(jsta),y_sta(jsta),z_sta(jsta),xs,ys,zs,ier)
                   if(ier == 0) then    ! station is inside the grid (should be the case)
                      nray=nray+1
                      lut_ray(idata)=nray              ! on stocke le numero du rayon pour la donnee idata
                      !==================================================== trace du rai nray nbre of points kr
                      call sub_rai(t,x_src(isrc),y_src(isrc),z_src(isrc),   &
                           x_sta(jsta),y_sta(jsta),z_sta(jsta),   &
                           rai,rmax,xg,yg,zg,n1,n2,n3,h1,h2,h3,rap,kr,ier)
                      kmax=max(kmax,kr)    ! recherche du maximum
                      if(ier .ne. -2) then
                         kmin=min(kmin,kr)    ! recherche du minimum 
                      else
                         iray_missing=iray_missing+1
                      endif
                      if(debug >0) write(*,*) ' OK SUCCES ',x_sta(jsta),y_sta(jsta),z_sta(jsta),jsta
                      if(debug >0) write(*,*) ' OK SUCCES ',x_src(isrc),y_src(isrc),z_src(isrc),isrc
                      !===================================  ecriture du rai sur le disque
                      !===================================  output is always from TRUE source to REAL station
                      call st_rai(pri,uray,utab,rai,rmax,kray_tot,kr,idsrc,idsta,nray)! always from source to station
                      ! JEAN write(*,*) ' tracing kr ',kr,rai(1,kr),rai(2,kr),rai(3,kr)
                   endif   ! on ier for station inside the grid
                endif   ! lut_src   (ID de la source)
             enddo  ! idata sur noff+1,noff+nt_cur        we have spanned the entire dataset ... et grab what is related to the source idsrc
          else 
             lut_ray(idata)=0                   !STEPHANIE: on met le numero du rayon  zero quand la source est OUT.
          endif    ! ki_m(isrc) related to the source isrc
       enddo       ! isrc on 1,nsrc
       deallocate(lut_inv)
    endif

    deallocate(rai)

    write(*,*) ' %%%%%%%%%%%%%%%%%%% missing rays ',iray_missing
    write(*,*) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    return
  end subroutine submain

  !##############################################################################
  !     Sauvegarde du rai courant sur le disque.
  !
  !     input :  pri (=1 => priorite station), (=0 => priorite source)
  !              uray    : fichier contenant les coordonnees des rais
  !              utab    : fichier contenant le pointeur sur le fichier des rais
  !              kr      : nombre de points du rai courant   local index
  !              kray_tot: pointeur de debut du rai courant dans le fichier uray  global index
  !              iray    : numero du rai courant
  !              idsrc   : ID OF THE SOURCE
  !              idsta   : ID OF THE STATION
  !
  !     output : kray_tot: pointeur de debut du rai suivant dans le fichier uray
  !
  subroutine st_rai(pri,uray,utab,rai,rmax,kray_tot,kr,idsrc,idsta,iray)
    implicit none
    integer(kind=4) :: ir,pri,uray,utab,kray_tot,kr,idsrc,idsta,iray,rmax
    real(kind=4), dimension(:,:) :: rai(3,rmax)
    !     ecriture pointeur rai
    if(iray == 1) then     ! first ray
       write(utab,rec=iray+1) 1,idsrc,idsta   ! si iray=1, on a kray_tot=0 ... d'ou drole de test ??
    else
       write(utab,rec=iray+1) kray_tot+1,idsrc,idsta !we write the global index kray_tot+1 
       !(just where is the first point of the ray)
    endif                          ! in the index file
    !     ecriture du rai

    if(pri == 1) then
       do ir=1,kr    ! number of points of this ray!            from the true source towards the station : normal backraytracing
          write(uray,rec=kray_tot+ir)  rai(1,ir),rai(2,ir),rai(3,ir)
       enddo
    else
       do ir=1,kr    ! number of points of this ray!            from the true source towards the station : normal backraytracing
          !          from the station towards the true source : do inverse backraytracing
          write(uray,rec=kray_tot+ir) rai(1,kr-ir+1),rai(2,kr-ir+1),rai(3,kr-ir+1)
       enddo
    end if
    !
    kray_tot=kray_tot+kr        ! global index
    !JEAN write(*,*) ' kray_tot ',kray_tot,kr
    !
    return
  end subroutine st_rai
  !
  !##############################################################################
  !     Tracing one ray
  !
  !     input : t(n1,n2,n2): table des temps de propagation (calcul par P&L)
  !             xo,yo,zo   : origine du rai
  !             xe,ye,ze   : extremite du rai
  !             rai(3,rmax): coordonnees des pts du rayon
  !             rmax       : nombre max de points sur le rai
  !
  !             xg,yg,zg   : origine de la grille P&L
  !             n1,n2,n3   : grille P&L
  !             h1,h2,h3   : pas dans les trois directions h1=h2=h3 en general
  !             kr         : index courant sur les points du rai

  !     output : kr : nombre de point sur le rai trace
  !              rai(3,kr) : coordonnees du rai trace (seulement kr elements utilises)
  !
  ! =============================================================================
  !  on trace le rai depuis l'extremite jusqu'a l'origine (point emettant les 
  !  ondes dans le calcul P&L) en descendant le long de la plus grande pente
  !  de la fonction t (cf equation de l'Eikonale). 
  ! #############################################################################
  !
  subroutine sub_rai(t,xo,yo,zo,xe,ye,ze,rai,rmax,xg,yg,zg,n1,n2,n3,h1,h2,h3,rap,kr,ier)
    implicit none
    ! grid
    integer(kind=4) :: n1,n2,n3
    real(kind=4) xg,yg,zg,h1,h2,h3,rap     ! rap is the ratio sampling with respect to h
    real(kind=4), dimension(:,:,:) :: t(n1,n2,n3)
    ! ray
    integer(kind=4) :: rmax
    real(kind=4) :: xo,yo,zo,xe,ye,ze,eps,h,d1,d2,d3,xr,yr,zr,xinorm
    real(kind=4), dimension(:,:) :: rai(3,rmax)

    integer(kind=4) :: kr,i,ismooth,ier,iter,iter_tol,iter_fail,iflag_38
    real(kind=4) :: dist_old_old,dist_old,dist,attraction
    real(kind=4) :: off1,off2,off3,xl
    data ismooth/5/     ! kr must be at least equal to 5
    data iflag_38/1/
    integer(kind=4) :: debug    ! 0 or 1 or 2
    save iflag_38

    !=============================== setup debug to higher values for verbose output
    debug=0

    !h=(h1*h2*h2)**0.3333
    !h=h1*rap
    h=rap     ! better absolute value rather than a ratio with respect to grid step

    eps=2*h    !value to be tuned ... no systematic checking ....   or 2*h1
!    eps=2*h1   ! to be related to the grid step
    !     extremite du rai : stockage du premier point
    kr=1
    rai(1,kr)=xe
    rai(2,kr)=ye
    rai(3,kr)=ze

    if(debug >0) write(*,*) ' STARTING POINT xe,ye,ze ',xe,ye,ze  !JEAN
    if(debug >0) write(*,*) ' ENDING POINT   xo,yo,zo ',xo,yo,zo  !JEAN

    !     changement de coordonnees  vers les positions dans la grille PL xy,yr,zr
    xr=0.;yr=0.;zr=0.
    call ccoo(xg,yg,zg,n1,n2,n3,h1,h2,h3,xe,ye,ze,xr,yr,zr,ier)
    if(ier.ne.0) then
       write(*,*) ' fatal error in sub_rai at the initial point = should not be as tested in submain '
       write(*,*) ' ier,rai_x,rai_y,rai_z ',ier,xe,ye,ze
       stop
    endif
    !     calcul du segment de rai : step length to be done in the deepest descent
    call seg_rai(xr,yr,zr,t,n1,n2,n3,d1,d2,d3)

    !     boucle sur les segments de rais

    dist=sngl(dsqrt(dble((rai(1,kr)-xo)**2)+dble((rai(2,kr)-yo)**2)+dble((rai(3,kr)-zo)**2)))

    dist_old=dist
    iter=0
    iter_tol=30               ! could be modified depending on the model
    iter_fail=60              ! could be modified dependins on the model
    attraction=10.            ! zone d'attraction pour test d'approche

    if(iflag_38 == 1) then    ! write only once this information
       iflag_38=0
       write(*,*) ' RAY tracing between source/receiver: adapt the tuning parameter in model.head '
       write(*,*) ' sampling step along ray '
       write(*,*) ' other parameters inside the subroutine submain.f90 in the directory RAYTRACING '
       write(*,*) ' attraction zone distance (meters) (related to the forward grid size)',attraction*eps
       write(*,*) ' iter_tol: number of steps moving away from the end point above warning writing',iter_tol
       write(*,*) ' iter_fail: number of steps moving away from the end point above stopping the ray tracing',iter_fail
       write(*,*) ' it is somehow connected with the subsampling along the ray ... see inversion.head file in the directory DATA '
    endif

    do while (dist >  eps)
       !------------------------- buffer de trois distances
       dist_old_old=dist_old
       dist_old=dist
       dist=sngl(dsqrt(dble((rai(1,kr)-xo)**2)+dble((rai(2,kr)-yo)**2)+dble((rai(3,kr)-zo)**2)))

!=====================================================================
       ! checking the distance between the current point on the ray and the receiver position
       ! should be performed only in the vicinity of the receiver  10 times eps
       !=====================================================================
       if(dist_old_old < dist_old .or. dist_old < dist .or. dist_old_old < dist .and. iter > iter_tol) then
         if(dist < attraction*eps) then
             !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ arriving at the attraction distance of the receiver
             write(*,*) ' WARNING : 2 PTS RAY TRACING convergence problem RAYTRACING>SUBMAIN'
             write(*,*) ' ray starts moving away the receiver ... '
             write(*,*) ' already ',iter_tol,' pts have been identified ... '
             write(*,*) ' dist A > B > C ',dist_old_old,dist_old,dist,eps
             iter=iter+1
             if(iter > iter_fail) then
                write(*,*) ' ======================== '
                write(*,*) ' KO FAILURE',xe,ye,ze, 'starting'  !JEAN
                write(*,*) ' KO FAILURE',xo,yo,zo, 'ending'    !JEAN
                write(*,*) ' Warning: we remove this ray for the tomography as the ray goes away the end point'
                write(*,*) ' ======================== '
                kr=1    ! erase computed values in order to cancel this ray
                rai(1,kr)=-999.;rai(2,kr)=-999.;rai(3,kr)=-999.
                ier=-2
                return
             endif ! end of the iteration
         endif    ! end of the attraction zone
       else
          iter=0
       endif

       kr=kr+1

       if(kr > rmax) then
          write(*,*) ' not enough space for storing the ray coordinates '
          write(*,*) ' increase value rmax in the input file timefd.par second line',rmax,kr
          write(*,*) ' and re-run the model code'
          stop
       endif

       xinorm=sngl(1.d00/dsqrt(dble(d1**2)+dble(d2**2)+dble(d3**2)))   !/sqrt(2.)
 !===================================== il faut un gradient qui ne soit pas nul ... et donc on teste
       !===================================== si eikonal bien passe, cela n'arrive pas mais ... on sait jamais
       if(xinorm > 1.e+23) then
          write(*,*) ' problem doing ray tracing, all components of gradient equal to zero ', d1,d2,d3
          write(*,*) ' report the error as unexpected result from eikonal solver'
          stop
       endif
       !JEAN write(*,*) ' d1,d2,d3 ', d1*xinorm*h,d2*xinorm*h,d3*xinorm*h,xinorm,d1,d2,d3,h
       !     calcul du point suivant du rai grace au gradient (d1,d2,d3)
       rai(1,kr)=rai(1,kr-1)-h*d1*xinorm
       rai(2,kr)=rai(2,kr-1)-h*d2*xinorm
       rai(3,kr)=rai(3,kr-1)-h*d3*xinorm

       if(debug > 1) write(*,*) ' kr, rai ' ,kr,rai(1,kr),rai(2,kr),rai(3,kr),debug,h,d1,d2,d3,xinorm
       !     changement de coordonnees       
       xr=0.;yr=0.;zr=0.
       call ccoo(xg,yg,zg,n1,n2,n3,h1,h2,h3,rai(1,kr),rai(2,kr),rai(3,kr),xr,yr,zr,ier)
       if(ier.ne.0) then
          write(*,*) ' warning error in sub_rai when tracing ray: point outside the box '
          write(*,*) ' warning removing the related data '
          write(*,*) ' ier,rai_x,rai_y,rai_z ',ier,rai(1,kr),rai(2,kr),rai(3,kr)
          kr=1    ! erase computed values in order to cancel this ray
          rai(1,kr)=-999.;rai(2,kr)=-999.;rai(3,kr)=-999.
          ier=-2
          return
       endif
       !     calcul du segment de rai : pas a effectuer selon la plus grande pente
       call seg_rai(xr,yr,zr,t,n1,n2,n3,d1,d2,d3)
    enddo   ! end of the while   bounded loop by rmax
    !----------------------- decalage entre dernier point et station
    off1=xo-rai(1,kr)
    off2=yo-rai(2,kr)
    off3=zo-rai(3,kr)
    !----------------------- deformation du rai pour coincidence finale
    if(kr > ismooth) then  ! si rayon assez long
       do i=1,ismooth
          xl=float(i)/float(ismooth)
          rai(1,kr-ismooth+i)=rai(1,kr-ismooth+i)+xl*off1
          rai(2,kr-ismooth+i)=rai(2,kr-ismooth+i)+xl*off2
          rai(3,kr-ismooth+i)=rai(3,kr-ismooth+i)+xl*off3
       enddo
    else
       ! do not 
       write(*,*) ' warning ####################IN '
       write(*,*) ' less than 5 points on this ray: increase the number of points',kr
       write(*,*) ' decrease the ray sampling (model.head)'
       write(*,*) ' and re-run the code model'
       write(*,*) ' initial position ',xe,ye,ze
       write(*,*) ' final position ',xo,yo,zo
       write(*,*) ' warning ####################OUT '
    endif
    return
  end subroutine sub_rai

  !##############################################################################
  !     conversion des coordonnees dans la grille P&L
  !     
  !     input : xo,yo,zo  : origine de la grille P&L
  !             n1,n2,n3  : grille P&L
  !             h1,h2,h3  : pas P&L   h1=h2=h3
  !             xanc,yanc,zanc  : point exprime dans le repere physique

  !     output :xnew,ynew,znew  : point exprime dans le repere P&L    
  !
  subroutine ccoo(xo,yo,zo,n1,n2,n3,h1,h2,h3,xanc,yanc,zanc,xnew,ynew,znew,ier)
    implicit none
    integer(kind=4) :: n1,n2,n3,ier
    real(kind=4) :: xo,yo,zo,h1,h2,h3
    real(kind=4) :: xanc,yanc,zanc,xnew,ynew,znew
    ier=0
    !     changement de repere   ! still real values although in the grid units
    xnew=(xanc-xo)/h1+1
    ynew=(yanc-yo)/h2+1
    znew=(zanc-zo)/h3+1
    !     
    if(xnew > float(n1) .or. ynew > (float(n2)) .or. znew > (float(n3)) .or. &
         xnew <    1      .or. ynew <     1       .or. znew < 1            ) then ! le point est hors du domaine P&L
       write(*,*) ' point out of domain'
       write(*,*) ' Cartesian exact position ',xanc,yanc,zanc
       write(*,*) ' Cartesian x segment ',xo,xo+(n1-1)*h1
       write(*,*) ' Cartesian y segment ',yo,yo+(n2-1)*h2
       write(*,*) ' Cartesian z segment ',zo,zo+(n3-1)*h3
       write(*,*) ' Grid position (fractional indexes) ',xnew,ynew,znew
       write(*,*) ' Grid minimal indexes above (2,2,2): 1 1 1 ',1,1,1
       write(*,*) ' Grid maximal indexes below (n1-1,n2-1,n3-1): n1,n2,n3 ',n1,n2,n3
       write(*,*) 'WARNING: something is going wrong in subroutine ccoo'
       ier=999    ! flag something is going wrong in ccoo 
       return
    endif
    return
  end subroutine ccoo

  !###############################################################################
  !     calcul du gradient de la table des temps sur la grille P&L 
  !     au point (xr,yr,zr).
  !
  !     input : n1,n2,n3     : grille P&L 
  !             t(n1,n2,n3)  : table des temps
  !             xr,yr,zr     : point courant du rai
  !
  !
  !     output : d1,d2,d3 : gradient de t au point xr,yr,zr.
  !
  ! 
  subroutine seg_rai(xr,yr,zr,t,n1,n2,n3,d1,d2,d3)
    implicit none
    integer(kind=4) n1,n2,n3,i,j,k
    real(kind=4) t(n1,n2,n3)
    real(kind=4) xr,yr,zr
    real(kind=4) d1,d2,d3
    real(kind=4) d1x,d1y,d1z 
    real(kind=4) d2x,d2y,d2z
    real(kind=4) d3x,d3y,d3z  
    real(kind=4) d4x,d4y,d4z
    real(kind=4) d5x,d5y,d5z 
    real(kind=4) d6x,d6y,d6z  
    real(kind=4) d7x,d7y,d7z
    real(kind=4) d8x,d8y,d8z
    !
    i=int(xr)
    j=int(yr)
    k=int(zr)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ double checking here
if(i < 1 .or. i > n1 .or. j < 1 .or. j > n2 .or. k < 1 .or. k > n3) then
  write(*,*) ' internal error in subgrad for ray tracing: please report ',i,j,k
  stop
endif

    !-------------------1. Calcul du gradient aux 8 sommets de la cellule
    call subgrad(i  ,j  ,k  ,d1x,d1y,d1z,t,n1,n2,n3)
    !--
    if(i /= n1) then
       call subgrad(i+1,j  ,k  ,d2x,d2y,d2z,t,n1,n2,n3)
    else
       d2x=d1x
       d2y=d1y
       d2z=d1z
    endif
    !--
    if(i /= n1 .and. j /= n2) then
       call subgrad(i+1,j+1,k  ,d3x,d3y,d3z,t,n1,n2,n3)
    elseif (i == n1 .and. j /= n2) then
       call subgrad(i  ,j+1,k  ,d3x,d3y,d3z,t,n1,n2,n3)
    elseif (i /= n1 .and. j == n2) then
       call subgrad(i+1,j  ,k  ,d3x,d3y,d3z,t,n1,n2,n3)
    elseif (i == n1 .and. j == n2) then
       call subgrad(i  ,j  ,k  ,d3x,d3y,d3z,t,n1,n2,n3)
    else 
       write(*,*) 'warning: bug seg_rai in SUBMAIN OF RAYTRACING'
    endif
    !--         
    if(j /= n2) then
       call subgrad(i  ,j+1,k  ,d4x,d4y,d4z,t,n1,n2,n3)
    else
       call subgrad(i  ,j  ,k  ,d4x,d4y,d4z,t,n1,n2,n3)
    endif
    !--
    if(k /= n3) then
       call subgrad(i  ,j  ,k+1,d5x,d5y,d5z,t,n1,n2,n3)
    else
       call subgrad(i  ,j  ,k  ,d5x,d5y,d5z,t,n1,n2,n3)
    endif
    !--
    if(i /= n1 .and. k /= n3) then
       call subgrad(i+1,j  ,k+1,d6x,d6y,d6z,t,n1,n2,n3)
    elseif (i == n1 .and. k /= n3) then
       call subgrad(i  ,j  ,k+1,d6x,d6y,d6z,t,n1,n2,n3)
    elseif (i /= n1 .and. k == n3) then
       call subgrad(i+1,j  ,k  ,d6x,d6y,d6z,t,n1,n2,n3)
    elseif (i == 1 .and. k == n3) then 
       call subgrad(i  ,j  ,k  ,d6x,d6y,d6z,t,n1,n2,n3)
    else
       write(*,*) 'warning: bug seg_rai in SUBMAIN OF RAYTRACING'
    endif
    !--
    if(i /= n1 .and. j /= n2 .and. k /= n3) then
       call subgrad(i+1,j+1,k+1,d7x,d7y,d7z,t,n1,n2,n3) !7
    else if (i == n1 .and. j /= n2 .and. k /= n3) then
       call subgrad(i  ,j+1,k+1,d7x,d7y,d7z,t,n1,n2,n3) !8 
    else if (i == n1 .and. j == n2 .and. k /= n3) then
       call subgrad(i  ,j  ,k+1,d7x,d7y,d7z,t,n1,n2,n3) !5
    else if (i == n1 .and. j == n2 .and. k == n3) then 
       call subgrad(i  ,j  ,k  ,d7x,d7y,d7z,t,n1,n2,n3) !1
    else if (i /= n1 .and. j /= n2 .and. k == n3) then
       call subgrad(i+1,j+1,k ,d7x,d7y,d7z,t,n1,n2,n3)  !3
    else if (i /= n1 .and. j == n2 .and. k /= n3) then
       call subgrad(i+1,j  ,k+1,d7x,d7y,d7z,t,n1,n2,n3) !6
    else if (i == n1 .and. j /= n2 .and. k == n3) then
       call subgrad(i  ,j+1,k  ,d7x,d7y,d7z,t,n1,n2,n3) !4
    else if (i /= n1 .and. j == n2 .and. k == n3) then
       call subgrad(i+1,j  ,k  ,d7x,d7y,d7z,t,n1,n2,n3) !2
    else
       write(*,*) 'warning: bug seg_rai in SUBMAIN OF RAYTRACING'
    endif
    !--
    if(j /= n2 .and. k /= n3) then
       call subgrad(i  ,j+1,k+1,d8x,d8y,d8z,t,n1,n2,n3)
    else if (j /= n2 .and. k == n3) then
       call subgrad(i  ,j+1,k  ,d8x,d8y,d8z,t,n1,n2,n3)
    else if (j == n2.and.k /= n3) then
       call subgrad(i  ,j  ,k+1,d8x,d8y,d8z,t,n1,n2,n3)
    else if (j == n2.and.k == n3) then
       call subgrad(i  ,j  ,k  ,d8x,d8y,d8z,t,n1,n2,n3)
    endif
    !-------------------2. Calcul du gradient en (xr,yr,zr)
    !                      (par interpolation)
    call interpol3ds(float(i),float(i+1),float(j),float(j+1),float(k),float(k+1),d1x,d2x,d3x,d4x,d5x,d6x,d7x,d8x,xr,yr,zr,d1)
    !      
    call interpol3ds(float(i),float(i+1),float(j),float(j+1),float(k),float(k+1),d1y,d2y,d3y,d4y,d5y,d6y,d7y,d8y,xr,yr,zr,d2)
    !      
    call interpol3ds(float(i),float(i+1),float(j),float(j+1),float(k),float(k+1),d1z,d2z,d3z,d4z,d5z,d6z,d7z,d8z,xr,yr,zr,d3)
    !
    return
  end subroutine seg_rai

  !###############################################################################
  !     Calcul du gradient de t(n1,n2,n3) au noeud (i,j,k) par 
  !     differences finies.
  !
  !     input : n1,n2,n3   : grille P&L
  !             i,j,k      : noeud considere 
  !             t(n1,n2,n3): table des temps 
  !
  !     output : d1,d2,d3 : grad[t(i,j,k)]
  !
  subroutine subgrad(i,j,k,d1,d2,d3,t,n1,n2,n3)
    implicit none
    integer(kind=4) :: i,j,k,n1,n2,n3
    real(kind=4) :: d1,d2,d3,t(n1,n2,n3),t1,t2,div
    !
    if(i > 1) then 
       t1=t(i-1,j,k)
    else
       t1=t(i,j,k)
    endif
    !
    if(i < n1) then
       t2=t(i+1,j,k)
    else
       t2=t(i,j,k)
    endif
    !
    if(i == 1 .or. i == n1) then
       div=1.
    else
       div=2.
    endif
    d1=(t2-t1)/div

!!!!!!!!!!
    if(j > 1) then 
       t1=t(i,j-1,k)
    else
       t1=t(i,j,k)
    endif
    !
    if(j < n2) then
       t2=t(i,j+1,k)
    else
       t2=t(i,j,k)
    endif
    if(j == 1 .or. j == n2) then
       div=1.
    else
       div=2.
    endif
    d2=(t2-t1)/div

!!!!!!!!!!
    if(k > 1) then 
       t1=t(i,j,k-1)
    else
       t1=t(i,j,k)
    endif
    !
    if(k < n3) then
       t2=t(i,j,k+1)
    else
       t2=t(i,j,k)
    endif
    !
    if(k == 1 .or. k == n3) then
       div=1.
    else
       div=2.
    endif
    d3=(t2-t1)/div

    return
  end subroutine subgrad

END MODULE s_ray


