MODULE s_new_parameter
contains

! ###########################################################################
!----------------------------------------------------------------------------
!                 model update
!    
!
!     input : n1,n2,n3 
!             ni1,ni2,ni3
!             ixo,iyo,izo
!             s(ni1,ni2,ni3) : solution( perturbation)  lentP ou lentS ou U ou V
!             v(n1,n2,n3)    : modele vitP ou vitS ou 1./U ou 1./V
!
!     output : v(n1,n2,n3)  : nouveau modele P ou S ou 1./U ou 1./V
!
!
!     remarque : le nouveau modele est:
!
!     pour les vitesses: 
!
!            Vit=1/(lent+delta_lent)
!
!     pour les paramtres U ou V:
!
!            1./U=1/(U+delta_U)
! 
subroutine sub_new_model(v,s,n1,n2,n3,ni1,ni2,ni3,ixo,iyo,izo,d_max,iopt)
implicit none
integer(kind=4) :: i,j,k,iopt,noffset,itot
integer(kind=4) :: n1,n2,n3,ni1,ni2,ni3,ixo,iyo,izo
real(kind=4), dimension(:,:,:) :: v(n1,n2,n3)
real(kind=4), dimension(:) :: s(:)      ! size depends on option P or (P & S)
real(kind=4) :: d_max,rdum,signe

noffset=0
if(iopt /= 0) then       ! S velocity or V parameter   uvs is used
  noffset=ni1*ni2*ni3
endif


signe=0. 
do k=1,ni3               ! loop over invertable nodes
  do j=1,ni2
    do i=1,ni1
      rdum=v(i+ixo-1,j+iyo-1,k+izo-1) ! save previous value
!          indice 1D  for the solution either P or S       noffset assure l'offset dans tableau s(*)
! if iopt=0 (1,1,1) itot=ni1+1+(ni1-1)*ni2+(ni1-1)*(ni2-1)*ni3
!                   itot=ni1*ni2*ni3-ni1*ni3-ni2*ni3+ni3+ni1*ni2-ni2+ni1
!                   itot=ni1*ni2*ni3+1                    
! if iopt/=0 '1,1,1' itot=1  OK
      itot=i+ni1*(j-1)+ni1*ni2*(k-1)+noffset
      v(i+ixo-1,j+iyo-1,k+izo-1)= 1. / ((1./v(i+ixo-1,j+iyo-1,k+izo-1) )+s(itot))  ! slowness !
      if(abs(v(i+ixo-1,j+iyo-1,k+izo-1)-rdum) > d_max) then
        signe=sign(1.,v(i+ixo-1,j+iyo-1,k+izo-1)-rdum)
        v(i+ixo-1,j+iyo-1,k+izo-1)=rdum+signe*d_max
      endif
!      write(56,*) i,j,k,rdum,s(itot),itot,v(i+ixo-1,j+iyo-1,k+izo-1)
    enddo
  enddo
enddo
return
end subroutine sub_new_model

! #########################################################################
!----------------------------------------------------------------------------
!
!                     fsrc update
!
!   sources parameters updates : position, origin time 
!   and index (for next derivative)
!   Rewrite fsrc with updates.
!
!
subroutine sub_new_fsrc(p_x,p_y,p_z,p_to,ki_pos,ki_to,ki_m,id_src,sol, &
                        nsrc,usrc,neqks,nshot,nblast,neqks_out,     &
                        xo,yo,zo,n1,n2,n3,h1,h2,h3,ni1,ni2,ni3, &
                        dhori_max,dvert_max,dto_max,iopt,flog)
!.........................................................................
!
!   nsrc        : sources number
!   neqks       : earthquakes number
!   usrc        : source file
!   xo,yo,zo    : model grid origin 
!   h1,h2,h3    : model grid steps
!   n1,n2,n3    : model grid nodes number 
!   p_x(nsrc)   : source position 
!   p_y
!   p_z
!   p_to        : source origin time
!   sol(nc)     : position perturbations if offset by nvit=n1*n2*n3 (only P) or 2*n1*n2*n3 (P&S)
!               : (x,y,z) multiplexed by 3*(ki_pos(isrc)-1)
!   sol(nc)     : origin time perturbations if offset by nvit+3*neqks
!   id_src(nsrc): id de la source (defini pour toute la manip)
!   ki_pos(nsrc): used for multiplexed solutions (only incremented for quakes inside the domain
!   ki_to(nsrc) : used for indexed solutions
!   ki_m(nsrc)  : source index
!                  ki_m(nsrc)  = 0 => source outside grid
!                  ki_m(nsrc)  = 1 => source inside grid
!                  ki_to(nsrc) <> 0 => origin time parameter is inversible (increment)
!                  ki_pos(nsrc) <> 0 => position parameters are inversible  (increment)
!   iopt        : 0   for P  and /=0 for P & S
!   flog        :   output
!.........................................................................
implicit none

integer(kind=4), dimension(:)  :: ki_pos(nsrc),ki_to(nsrc),ki_m(nsrc)
integer(kind=4), dimension(:) :: id_src(nsrc)

integer(kind=4) :: noffset_v,noffset_pos,iout,iopt,flog
real(kind=4) :: xo,yo,zo,h1,h2,h3,xp,yp,zp
integer(kind=4) :: n1,n2,n3,ni1,ni2,ni3

real(kind=4), dimension(:) :: p_x(:),p_y(:),p_z(:)
real(kind=4), dimension(:) :: p_to(:)
real(kind=4), dimension(:) :: sol(:)    ! alias dp(3*eqks), dt(neqks+nblasts)

real(kind=4) :: dp_x,dp_y,dp_z,dt

integer(kind=4) :: nsrc,usrc,irec
integer(kind=4) :: kpos,kto,neqks,neqks_out,nshot,nblast,isrc
integer(kind=4) :: neqks_tot,nshot_tot,nblast_tot
real(kind=4) :: dto_max,dhori_max,dvert_max,signe 

integer(kind=4) :: nl

! --------------------- offset for velocity parameters
noffset_v=ni1*ni2*ni3     ! P velocity field
if(iopt /= 0) then        ! iopt alias uvs   = 0 for only P
  noffset_v=noffset_v+noffset_v      ! P & S
endif
! --------------------- offset for position parameters  (quakes only)
noffset_pos=noffset_v+3*neqks ! number of eqks

!     update new sources parameters
do isrc=1,nsrc ! loop over sources

  if(ki_pos(isrc) /= 0) then ! need to update position
    dp_x=sol(noffset_v+3*(ki_pos(isrc)-1)+1)
    dp_y=sol(noffset_v+3*(ki_pos(isrc)-1)+2)
    dp_z=sol(noffset_v+3*(ki_pos(isrc)-1)+3)
!     update position 
    signe=0.
!    write(*,*) isrc,dp_x,dp_y,dp_z
    IF(abs(dp_x) <= dhori_max ) then
      p_x(isrc)=p_x(isrc)+dp_x
    ELSE
      signe=sign(1.,dp_x)
      p_x(isrc)=p_x(isrc) + signe*dhori_max
    ENDIF
    IF(abs(dp_y) <= dhori_max ) then
      p_y(isrc)=p_y(isrc)+dp_y
    ELSE
      signe=sign(1.,dp_y)
      p_y(isrc)=p_y(isrc) + signe*dhori_max
    ENDIF
    IF(abs(dp_z) <= dvert_max) then
      p_z(isrc)=p_z(isrc)+dp_z
    ELSE
      signe=sign(1.,dp_z)
      p_z(isrc)=p_z(isrc) + signe*dhori_max
    ENDIF
  endif

  if(ki_to(isrc) /= 0) then ! need to update origin time (quakes and blasts)
    dt=sol(noffset_pos+ki_to(isrc))
    IF(abs(dt) <= dto_max) then
      p_to(isrc)=p_to(isrc)+dt
    ELSE
      signe=sign(1.,dt)
      p_to(isrc)=p_to(isrc) + signe*dto_max
    ENDIF
  endif
!     check if updated position is inside or outside model grid
!     ki_m(isrc)=0 => source outside domain
!     ki_m(isrc)=1 => source inside domain
!     call in_out
  xp=p_x(isrc)           ! new position
  yp=p_y(isrc)
  zp=p_z(isrc)
  call in_out(xo,yo,zo,h1,h2,h3,n1,n2,n3,xp,yp,zp,iout)
  ki_m(isrc)=iout
  IF(ki_m(isrc) == 0) write(*,*) 'WARNING LOST ONE SOURCE',ki_m(isrc) 
end do

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     analyzing quakes, shots and blasts
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

neqks_tot=neqks+neqks_out     ! normalement reste le meme
nshot_tot=nshot
nblast_tot=nblast

kpos=0
kto=0
neqks=0
nshot=0
nblast=0
neqks_out=0
irec=0

do isrc=1,neqks_tot                  ! loop over sources for quakes
  if(ki_pos(isrc) /= 0 .and. ki_to(isrc) /= 0) then
!     earthquakes
    if(ki_m(isrc) == 1) then ! inside
      kpos=kpos+1            ! increment for further inversion  
      kto=kto+1
      neqks=neqks+1
      ki_pos(isrc)=kpos      ! put index for further inversion (same as neqks)
      ki_to(isrc)=kto        ! put index for further inversion (same as neqks)
      irec=irec+1
      write(usrc,rec=irec+1) p_x(isrc),p_y(isrc),p_z(isrc),p_to(isrc),&
         ki_pos(isrc),ki_to(isrc),ki_m(isrc),id_src(isrc)
    else                     !outside
      neqks_out=neqks_out+1
      ki_pos(isrc)=0
      ki_to(isrc)=0  
      write(usrc,rec=irec+1) p_x(isrc),p_y(isrc),p_z(isrc),p_to(isrc),&
              ki_pos(isrc),ki_to(isrc),ki_m(isrc),id_src(isrc)
    endif
  endif
enddo

write(*,*) ' total declared earthquakes ',neqks_tot
write(*,*) ' found earthquakes ',neqks+neqks_out
write(*,*) ' earthquakes inside the box ',neqks
write(*,*) ' earthquakes outside the box ',neqks_out

write(flog,*) ' total declared earthquakes ',neqks_tot
write(flog,*) ' found earthquakes ',neqks+neqks_out
write(flog,*) ' earthquakes inside the box ',neqks
write(flog,*) ' earthquakes outside the box ',neqks_out

do isrc=neqks_tot+1,neqks_tot+nshot_tot                 ! loop over sources for shots
!     shot
  if(ki_m(isrc) == 0) then
    write(*,*) ' shot ',id_src(isrc),' is declared outside ?: should not happen'
  endif
  if(ki_pos(isrc) /= 0 .or. ki_to(isrc) /= 0) then
    write(*,*) ' should not be observed ',ki_pos(isrc),ki_to(isrc),' for shot ',id_src(isrc)
    ki_pos(isrc)=0         ! erase them
    ki_to(isrc)=0
  endif
  nshot=nshot+1
  irec=irec+1              ! write it
  write(usrc,rec=irec+1) p_x(isrc),p_y(isrc),p_z(isrc),p_to(isrc), &
       ki_pos(isrc),ki_to(isrc),ki_m(isrc),id_src(isrc)
enddo

write(*,*) ' total declared shots ',nshot_tot
write(*,*) ' found shots ',nshot

write(flog,*) ' total declared shots ',nshot_tot
write(flog,*) ' found shots ',nshot

do isrc=neqks_tot+nshot+1,neqks_tot+nshot_tot+nblast_tot               ! loop over sources for blasts
!     blast
  if(ki_m(isrc) /= 0) then
    write(*,*) ' blast ',id_src(isrc),' is declared outside ?: should not happen'
  endif
  if(ki_pos(isrc) /= 0) then
    write(*,*) ' blast ',id_src(isrc),' has been declared having a moving position '
    ki_pos(isrc)=0    ! on corrige
  endif
  if(ki_to(isrc) == 0) then
    write(*,*) ' blast ',id_src(isrc),' has been declared having a fixed origin time '
  endif    ! on corrige apres
  nblast=nblast+1
  kto=kto+1
  ki_to(isrc)=kto       ! put index for further inversion (same as nblast (+neqks) )
  irec=irec+1
  write(usrc,rec=irec+1) p_x(isrc),p_y(isrc),p_z(isrc),p_to(isrc),   &
            ki_pos(isrc),ki_to(isrc),ki_m(isrc),id_src(isrc)
enddo

write(*,*) ' total declared blasts ',nblast_tot
write(*,*) ' found shots ',nblast

write(flog,*) ' total declared blasts ',nblast_tot
write(flog,*) ' found shots ',nblast

!
!  verification si toutes les sources
!

if(nsrc == neqks+nblast+nshot+neqks_out) then
  write(usrc,rec=1) neqks+nblast+nshot+neqks_out,neqks,nshot,nblast,neqks_out,0,0,0
else
  write(*,*) ' ERROR in the sub_new_fsrc subroutine from new_parameters '
  write(*,*) ' we have lost events ! '
  stop
endif
return
end subroutine sub_new_fsrc

! ##############################################################
!------check if source is inside domain
!     input :
!     xo,yo,zo : domain origin
!     n1,n2,n3 : nodes numbers
!     h1,h2,h3 : steps
!     p : source position
!
!     output :
!     in_out = 1 if inside domain
!     in_out = 0 if outside domain
!
subroutine in_out(xo,yo,zo,h1,h2,h3,n1,n2,n3,xp,yp,zp,iout)
implicit none
integer(kind=4) :: n1,n2,n3,iout
real(kind=4) :: xo,yo,zo,h1,h2,h3,xp,yp,zp
 
if( xp > xo .and. yp > yo .and. zp > zo .and. &
    xp < float(n1-1)*h1+xo .and. yp < float(n2-1)*h2+yo .and. zp < float(n3-1)*h3+zo) then 
  iout=1
else
  iout=0
end if
return
end subroutine in_out

END MODULE s_new_parameter
