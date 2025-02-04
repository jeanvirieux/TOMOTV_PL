!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! from slowness vector, provide azimuth wrt North (py) and dip wrt horizon
! lower-sphere projection (pz >0 downward)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine slowness2angle(px,py,pz,azimuth,dip)
  implicit none

  real(kind=4) :: px,py,pz,azimuth,dip
  real(kind=4) :: pi=3.141592653

  real(kind=8) :: px1,py1,pz1,pnorm

  !
  ! make sure that we live in the lower hemisphere
  !
  if(pz <0.0) then
     px1=-px;py1=-py;pz1=-pz
  else
     px1=px;py1=py;pz1=pz
  endif
  !
  ! now we leave in a safe domain
  !
  if(pz1 > 0.0) then
     pnorm=sqrt(px1**2+py1**2+pz1**2)
     px1=px1/pnorm;py1=py1/pnorm;pz1=pz1/pnorm
     azimuth=sngl(atan(px1/py1))*180./pi             ! droite haut par defaut
     !JEAN     write(*,*) px1,py1,pz1,pnorm
     if(py1 < 0. .and. px1 > 0.) then          ! gauche haut
        azimuth=180.00+azimuth
     elseif(py1 > 0. .and. px1 < 0.) then      ! droite bas
        azimuth=360.+azimuth
     elseif(py1 < 0. .and. px1 < 0.) then      ! gauche bas
        azimuth=180.00+azimuth
     endif
     dip=sngl(asin(pz1))*180./pi
     !JEAN    write(*,*) px1,py1,pz1,dip,azimuth
  else
     ! pz=0 dip=0°   ! horizontal angle
     dip=0.00
     pnorm=sqrt(px1**2+py1**2)
     px1=px1/pnorm;py1=py1/pnorm
     azimuth=sngl(atan(px1/py1))*180.00/pi  ! between -pi/2 to pi/2
     if(py1 < 0. .and. px1 > 0.) then          ! gauche haut
        azimuth=180.00+azimuth
     elseif(py1 > 0. .and. px1 < 0.) then      ! droite bas
        azimuth=360.+azimuth    !;write(*,*) ' az ',azimuth
     elseif(py1 < 0. .and. px1 < 0.) then      ! gauche bas
        azimuth=180.00+azimuth
     endif
  endif
  return
end subroutine slowness2angle
