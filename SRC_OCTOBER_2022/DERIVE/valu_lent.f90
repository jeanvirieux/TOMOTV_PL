MODULE s_valu_lent
contains

! ############################################################################
!
!
!     READ SLOWNESS
!
!     Determination de la lenteur aux hypocentres.
!
!     input:
!     flog     : unite logique fichier log
!     xo,yo,zo : origine modele
!     n1,n2,n3 : grille modele
!     h1,h2,h3 : 
!     nsrc     : nombre de sources
!     uvp,fvp  : P velocity unit and file
!     uvs,fvs  : S velocity unit and file (used as a flag uvs=0)
!     
!     output :
!     plent(nsrc) : lenteurs P aux foyers     (allocated HERE)
!     slent(nsrc) : lenteurs S aux foyers     (allocated HERE)
!
! ###########################################################################
subroutine valu_lent(uvp,uvs,fvp,fvs,sr_x,sr_y,sr_z,ki_pos,nsrc,   &
                     xo,yo,zo,n1,n2,n3,h1,h2,h3, &
                     plent,slent,flog)
implicit none
integer(kind=4) :: flog,uvp,uvs
integer(kind=4) :: isrc,nsrc
integer(kind=4) :: n1,n2,n3
real(kind=4) :: xo,yo,zo,h1,h2,h3
integer(kind=4) :: kp
real(kind=4) :: xs,ys,zs
character(len=*) :: fvp,fvs

real(kind=4), dimension(:) :: sr_x(nsrc),sr_y(nsrc),sr_z(nsrc)
integer(kind=4), dimension(:) :: ki_pos(nsrc)

real(kind=4), dimension(:) :: plent,slent
real(kind=4), allocatable,dimension(:,:,:) :: vit

write(flog,'(10x,a)') ' setting hypocenters P slowness'

plent(:)=0.0

allocate(vit(n1,n2,n3))

! should read the velocity file fvp
open(uvp,file=fvp,access='direct',recl=4*n1*n2*n3)
read(uvp,rec=1) vit
close(uvp)

do isrc=1,nsrc
! position de la source
  xs=sr_x(isrc)
  ys=sr_y(isrc)
  zs=sr_z(isrc)
  kp=ki_pos(isrc)
  if(kp /= 0) then      ! quakes
    call reads(xs,ys,zs,xo,yo,zo,n1,n2,n3,h1,h2,h3,plent,vit,nsrc,isrc)            ! avec n*
  endif
enddo

if(uvs.ne.0) then
  write(flog,'(10x,a)') ' setting hypocenters S slowness'
  slent(:)=0.0
                         ! should read the velocity file fvs
  open(uvs,file=fvs,access='direct',recl=4*n1*n2*n3)
  read(uvs,rec=1) vit
  close(uvs)
  do isrc=1,nsrc
    xs=sr_x(isrc)
    ys=sr_y(isrc)
    zs=sr_z(isrc)
    kp=ki_pos(isrc)
    if(kp /= 0) then    ! quakes
      call reads(xs,ys,zs,xo,yo,zo,n1,n2,n3,h1,h2,h3,slent,vit,nsrc,isrc)        ! A COMPRENDRE
    endif
  enddo
endif
!============= free memory   (vit but not plent and slent)
deallocate(vit)
!     
return
end subroutine valu_lent

! ##############################################################################
!...............................................................................
!     Evaluation de la lenteur au foyer par interpolation. 
!
!     input :
!     x,y,z  : coordonnees de la source
!     xo,yo,zo : origine de la grille du modele
!     n1,n2,n3 : grille du modele
!     h1,h2,h3 : pas
!     src : nombre total de sources
!     isrc : numero de la source courante
!     v(n1,n2,n3) : vitesses 
!
!     output :
!     slow(isrc) : lenteur a la source isrc
! ##############################################################################

subroutine reads(x,y,z,xo,yo,zo,n1,n2,n3,h1,h2,h3,slow,v,nsrc,isrc)
use s_interpol3ds

implicit none
integer(kind=4) :: n1,n2,n3,nsrc,isrc
real(kind=4) :: x,y,z,xs,ys,zs,xo,yo,zo
real(kind=4), dimension(:) :: slow(nsrc)
real(kind=4), dimension(:,:,:) :: v(n1,n2,n3)
real(kind=4) :: s1,s2,s3,s4,s5,s6,s7,s8,slow_s
integer(kind=4) :: ix,iy,iz
real(kind=4) :: h1,h2,h3
!
      ix=int((x-xo)/h1)+1
      iy=int((y-yo)/h2)+1
      iz=int((z-zo)/h3)+1
! ============================== slowness from velocity : "slowness" interpolation OK
      s1=1./v(ix  ,iy  ,iz  )
      s2=1./v(ix+1,iy  ,iz  )
      s3=1./v(ix+1,iy+1,iz  )
      s4=1./v(ix  ,iy+1,iz  )
      s5=1./v(ix  ,iy  ,iz+1)
      s6=1./v(ix+1,iy  ,iz+1)
      s7=1./v(ix+1,iy+1,iz+1)
      s8=1./v(ix  ,iy+1,iz+1)
!
      xs=(float(ix-1)*h1)+xo
      ys=(float(iy-1)*h2)+yo
      zs=(float(iz-1)*h3)+zo
!
call interpol3ds(xs,xs+h1,ys,ys+h2,zs,zs+h3,s1,s2,s3,s4,s5,s6,s7,s8,x,y,z,slow_s)
slow(isrc)=slow_s
return
end subroutine reads

END MODULE s_valu_lent

