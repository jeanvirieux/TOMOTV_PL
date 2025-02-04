MODULE s_read_input
contains

!################################################################
!
!
!################################################################
  subroutine read_par(vip,vis,se,to,ni1,ni2,ni3,flog)
    implicit none
    integer(kind=4) :: vip,vis,se,to
    integer(kind=4) i,ni1,ni2,ni3,ixo,iyo,izo,n1,n2,n3,uvp,uvs
    real(kind=4) :: xo,yo,zo,h1,h2,h3
    integer(kind=4) :: src,eqks,blast,shot,eqks_out,idum1,idum2,idum3,flog
    integer(kind=4) :: chx
    character(len=132) name
    i=0
    open(49,file='inversion.par',status='old',err=1000)
    read(49,*) chx
    read(49,*) xo,yo,zo
    read(49,*) n1,n2,n3
    read(49,*) h1,h2,h3
    read(49,'(a)') name
    uvp=i
    i=i+1
    read(49,*) ixo,iyo,izo
    read(49,*) ni1,ni2,ni3
    vip=ni1*ni2*ni3
    !============================ if we have also Vs ... write the name of the file
    if(chx == 2) then
       uvs=i
       i=i+1
       read(49,'(a)') name
       vis=vip
    else
       uvs=0
       vis=0
    endif
    !
    close(49)

    open(49,file='fsrc',access='direct',recl=4*8)
    read(49,rec=1) src,eqks,shot,blast,eqks_out,idum1,idum2,idum3
    close(49)
    !
    se=(3*eqks)      !STEPHANIE DD    CONSERVE   JEAN
    to=(eqks+blast)  !STEPHANIE DD    CONSERVE   JEAN   MAIS PAS DD
    !
    return
1000 continue
    write(flog,*) ' error in reading the file inversion.par '
    write(*,*) ' error in reading the file inversion.par '
    stop
  end subroutine read_par
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
!     lecture des parametres l0,l1,r0,r1 pour calcul D1
!     mais aussi des echelles arbitraires entre classes de parametres
!     cp sur Up ou Vp
!     cs sur Us ou Vs
!     cpo sur xo,yo,zo
!     cto sur to
!     stat_res ...
!
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
subroutine read_ent(cp,cs,cpo,cto,lx,ly,lz,r0,r1,l0,l1,     &
                           stat_res,lap,lap1,irep,flog) 
implicit none
real(kind=4) :: cp,cs,cpo,cto,lx,ly,lz,l1,l0,r0,r1
integer(kind=4) :: stat_res,flog,lap,lap1,irep
character(len=1) :: carac

open(48,file='inversion.head',status='old',err=100)
read(48,'(a)') carac
read(48,*) r0,r1
write(flog,*) ' residues (sec) in picked time for weighting starting at ', r0,' and ', r1, 'zero after that'
read(48,'(a)') carac
read(48,*) l0,l1
write(flog,*) ' residues (meters) in ray length for weighting starting at ', l0,' and ', l1, 'zero after that'
read(48,'(a)') carac
read(48,*) lap     ! flag=1 for smoothing penalty 
write(flog,*) ' smoothing penalty option (1=yes) ',lap
read(48,'(a)') carac
read(48,*) lap1     ! option from 1 to 3 for smoothing 
write(flog,*) ' option for total laplacian, horiz/verti laplacian or Cartesian laplacian '
read(48,'(a)') carac
read(48,*) lx      ! coefficient of smoothing along x
read(48,'(a)') carac
read(48,*) ly      ! coefficient of smoothing along y
read(48,'(a)') carac
read(48,*) lz      ! coefficient of smoothing along z
write(flog,*) ' smoothing x,y,z coefficients ',lx,ly,lz
read(48,'(a)') carac
read(48,*) irep,stat_res    ! option for residues
write(flog,*) ' option for setting residues weights -0 set weight to 1 and -1 read fwei:',irep
write(flog,*) ' setting statistics option (not yet) ',stat_res
read(48,'(a)') carac
read(48,*) cp,cs
read(48,'(a)') carac
read(48,*) cpo,cto
write(flog,*) ' multiparameter balancing Vp, Vs, (xquake,yquake,zquake), t_orig '
write(*,*) cp,cs,cpo,cto
close(48)
write(*,*) ' end of the reading '
return
100 continue
write(flog,*) ' missing file inversion.head '
write(*,*) ' missing file inversion.head '
stop 
end subroutine read_ent
!--

END MODULE s_read_input
