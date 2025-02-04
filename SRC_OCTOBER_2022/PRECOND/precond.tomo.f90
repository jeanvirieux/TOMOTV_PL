!################################################################################
!  Correction of a bug as we need to increase the dimensions of the system      |
!                 when considering laplacian smoothing      2014   JEAN         |
!                                                                               |
!================================================================================
!                                                                               |
!                                                                               |
!        preconditioning of the least-square system  Ax=b                       |
!                                                                               |
!                                                                               |
!        Ax=b is replaced by the system (D1.A.D2.y=D1.b,x=D2.y) for             | 
!        many reasons and one of them is the numerical stability.               |
!        Smoothing option is possible through a matrix L                        |
!        We solve :                                                      |
!                                                                               |
!                        |D1.A.D2| y - |D1.b| , x=D2.y                          |
!                        |   L   |   - | 0  |                                   |
!                                                                               |
!        On top of that, one may consider adding the damping used for example by|
!        an iterative solver as LSQR.                                           |
!                                                                               |
! ##############################################################################|
!        !!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!                     |
!  in this module, files fmat.* are erased <<<< Be careful                      |
!          ====> and save these files if you plan to reuse them                 |
!        !!!!!!!!!!!!!!!!! END OF WARNING   !!!!!!!!!!!!!!!!!!                  |  
! ##############################################################################|
!  input -----------------------------------------------------------------------|
!        r0 residue limit: weighting=1 below this value                         |
!        r1 residue limit: weighting=0 beyond this value                        |  
!        cP coefficient for P wave                                              | 
!        cS coefficient for S wave                                              |
!        cpo coefficient for X,Y,Z positions                                    |
!        cto coefficient for origin time                                        |
!        l0 ray length limit: weighting=1 below this value                      |
!        l1 ray length limit: weighting=9 beyond this value                     |
!        poids file for external weighting                                      |
!                                                                               |
!--input (stdin)----------------------------------------------------------------| 
!                                                                               |
! cp cs cpo cto lx ly lz (real)parametres of scaling (use in the D2 matrix)     | 
!                         (real)weighting of smoothing along the x,y,z direction|
!                                                                               | 
! lap           (integer) =0 no Laplacian smoothing,                            | 
!                          =1 should use Laplacian smoothing                    |
!                                                                               |
!    if lap=1: rep (integer) (smoothing definition through L)                   |
!                                                                               |
!            rep=1 lx*laplacien (x) + ly*laplacien(y) + lz*laplacien(z) = 0     |
!            rep=2 lx*laplacien (x) + ly*laplacien(y)=0,lz*laplacien(z) = 0     |
!            rep=3 lx*laplacien (x)=0,ly*laplacien(y)=0,lz*laplacien(z) = 0     |
!                                                                               |
!                                                                               |
! stat_res      (integer) =0, 1, ou 2                (for computing D1)         | 
!                                                                               |
!   rep                                              (for computing D1)         |
!                                                                               |
!     if stat_res=1 : C       (real) (between 0. et 10.)                        |
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@|
!                     coeff   (real)                                            |
!                     adviced values :                                          |
!                               coeff=1e3 for arrival times                     | 
!                               coeff=1e4 for delays                            |
!                                                                               |
!                                                                               | 
!----- definition of the D1 matrix  --------------------------------------------|
!                                                                               |
!       D1 matrix can be red from the file poids or from the file Rpoids or     |
!       in both. The total weighting will be the product of these two numbers   |
!                                                                               |
!       values of stat_res:                                                     |
!                                                                               | 
!          stat_res = 0 : nothing                                               | 
!          stat_res = 1 : statistics over residues and file Rweight is created  |
!          stat_res = 2 : reading the Rweight file                              |
!                                                                               |
!  D1      rep =  1   :   read file fwei (binaire,real*4,rec=4*nt) nt all data |
!                         as defined by the user                                |
!          rep =  0   :   disregard the file fwei                              |
!                                                                               |
!                      --------------------------------------                   |
!                      | stat_res=0 |stat_res=1 |stat_res=2 |                   |
!          |-----------+------------+-----------+-----------+                   |
!          |rep=1      |     1      |    2      |   3       |                   |
!          |-----------+------------+-----------+-----------+                   |
!          |rep=0      |     4      |    5      |   6       |                   |
!          |-----------+-------------------------------------                   |
!                                                                               |
!         1:  : D1 red in a file named fwei (binary format,real*4,rec=4*nt)    |
!               this file has to be built by the used                           |
!               D1=weight                                                       |
!                                                                               |
!         2.  : Reading the file fwei and statistics over residues leading to  |
!               another weight Rweight, giving us the                           |
!               D1=weight*Rweight                                               |
!                                                                               |
!         3.  : Reading the file fwei and reading the file Rweight             |
!               D1=weight*Rweight                                               |
!                                                                               |
!         4.  : D1=Identity                                                     |
!                                                                               |
!         5.  : Statistics over residues leading to a weight file Rweight and   |
!               D1=Rweight                                                      |
!                                                                               | 
!         6.  : Reading the file named Rweight                                  |
!               D1=Rweight                                                      |
!                                                                               |
!   on top, there is the weighting related to residues and ray lengths          |
!                                                                               |
!                                                                               |
!------- definition of D2 matrix -----------------------------------------------| 
!                                                                               |
!                                                                               |
!         1. computing the L2 norm of the columns of the sensitivity matrix A   |
!         2. computing max value of these norms for each kind of parameters     |
!            mp for P wave velocities,                                          |
!            ms for S wave velocities,                                          |
!            mh for hypocenters,                                                |
!            mt for origin times,                                               |
!         3. evaluate mp=cp/mp, ms=cs/ms , mh=cpo/mh , mt=cto/mt                |
!         4. compute D2=diag(mp...,ms...,mh...,mt...)                           |
!                                                                               |
!         adviced values : cp=0.1, cs=0.1, cpo=1., cto=0.5.                     | 
!                                                                               |
!------ definition de L---------------------------------------------------------|
!                                                                               |
!       For smoothing the tomographic image, we add complentary equations acting|
!       as contraints on the least-square system to be solved. We define three  |
!       smoothing possibilities                                                 | 

!      1. for each node, we compute  : lx*lap(x)+ly*lap(y)+lz*lap(z)=0          |
!      2. for each node, we compute  : lx*lap(x)+ly*lap(y)=0, lz*lap(z)=0       |
!      3. for each node, we compute  : lx*lap(x)=0, ly*lap(y)=0, lz*lap(z)=0    |
!                                                                               |
!     (notation : lap(u) = laplacien(u) : second derivatives in the direction u |
!------------------------------------------------------------------------------ |
!  modification we read fwei and not the poids file                            |
!################################################################################ 
program precond_smooth
use s_read_input
use s_read_obs
use s_precond
use s_lap
use s_read_matrix
use s_write_matrix
use s_write_rms
implicit none

integer(kind=4) :: nnz,nnz1,nl,nl1,nc,iunit,nnz_off,nl_off,nnz_tot,nl_tot
integer(kind=4) :: nix,niy,niz
integer(kind=4) :: vip,vis,sei,to
integer(kind=4) :: umat,uic,uid,ures,lumat,luic,luid
integer(kind=4) :: urp,uhist,flog,uw
real(kind=4)    :: cp,cs,cpo,cto,lx,ly,lz,l0,l1,r0,r1,coeff
integer(kind=4) :: stat_res,lap,lap1,irep
integer(kind=4) :: nlp,i,nt,ntp,nts,i1,i2,i3,idata,irhs

integer(kind=4), allocatable, dimension(:) :: ic,id,lic,lid
real(kind=4), allocatable, dimension(:) :: mat,lmat

real(kind=4), allocatable, dimension(:) :: d1,d2

real(kind=4), allocatable, dimension(:) :: temps,dtemps
integer(kind=4), allocatable, dimension(:) :: id_dat,lut_src,lut_sta,lut_ray
integer(kind=4), allocatable, dimension(:) :: idnl
real(kind=4), allocatable, dimension(:) :: res,lres
real(kind=4), allocatable, dimension(:) :: poids,xl_ray,weight

character(len=132) :: fmat,fic,fid,lfmat,lfic,lfid

logical :: ondesS
!
ondesS=.false.
nnz1=0
nl1=0
iunit=10
flog=iunit
iunit=iunit+1
uhist=iunit
iunit=iunit+1
urp=iunit
iunit=iunit+1
uw=iunit
iunit=iunit+1

open(flog,file='flog.precond')
!
!--------------lecture des parametres-------------------------------
!              on utilise le stdin 
call read_par(vip,vis,sei,to,nix,niy,niz,flog)
call read_ent(cp,cs,cpo,cto,lx,ly,lz,r0,r1,l0,l1,stat_res,lap,lap1,irep,flog)

!=========================================
!   read fobs for tracking used data       we need nt and id_dat
!=========================================
open(iunit,file='fobs',access='direct',recl=6*4)
read(iunit,rec=1) nt,ntp,nts,i1,i2,i3   ! read the total number of data (split into P & S data)

allocate(id_dat(nt))    ! id of the data
allocate(temps(nt))      ! both P and S waves   ntp first and then nts
allocate(dtemps(nt))
allocate(lut_src(nt))    ! id of the source
allocate(lut_sta(nt))    ! id of the station
allocate(lut_ray(nt))    ! lut for each data (=0 if no ray or iray if one ray)

call read_fobs(iunit,id_dat,temps,dtemps,lut_src,lut_sta,lut_ray,nt,ntp,nts)
close(iunit)
deallocate(temps);deallocate(dtemps);deallocate(lut_src);deallocate(lut_sta);deallocate(lut_ray)
!======== we keep id_dat because we need it for association through the idnl (fidnl)

if (vis /= 0) ondesS=.true.

!
!     matrice de lissage
!     0. :pas de matrice de lissage
!     1. :calcul d'une matrice de lissage
!
  if(lap == 1) then
     write(flog,*) ' smoothing penalty '
     !     1. laplacien=0
     !     2. laplacien horizontal=0 et laplacien vertical=0
     !     3. laplacien x=0, laplacien y=0 et laplacien z=0
     if(lap1 == 1) then
        write(flog,*) 'constraint on the laplacian =0'
        nnz1=7*(vip+vis)
        nl1=vip+vis
     elseif (lap1 == 2) then
        write(flog,*) 'constraint on the horizontal laplacian on one side =0'
        write(flog,*) ' and on the vertical laplacian on the other side =0 '
        nnz1=8*(vip+vis)
        nl1=2*(vip+vis)
     elseif (lap1 == 3) then
        write(flog,*) ' constraint on each second derivative along x,y and z =0 '
        nnz1=9*(vip+vis)
        nl1=3*(vip+vis)
     else
        write(flog,*) ' if smoothing, only options 1,2,3 are allowed '
        write(*,*) ' unknown option should be (1,2,3)',lap1
     endif
  else
     write(flog,*) ' no smoothing operation '
  endif
!
!-------------lecture du systeme lineaire a preconditionner-------------------
!
!     lecture du la matrice de la taille normale

!
! reading information for the inversion
!

write(flog,*) '---------------------------------------------------------------'
write(flog,'(35x,a,a15)') ' reading matrix.par  '
open(49,file='matrix.par')
read(49,'(a)') fmat
read(49,'(a)') fic
read(49,'(a)') fid
read(49,*) nnz
read(49,*) nl,nlp                    ! nlp ne sert a rien ici
read(49,*) nc                        ! column
close(49)

write(*,*) ' nbre of rows nl', nl
write(*,*) ' nbre of columns nc',nc
write(*,*) ' nbre of non-zero elements nnz',nnz

allocate(mat(nnz))
allocate(ic(nnz))
allocate(id(nl+1))

! logical units

umat=iunit
iunit=iunit+1
uic=iunit
iunit=iunit+1
uid=iunit
iunit=iunit+1

call read_matrix(mat,ic,id,nnz,nl,nc,umat,uic,uid,fmat,fic,fid,flog)
!====================================================
!     reading the RHS vector (residus only nl)
!====================================================
ures=iunit
iunit=iunit+1
allocate(res(nl))  ! but we may have nl+nl1 when considering laplacian smoothing
open(ures,file='fdift',access='direct',recl=4*nl)
read(ures,rec=1) res
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    ures
! ATTENTION on ne ferme pas cette unite pour l'instant
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    ures
!
!   read fidnl (dimension nl) which provides the id_dat=idnl(*) of the related data
!
allocate(idnl(nl))
open(iunit,file='fidnl',access='direct',recl=4*nl)
read(iunit,rec=1) idnl
close(iunit)
!                            ==============================================================
! read user-defined weight   if(idnl(irhs)=id_dat(idata)) then weight(idata) <--> res(irhs)
!                            ==============================================================
allocate(weight(nt))
allocate(poids(nl))        ! we must define the weight for rhs term
!  ponderation des equations: poids des fichiers
write(flog,*) '0 no reading of user weight'
write(flog,*) '1 reading of a user weight '
!write(flog,*) '2 automatic weight design through histrogram NOT IMPLEMENTED'
poids(:)=0.
!
!  option for weight management
!
if(irep == 1) then
  open(iunit,file='fwei',access='direct',recl=4*nt)
     read(iunit,rec=1) weight
     close(iunit)
!@@@@@@@@@@@@@@@@@@@@@@@@@ association between data and rhs
  do irhs=1,nl   ! loop over rhs vector for association 
    do idata=1,nt
      if(id_dat(idata) == idnl(irhs)) then
        poids(irhs)=weight(idata)
        goto 2000
      endif
    enddo   ! end of loop over data idata
    write(*,*) ' error in PRECOND: unknown id observation ',id_dat(idata),idnl(irhs)
    stop
2000 continue
  enddo
elseif(irep == 0) then
  write(*,*) ' we assume weight equal to 1: we disregard file fwei'
  poids(:)=1.
else    ! other cases  ... works for values different from 0 and 1
  write(flog,*) ' we assume that there is not user-defined weight '
  poids(:)=1.
endif
!
! make up the residues with the user-weight into the file fresv
!
!============================= we want to save weighted residues for all data
do idata=1,nt    ! loop over data
  weight(idata)=-999.
  do irhs=1,nl   ! loop over rhs vector for association 
    if(id_dat(idata) == idnl(irhs)) then
      weight(idata)=res(irhs)*poids(irhs)
!      write(77,*) ' idata, irhs ',idata,irhs,id_dat(idata),idnl(irhs)
      goto 3000
    endif
  enddo   ! end of loop over rhs
3000 continue
enddo     ! end of loop over idata
open(uw,file='fresv',access='direct',recl=4*nt)
write(uw,rec=1) weight    ! we have put weighted residues inside this vector
close(uw)
!
! longueurs des rais ... we do expect a variable weighting between rays   
!
allocate(xl_ray(nl))
open(uw,file='frays.xlen',access='direct',recl=4*nl)
read(uw,rec=1) xl_ray
close(uw)
!
!-------------- preconditioning matrix
!     1. diagonal matrices : D1 et D2
!        renormalisation matrix (right multiplication)
!        D2   (dimension nc)    related to the dimension of the model space
!
!        weight (left multiplication)
!        D1   (dimension nl)    related to the dimension of the current dataset
!                     WE MAY HAVE DIFFERENT OPTIONS
!                            1. reading the weight file
!                            2. making automatic statistics on residues
!                            3. both by multiplying contribution
!                            4. weight over values of residues r0 and r1
!                            5. weight over values of rays l0 and l1
!=====================================
!
!---------------------calcul de D1
allocate(d1(nl))                        ! coeff is not yet used dd
call sub_D1(d1,res,poids,nl,r0,r1,xl_ray,l0,l1,coeff)
!
!-------------------CALCUL DE LA MATRICE PRECONDITIONNEE---------------------------
!------------------calcul de D2
allocate(d2(nc))
d2(:)=0.
call sub_D2(d2,nc,vip,vis,sei,to,mat,ic,id,nnz,nl,cp,cs,cpo,cto)
!------------------calcul de D1.A
call mut_d_m(d1,mat,id,nnz,nl)
!------------------calcul de D1.b   ! new updated residues on the fdift file
!                               b=res
call mut_d_v(d1,res,nl)
write(ures,rec=1) res       ! we overwrite these values
close(ures)                 ! on ferme cette unite
!============================= we want to save that for all data
do idata=1,nt    ! loop over data
  weight(idata)=-999.
  do irhs=1,nl   ! loop over rhs vector for association 
    if(id_dat(idata) == idnl(irhs)) then
      weight(idata)=res(irhs)     ! new weighted residues (ray length & residue values)
    endif
  enddo   ! end of loop over rhs
enddo     ! end of loop over idata
open(uw,file='fresw',access='direct',recl=4*nt)
write(uw,rec=1) weight    ! we have put weighted residues inside this vector
close(uw)
!  rms ponderee sur les datas : on apprend le rms apres son calcul
call write_rms(res,d1,nl)   ! toujours l'unite 49 dans cette subroutine
!------------------calcul de (D1.A).D2
call mut_m_d(d2,mat,ic,nnz,nc)
! il faut sauvegarder D2 pour ensuite mettre a l'echelle ou renormaliser lors de l'inversion
open(uw,file='fnorm',access='direct',recl=4*nc)
write(uw,rec=1) d2
close(uw)
!
! write the normalized sensitive matrix (without laplacian): we delete the previous one
!
call write_matrix(mat,ic,id,nnz,nl,umat,uic,uid,fmat,fic,fid,flog)
!
! free memory
!
deallocate(mat);deallocate(ic);deallocate(id);deallocate(res);deallocate(weight)
deallocate(id_dat);deallocate(idnl)

if(lap == 0) then
  write(*,*) ' no laplacian smoothing '
  stop       ! fin des modifications sur la matrice A   D2 sera necessaire pour la vraie sol.
endif

!######################### C'est fini avec le scaling/normalisation : on continue seulement si laplacien

!------------------calcul de la matrice de lissage (laplacien)
!================== on alloue les tailles maximales
allocate(lmat(nnz1))
allocate(lic(nnz1))
allocate(lid(nl1+1))
allocate(lres(nl1))

call sub_lap(lap,lmat,lic,lid,lres,nnz1,nl1,lx,ly,lz,nnz,nl,nix,niy,niz,ondesS)

lumat=iunit
iunit=iunit+1
luic=iunit
iunit=iunit+1
luid=iunit
iunit=iunit+1
lfmat='lfmat.x'
lfic='lfmat.ic'
lfid='lfmat.id'
!------------------ comme on a modifie nnz1 et nl1   pas de coherence et donc via _check
call write_matrix_check(lmat,lic,lid,nnz1,nl1,lumat,luic,luid,lfmat,lfic,lfid,flog)
deallocate(lmat)
deallocate(lic)
deallocate(lid)

open(ures,file='lfdift',access='direct',recl=4*nl1)
write(ures,rec=1) (lres(i),i=1,nl1)    
close(ures)
deallocate(lres)

!================= on va allouer maintenant les bonnes tailles totales
nnz_tot=nnz+nnz1
nl_tot=nl+nl1
allocate(mat(nnz_tot))
allocate(ic(nnz_tot))
allocate(id(nl_tot+1))
allocate(res(nl_tot))
! on lit d'abord la premiere partie
nnz_off=0
nl_off=0
call read_matrix_check(mat,ic,id,nnz,nl,nc,nnz_off,nl_off,umat,uic,uid,fmat,fic,fid,flog)

open(ures,file='fdift',access='direct',recl=4*nl)
read(ures,rec=1) (res(i),i=1,nl)    
close(ures)

! on lit la deuxieme partie (laplacien)
nnz_off=nnz
nl_off=nl
call read_matrix_check(mat,ic,id,nnz1,nl1,nc,nnz_off,nl_off,umat,uic,uid,lfmat,lfic,lfid,flog)

open(ures,file='lfdift',access='direct',recl=4*nl1)
read(ures,rec=1) (res(i),i=nl+1,nl_tot)    
close(ures)

! on ecrit l'ensemble sur les fichiers standards

call write_matrix(mat,ic,id,nnz_tot,nl_tot,umat,uic,uid,fmat,fic,fid,flog)

open(ures,file='fdift',access='direct',recl=4*(nl+nl1))    ! on augmente la longueur du record
write(ures,rec=1) res    
close(ures)

! updating the matrix.par file

  write(flog,*) '---------------------------------------------------------------'
  write(flog,'(35x,a)') ' writing matrix.par  '
  open(49,file='matrix.par')
  write(49,'(a)') fmat
  write(49,'(a)') fic
  write(49,'(a)') fid
  write(49,*) nnz_tot
  write(49,*) nl_tot,nlp                    ! nlp ne sert a rien ici
  write(49,*) nc                            ! column
  close(49)

deallocate(mat)
deallocate(ic)
deallocate(id)
deallocate(res)

stop
end program precond_smooth









