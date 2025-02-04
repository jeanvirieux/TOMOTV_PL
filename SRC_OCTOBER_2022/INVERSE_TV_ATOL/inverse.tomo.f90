!     
!     ##########################################################
!     #                  Program INVERSE                       #
!     #      Least-square traveltime residual inversion        #
!     #                   using LSQR method                    #
!     ##########################################################
!     
!     Input parameters :
!     
!     - Parameters of file: inversion.par
!     - Damping parameter (damp)
!
!     Input files :
!      
!     - velocity models or 1./U,1./V parameters files 
!     - fsrc
!     - matrix files & residuals files 
!
!     Output files :
!
!     - updated fsrc & velocity models (or 1./U, 1./V parameters files)
!
!     
!               Least-squared resolution of 
!
!                             Ax=b
!
!      for the seismic travel-time tomography
!      
!      A is the sensitivity matrix (Fréchet derivative matrix)
!      b is the data vector containing travel-time delays to be inverted      
!      x is the parameter vector containing perturbation solution to be
!                              inserted to previous values
! 
!            x = (dU,dV,dposx,dposy,dposz,dto), 
!
!            dU   is the perturbation of parameter U (often slowness or velocity),
!            dV   is the perturbation of parameter V,
!            pos  is the perturbation of quakes positions,
!            to   is the perturbation of origin times of quakes.
!
!     with (depending on the choise of selected parameters done in vpvs2uv.f) : 
!
!                U=Vp or Vp*Vs
!                V=Vs or Vp/Vs
!
!--------------------DIMENSIONS-----------------------------------------------    
!  nnz      :number of non-zero elements of the matrix A                      
!  nc       :number of columns of A                                           
!  nl       :number of lines of A = equal to the number of rays               
!  ni       :number of parameters for velocities or U and V             
!  np       :number of parameters for positions                          
!  nt       :number of parameters for origin times                           
!  neqks    :number of sources to be inverted                                
!                         TO BE CHECKED BECAUSE BLASTS SHOULD BE CONSIDERED
!  nsrc     :number total de sources                                         

!----------------------------------------------------------------------------
!############################################################################
program inversion
use s_matrix
use s_new_parameter
use s_read
use s_lsqr

implicit none

!         nombre de parametres
!         grille d'inversion des vitesses
integer(kind=4) :: n1,n2,n3
integer(kind=4) :: ni1,ni2,ni3,nvit
integer(kind=4) :: ixo,iyo,izo
real(kind=4) :: xo,yo,zo,h1,h2,h3

!  sources
integer(kind=4) :: nsrc
real(kind=4), allocatable,dimension(:)  :: x_src,y_src,z_src      ! geometric position
real(kind=4), allocatable,dimension(:) :: to_src         ! origin time
integer(kind=4), allocatable,dimension(:) :: ki_pos     ! 0 except for eqks where incremented 
integer(kind=4), allocatable,dimension(:) :: ki_to      ! 0 except for eqks and blasts where incremented
integer(kind=4), allocatable,dimension(:) :: ki_m       ! 1 when in and 0 when out the model
integer(kind=4), allocatable,dimension(:) :: id_src     ! id on the sources (FIXED)

!  linear system   Ax=b  sous forme space    b=res    x=sol    se is the standard error

real(kind=4), allocatable,dimension(:) :: mat
integer(kind=4), allocatable,dimension(:) :: ic,id
integer(kind=4) :: nc,nl,nll,nnz,nlp,nnzz

real(kind=4), allocatable,dimension(:) :: res
real(kind=4), allocatable,dimension(:) :: sol,se
character(len=132) :: fmat,fic,fid

! normalisation
real(kind=4), allocatable,dimension(:) :: rnorme

!  model
real(kind=4) ,allocatable,dimension(:,:,:) :: vit

!         parametres LSQR
integer(kind=4) :: flog_lsqr,flog
integer(kind=4) :: istop
real(kind=4) :: atol,btol,conlim,anorm,damp
real(kind=4) :: acond,rnorm,xnorm,arnorm
integer(kind=4) :: itn,itnlim ! parametres LSQR en input
real(kind=4), allocatable,dimension(:) :: work1,work2
logical :: wantse
!         fonction esclave de lsqr
!   external aprod

!      fichiers 
integer(kind=4) :: iunit,uvp,uvs,usrc,irec,i
character(len=132) :: fvp,fvs
integer(kind=4) :: neqks,nblast,nshot,neqks_out,idum1,idum2,idum3

real(kind=4) :: dvp_max,dvs_max,dhori_max,dvert_max,dto_max

real(kind=4) :: solmax,solmin

wantse= .false.

! ............    file units start at 10
iunit=10
flog=iunit
iunit=iunit+1
open(flog,file='flog.inverse.tomo')

!
!.............     lecture des entrees.
!     lecture de inversion.par
!     uvp (and possibly uvs) are stream file related to fvp and fvs
call read_par(iunit,xo,yo,zo,n1,n2,n3,h1,h2,h3,ixo,iyo,izo,ni1,ni2,ni3, &
              nvit,uvp,uvs,fvp,fvs,damp, &
              dvp_max,dvs_max,dhori_max,dvert_max,dto_max,itnlim,flog)
if(ni1 > n1 .or. ni2 > n2  .or. ni3 > n3) then
  write(*,*) ' inversion window should be smaller than the inversion box'
  write(*,*) ' please revise your input of the module MODEL'
  stop
endif
!
! reading information for the inversion
!
write(flog,*) '---------------------------------------------------------------'
write(flog,'(35x,a,a15)') 'LECTURE DE matrix.par'
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

!    lecture matrice sparse a inverser 
call read_matrix(fmat,fic,fid,mat,ic,id,nnz,nl,nc,iunit,flog)
write(*,*) ' success in reading the sparse sensitivity matrix'
write(*,*) ' nbre of rows nl', nl
write(*,*) ' nbre of columns nc',nc
write(*,*) ' nbre of non-zero elements nnz',nnz

!     lecture des residus (seulement ceux qui seront utilises dans l'inversion d'ou nl)
allocate(res(nl))
call read_fdift(iunit,res,nl)
write(*,*) ' success in reading the 2n member time residuals ',nl
write(*,*) ' the matrix A has a total of non-zero elements ',nnz

!     solution
allocate(sol(nc))  ! mixing slowness + positions + origin times
sol(:)=0.0
allocate(se(nc))   ! standard error
se(:)=0.0

!
!.........     inversion du systeme.
!
flog_lsqr=iunit
iunit=iunit+1
open(flog_lsqr,file='flog.lsqr')
write(flog,*) 'appel a lsqr'

nll=nl+1   ! avoid passing operated values
nnzz=nnz
atol   = 1.e-12    ! lsqr parameters    au lieu de 1.e-8
btol   = 1.e-6     ! idem
conlim = 1.e4      ! idem


allocate(work1(nc))
allocate(work2(nc))
work2(:)=0.0
work1(:)=0.0

!---------------- appel de lsqr

call lsqr(nl,nc,damp,wantse, &
          nnz,  nll,  nnzz, ic,id,mat,res,work1,work2,sol,se, &
          atol,btol,conlim,itnlim,flog_lsqr,istop,itn,anorm,acond,rnorm,arnorm,xnorm)

!subroutine lsqr(m,n,damp,wantse, &
!          leniw,lenjw,lenrw,iw,jw,rw, u,  v,    w,    x,  se,&
!         atol,btol,conlim,itnlim,nout     ,istop,anorm,acond,rnorm,arnorm,xnorm)

deallocate(work1)
deallocate(work2)

solmax=-1.e19
solmin=1.e+19
do i=1,nc
if(sol(i) > solmax) solmax=sol(i)
if(sol(i) < solmin) solmin=sol(i)
enddo

write(*,*) ' MAX ET MIN DE LA SOLUTION TOTALE ',solmax,solmin

!
!...........     ecriture des sorties.
!
!     lecture de la matrice de renormalisation
write(*,*) ' reading the original values for updating these values '
write(*,*) ' first renormalisation to be done ' 
allocate(rnorme(nc))
call read_renorm(iunit,rnorme,nc)
write(*,*)'fnorm:',rnorme(1),rnorme(nc)
!     renormalisation du resultat brut sorti de lsqr.  
call subrenorm(sol,rnorme,nc)
deallocate(rnorme)
write(flog,*) ' renormalisation ok'
!
! ................ attention aux offsets dans la solution sol (slowness, positions and origin times)
!
!     update the new model  - first read the current model
allocate(vit(n1,n2,n3))
open(uvp,file=fvp,access='direct',recl=4*n1*n2*n3)
read(uvp,rec=1) vit
!     update now P or U parameter 
if(nvit > 0) then      ! modele P ou U  attention il faut decipherer sol
  write(*,*) ' updating the P model'
  call sub_new_model(vit,sol,n1,n2,n3,ni1,ni2,ni3,ixo,iyo,izo,dvp_max,0)

!###############################################################################
! overlap boundary values
!###############################################################################
  vit(1,:,:)=vit(2,:,:)             ! left
  vit(n1,:,:)=vit(n1-1,:,:)         ! right
  vit(:,1,:)=vit(:,2,:)             ! front
  vit(:,n2,:)=vit(:,n2-1,:)         ! back
  vit(:,:,1)=vit(:,:,2)             ! top
  vit(:,:,n3)=vit(:,:,n3-1)         ! bottom
  write(uvp,rec=1) vit
  write(*,*) ' ending the updating of the P model'
endif
close(uvp)

!     update now S or V parameter
if(uvs /= 0 .and. nvit > 0) then ! modele S ou V   prendre la deuxieme partie de sol
  open(uvs,file=fvs,access='direct',recl=4*n1*n2*n3)
  read(uvs,rec=1) vit
  write(*,*) ' updating the S model '
  call sub_new_model(vit,sol,n1,n2,n3,ni1,ni2,ni3,ixo,iyo,izo,dvs_max,uvs)

!###############################################################################
! overlap boundary values
!###############################################################################
  vit(1,:,:)=vit(2,:,:)             ! left
  vit(n1,:,:)=vit(n1-1,:,:)         ! right
  vit(:,1,:)=vit(:,2,:)             ! front
  vit(:,n2,:)=vit(:,n2-1,:)         ! back
  vit(:,:,1)=vit(:,:,2)             ! top
  vit(:,:,n3)=vit(:,:,n3-1)         ! bottom

  write(uvs,rec=1) vit
  write(*,*) ' ending the updating of the S model'
  close(uvs)
endif
deallocate(vit)

!   reading fsrc before updating it
usrc=iunit
iunit=iunit+1
open(usrc,file='fsrc',access='direct',recl=8*4)
read(usrc,rec=1) nsrc,neqks,nshot,nblast,neqks_out,idum1,idum2,idum3

allocate(x_src(nsrc))      ! geometric position
allocate(y_src(nsrc))
allocate(z_src(nsrc))
allocate(to_src(nsrc))        ! origin time
allocate (ki_pos(nsrc))   ! flag on position 0 except for eqks where increment 
allocate (ki_to(nsrc))    ! flag on time 0 except for eqks and shots where increment
allocate (ki_m(nsrc))     ! flag on the in/out of the event
allocate (id_src(nsrc))   ! id of the event   (FIXED during the experiment)

do irec=1,nsrc
  read(usrc,rec=irec+1) x_src(irec),y_src(irec),z_src(irec),to_src(irec), &
                        ki_pos(irec),ki_to(irec),ki_m(irec),id_src(irec)
enddo

!     update fsrc ( should define new eqks inside the model)
write(*,*) ' updating quakes positions '
call sub_new_fsrc(x_src,y_src,z_src,to_src,ki_pos,ki_to,ki_m,id_src,sol, &
                  nsrc,usrc,neqks,nshot,nblast,neqks_out,     &
                  xo,yo,zo,n1,n2,n3,h1,h2,h3,ni1,ni2,ni3, &
                  dhori_max,dvert_max,dto_max,uvs,flog)
close(usrc)

deallocate(x_src)
deallocate(y_src)
deallocate(z_src)    
deallocate(to_src)    
deallocate(ki_pos) 
deallocate(ki_to)  
deallocate(ki_m)   
deallocate(id_src) 
stop
end program inversion

