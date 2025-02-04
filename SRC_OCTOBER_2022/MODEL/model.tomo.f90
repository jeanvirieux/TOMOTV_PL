!=======================================================================
!
!   January 2015: making a depth analysis in order to see 
!                 evolution of the inverse model
!                 aside the eikonal fine grid
!
!***********************************************************************
! module : model
! ==============
! purpose : interpolation the inversion model into the forward regular grid for
!           travel-time computation. The forward grid is a cubic grid.
!
! This program builds as well files "timefd.par" and "inversion.par"
!
! One may select only P waves or P & S waves. Grids for P and S velocities
! are identical.
!
! model description : the model is described on a rectangular grid for
!           the moment ..... 
!                      algorithms are selected for possible extension to
!                      unstructured finite element grids
!
! the log file name is "flog.model" 
! the input file is "model.head"      (march 2016) for Jenny Jacques
!
!-----------------------------------------------------------------------
! rebuilt by Jean Virieux in fortran 90 (summer 2009)
!***********************************************************************
program model
  use s_interpol_mod
  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modinvp,modinvs
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modforwp,modforws
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: modlayerp,modlayers
  !     integer pmodforwp,pmodforws
  !----------------------------- variables
  real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  real(kind=4) h1,h2,h3,x1_fwd,x2_fwd,x3_fwd       ! spatial step for FMM and origin position
  integer(kind=4) n1inv,n2inv,n3inv,n1,n2,n3,chx   ! dimensions for inversion grid and for FMM grid
  integer(kind=4) ixo,iyo,izo,nix,niy,niz          ! grid indexes for the inversion and dimensions of this grid
  integer(kind=4) flog,finput,foutput              ! flux entree et sortie

  !  character(len=132) name_acquiP,name_acquiS       ! acquisition files
  !  character(len=132) src_inv,src_tab               !
  !  character(len=132) name_modelp,name_forwardp     ! filename for model and FMM for P wave
  !  character(len=132) name_models,name_forwards     ! filename for model and FMM for S wave

  real(kind=4)                       :: rap

  real(kind=4)                       :: vpmin,vpmax,vsmin,vsmax,vmoy
  integer(kind=4)                    :: ix,iy,iz,itmoy
  character(len=1)                   :: carac

  integer(kind=4)                    :: iopt=1    !!! harmonic



  flog=51
  open(flog,file='flog.model')
  finput=flog+1
  ! #################################### choix des ondes
  open(finput,file='model.head',status='old',err=1000)
  read(finput,'(a)') carac
  read(finput,*) chx
  write(flog,*) ' flag 1 pour ondes P, 2 pour ondes P et S',chx
  ! ####################################  input
  read(finput,'(a)') carac
  read(finput,*) x1inv,x2inv,x3inv
  write(flog,*) ' (x1,x2,x3) origin of the inversion grid ',x1inv,x2inv,x3inv
  read(finput,'(a)') carac
  read(finput,*) n1inv,n2inv,n3inv
  write(flog,*) ' (n1,n2,n3) inversion grid ',n1inv,n2inv,n3inv
  read(finput,'(a)') carac
  read(finput,*) h1inv,h2inv,h3inv
  write(flog,*) ' (h1,h2,h3) grid sampling ', h1inv,h2inv,h3inv
  !  read(finput,'(a)') carac
  !  read(finput,*) x1_fwd,x2_fwd,x3_fwd
  x1_fwd=x1inv;x2_fwd=x2inv;x3_fwd=x3inv
  !  write(flog,*) ' (x1,x2,x3) origin of the forward grid ',x1_fwd,x2_fwd,x3_fwd
  read(finput,'(a)') carac
  read(finput,*) h1
  h2=h1; h3=h1     ! cube for the forward problem
  write(flog,*) ' cubic grid for the eikonal grid h1=h2=h3 ',h1,h2,h3
  n1=int(((x1inv+(n1inv-1)*h1inv)-x1_fwd)/h1)+1 
  n2=int(((x2inv+(n2inv-1)*h2inv)-x2_fwd)/h2)+1 
  n3=int(((x3inv+(n3inv-1)*h3inv)-x3_fwd)/h3)+1 
  write(flog,*) ' (n1,n2,n3) eikonal grid h1=h2=h3 ',n1,n2,n3
  read(finput,'(a)') carac
  read(finput,*) rap
  write(flog,*) ' ray sampling in meter ',rap
  read(finput,'(a)') carac
  read(finput,*) ixo,iyo,izo
  write(flog,*) ' start indexes of subgrid to be inverted ',ixo,iyo,izo
  read(finput,'(a)') carac
  read(finput,*) nix,niy,niz
  write(flog,*) ' end indexes of subgrid to be inverted ',nix,niy,niz
  !   ixo=1;iyo=1;izo=1;nix=n1inv;niy=n2inv;niz=n3inv
  close(finput)

  ! #####################################  arrays
  ALLOCATE(modinvp(n1inv,n2inv,n3inv))
  ALLOCATE(modlayerp(n1inv,n2inv,n3inv))
  ALLOCATE(modforwp(n1,n2,n3))
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     ALLOCATE(modinvs(n1inv,n2inv,n3inv))
     ALLOCATE(modlayers(n1inv,n2inv,n3inv))
     ALLOCATE(modforws(n1,n2,n3))
  endif

  ! #####################################  io
  open(7,file='modelP',access='direct',recl=4*n1inv*n2inv*n3inv,status='old',err=1000)
  open(17,file='average_layer_vp',access='direct',recl=4*n1inv*n2inv*n3inv)       ! JEAN new file
  open(8,file='modelP.fwd',access='direct',recl=4*n1*n2*n3)
  !------------ if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     open(9,file='modelS',access='direct',recl=4*n1inv*n2inv*n3inv,status='old',err=2000)
     open(19,file='average_layer_vs',access='direct',recl=4*n1inv*n2inv*n3inv)    ! JEAN new file
     open(10,file='modelS.fwd',access='direct',recl=4*n1*n2*n3)
  endif

  ! #####################################  reading the inverse model
  read(7,rec=1) modinvp             ! true grid in P
  ! make sure that the first two nodes are frozen (the acquisition should be such that iz > 2
  ! for example  -7000 m, -5000 m, -3000 m (stepping 2000 m): Mont Blanc at -4807m is OK
  modinvp(:,:,2)=modinvp(:,:,3)     ! always copy the two ghosts lines
  modinvp(:,:,1)=modinvp(:,:,2)     ! 
  write(*,*) ' max and min input P velocities ',maxval(modinvp),minval(modinvp)
  ! ----------- if we have selected P & S velocities simultaneously
  if(chx.eq.2) then
     read(9,rec=1) modinvs             ! true grid in S
     modinvs(:,:,2)=modinvs(:,:,3)     ! always copy the two ghosts lines
     modinvs(:,:,1)=modinvs(:,:,2)     ! 
     write(*,*) ' max and min input  S velocities ',maxval(modinvs),minval(modinvs)
  endif

  ! ####################################  make the interpolation from inversion to FMM
  write(*,*) ' P wave forward grid computation '
  call subinterpol(modinvp,modforwp,x1inv,x2inv,x3inv,n1inv,n2inv,n3inv, &
       h1inv,h2inv,h3inv,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vpmin,vpmax,iopt)
  ! ----------- if we have selected P & S velocities simultaneously
  write(*,*) 'Vp global min, max ',vpmin,vpmax
  if(chx.eq.2) then
     write(*,*) ' S wave forward grid computation '
     call subinterpol(modinvs,modforws,x1inv,x2inv,x3inv,n1inv,n2inv,n3inv, &
          h1inv,h2inv,h3inv,x1_fwd,x2_fwd,x3_fwd,n1,n2,n3,h1,h2,h3,vsmin,vsmax,iopt)
     write(*,*) 'Vs global min, max ',vsmin,vsmax
  endif

  ! ##############################################################
  ! Depth analysis for better understanding what is going on during inversion
  ! ##############################################################

  do iz=1,n3inv
     vpmin=1.e+29
     vpmax=-1.e+29
     vmoy=0.
     itmoy=0
     do iy=1,n2inv
        do ix=1,n1inv
           if(modinvp(ix,iy,iz) > vpmax) vpmax=modinvp(ix,iy,iz)
           if(modinvp(ix,iy,iz) < vpmin) vpmin=modinvp(ix,iy,iz)
           itmoy=itmoy+1
           vmoy=vmoy+modinvp(ix,iy,iz)
!          if(modinvp(ix,iy,iz) > 10000.) write(*,*) ix,iy,iz,modinvp(ix,iy,iz)
        enddo !iy
     enddo !ix
     vmoy=vmoy/float(itmoy-1)
     !============================= set the related average layer model: remove values at edges
     do iy=1,n2inv
        do ix=1,n1inv
           modlayerp(ix,iy,iz)=vmoy
        enddo  ! ix
        modlayerp(1,iy,iz)=modlayerp(2,iy,iz)
        modlayerp(n1inv,iy,iz)=modlayerp(n1inv-1,iy,iz)
     enddo  ! iy
     do ix=1,n1inv
        modlayerp(ix,1,iz)=modlayerp(ix,2,iz)
        modlayerp(ix,n2inv,iz)=modlayerp(ix,n2inv-1,iz)
     enddo ! ix

     write(flog,*) 'Vp: depth wrt the origin (layer)',(iz-1)*h3inv,iz*h3inv
     write(flog,*) 'vpmin, vpmax, vpmean',vpmin,vpmax,vmoy
  enddo ! iz
  modlayerp(:,:,1)=modlayerp(:,:,2)   ! always copy the second line into the first line
  modlayerp(:,:,n3inv)=modlayerp(:,:,n3inv-1) ! bottom is also copied

  ! ############################################################ end of the depth information
  ! #####################################  writing the forward model
!  modforwp(:,:,2)=modforwp(:,:,3)   ! always copy the second line into the first line
!  modforwp(:,:,1)=modforwp(:,:,2) ! bottom is also copied
  write(8,rec=1) modforwp
  write(17,rec=1) modlayerp               ! JEAN new input file ... if needed

  if(chx == 2) then
     ! ##############################################################
     ! Depth analysis for better understanding what is going on during inversion
     ! ##############################################################

     do iz=1,n3inv
        vsmin=1.e+29
        vsmax=-1.e+29
        vmoy=0.
        itmoy=0
        do iy=1,n2inv
           do ix=1,n1inv
              if(modinvs(ix,iy,iz) > vsmax) vsmax=modinvs(ix,iy,iz)
              if(modinvs(ix,iy,iz) < vsmin) vsmin=modinvs(ix,iy,iz)
              itmoy=itmoy+1
              vmoy=vmoy+modinvs(ix,iy,iz)
           enddo !iy
        enddo !ix
        vmoy=vmoy/float(itmoy-1)

        !============================= set the related average layer model remove values at edges
        do iy=2,n2inv-1
           do ix=2,n1inv-1
              modlayers(ix,iy,iz)=vmoy
           enddo ! ix
           modlayers(1,iy,iz)=modlayers(2,iy,iz)
           modlayers(n1inv,iy,iz)=modlayers(n1inv-1,iy,iz)
        enddo ! iy
        do ix=1,n1inv
           modlayers(ix,1,iz)=modlayers(ix,2,iz)
           modlayers(ix,n2inv,iz)=modlayers(ix,n2inv-1,iz)
        enddo ! ix

        write(flog,*) 'Vs: depth wrt origin (layer)',(iz-1)*h3inv,iz*h3inv
        write(flog,*) 'vsmin, vsmax, vsmean',vsmin,vsmax,vmoy
     enddo ! iz
     modlayers(:,:,1)=modlayers(:,:,2)             ! always copy the second line into the first line
     modlayers(:,:,n3inv)=modlayers(:,:,n3inv-1)   ! always copy the second line into the first line
     ! ############################################################ end of the depth information
     ! ----------- if we have selected P & S velocities simultaneously
     write(10,rec=1) modforws
     write(19,rec=1) modlayers            ! JEAN new input file ... if needed
  endif

  ! #################################  writing the parameter file "timefd.par"
  foutput=finput+1
  open(foutput,file='timefd.par',status='unknown')
  write(foutput,'(i1)') chx
  write(foutput,'(I8)') 15000    ! rmax : max number of points along a ray (should be tuned by hand)
  write(foutput,'(3f15.6)') x1_fwd,x2_fwd,x3_fwd
  write(foutput,'(3i6)') n1,n2,n3
  write(foutput,'(3f15.6)') h1,h2,h3
  write(foutput,*) rap,'                 ray sampling '
  write(foutput,'(a)') 'modelP.fwd'
  if(chx.eq.2) then                       ! at the end : for possible modification by hand 
     write(foutput,'(a)') 'modelS.fwd'        ! between P and P&S
  endif
  close(foutput)

  ! ################################# writing the inversion file "inversion.par"
  open(foutput,file='inversion.par',status='unknown')
  write(foutput,'(i1)') chx
  write(foutput,'(3f15.6)') x1inv,x2inv,x3inv
  write(foutput,'(3i4)') n1inv,n2inv,n3inv
  write(foutput,'(3f15.6)') h1inv,h2inv,h3inv
  write(foutput,'(a)') 'modelP'
  write(foutput,*) ixo,iyo,izo
  write(foutput,*) nix,niy,niz
  if(chx.eq.2) then                       ! at the end : for possible modification by hand
     write(foutput,'(a)') 'modelS'           ! between P and P&S
  endif
  close(foutput)
  close(flog)
  stop
  !============================ modelP is missing
1000 continue
  write(flog,*) ' error: missing Vp file '
  close(flog)
  stop
  !============================ modelS is missing
2000 continue
  write(flog,*) ' error: missing Vs file '
  close(flog)
  stop
end program model






