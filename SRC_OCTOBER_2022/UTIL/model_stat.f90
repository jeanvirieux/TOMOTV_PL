!=======================================================================
! From a suite of models, compute mean and rms
!***********************************************************************
! using a selection of best models if needed ...
!***********************************************************************
program model_mean

  implicit none

  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: model,model_ini
  REAL(kind=4) ,ALLOCATABLE, DIMENSION(:,:,:), TARGET :: mod_mean,mod_mean_ini,mod_rms,mod_rms_ini

  !----------------------------- variables
  real(kind=4) h1inv,h2inv,h3inv,x1inv,x2inv,x3inv ! spatial steps and origin position for the inversion grid
  integer(kind=4) n1inv,n2inv,n3inv,chx   ! dimensions for inversion grid 
  integer(kind=4) flog,finput             ! flux entree et sortie

  integer(kind=4)                    :: i1,i2,i3
  character(len=4)                   :: carac

  integer(kind=4) :: nmod,imod,iopt,iter,jmod

  real(kind=4) :: fcost,fcost_seuil,fcost_ini
  real(kind=4) :: misfit_ini,misfit_final

  logical :: file_exists

  finput=51
  ! #################################### choix des ondes
  open(finput,file='model.head',status='old')
  read(finput,'(a)') carac
  read(finput,*) chx
  ! ####################################  input
  read(finput,'(a)') carac
  read(finput,*) x1inv,x2inv,x3inv
  read(finput,'(a)') carac
  read(finput,*) n1inv,n2inv,n3inv
  read(finput,'(a)') carac
  read(finput,*) h1inv,h2inv,h3inv
  close(finput)

  write(*,*) ' enter number of models '
  read(*,*) nmod
  write(*,*) 'enter P (1) or S (2) '
  read(*,*) iopt
  write(*,*) ' enter the fcost threshold '
  read(*,*) fcost_seuil


  ALLOCATE(mod_mean(n1inv,n2inv,n3inv))
  ALLOCATE(mod_mean_ini(n1inv,n2inv,n3inv))
  ALLOCATE(mod_rms(n1inv,n2inv,n3inv))
  ALLOCATE(mod_rms_ini(n1inv,n2inv,n3inv))

  ALLOCATE(model(n1inv,n2inv,n3inv))
  ALLOCATE(model_ini(n1inv,n2inv,n3inv))

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extraction des misfits initiaux et misfits finaux
  open(9,file='rms.tot',status='unknown')
  open(10,file='misfit.ini',status='unknown')
  open(11,file='misfit.final',status='unknown')

  do imod=1,nmod
     write(carac,'(I4.4)') imod
     inquire(file='MODEL.FIN/rms.'//carac//'.pond',exist=file_exists)
     if(file_exists) then
        open(7,file='MODEL.FIN/rms.'//carac//'.pond',status='old')
        read(7,*) misfit_ini
10      continue
        read(7,*,end=11) misfit_final
        goto 10
11      continue
        close(7)
        write(9,*)  imod,misfit_ini,misfit_final
        write(10,*) imod,misfit_ini
        write(11,*) imod,misfit_final
     else
        write(9,*)  imod,-999.,-999.    ! data is missing
        write(10,*)  imod,-999.
        write(11,*)  imod,-999.
     endif
  enddo
  close(9)
  close(10)
  close(11)


!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ fonction cout par modele  
  open(9,file='rms.tot',status='old')

  mod_mean(:,:,:)=0.
  mod_mean_ini(:,:,:)=0.
  mod_rms(:,:,:)=0.
  mod_rms_ini(:,:,:)=0.
  jmod=0

  do imod=1,nmod

     read(9,*) iter,fcost_ini,fcost

     if(fcost < fcost_seuil) then
        jmod=jmod+1
        write(carac,'(I4.4)') imod

        if(iopt == 1) inquire(file='MODEL.FIN/modelP.'//carac//'.21',exist=file_exists)
        if(iopt == 2) inquire(file='MODEL.FIN/modelS.'//carac//'.21',exist=file_exists)
        
        if(file_exists) then
           ! #####################################  io
           if(iopt == 1) then
              open(7,file='MODEL.FIN/modelP.'//carac//'.21',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
              open(8,file='MODEL.INI/modelP.'//carac//'.ini',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
           endif
           !------------ if we have selected P & S velocities simultaneously
           if(iopt == 2) then
              open(7,file='MODEL.FIN/modelS.'//carac//'.21',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
              open(8,file='MODEL.INI/modelS.'//carac//'.ini',access='direct',recl=4*n1inv*n2inv*n3inv,status='old')
           endif

           ! #####################################  reading the inverse model
           read(7,rec=1) model             ! true grid in P
           read(8,rec=1) model_ini
           close(7)

           mod_mean(:,:,:)=mod_mean(:,:,:)+model(:,:,:)
           mod_mean_ini(:,:,:)=mod_mean_ini(:,:,:)+model_ini(:,:,:)
           mod_rms(:,:,:)=mod_rms(:,:,:)+model(:,:,:)**2
           mod_rms_ini(:,:,:)=mod_rms_ini(:,:,:)+model_ini(:,:,:)**2
        else
           write(*,*) ' data missing for ',imod
        endif

     endif

  enddo   ! imod

  close(9)

  write(*,*) ' number of selected models ',jmod

  if(jmod > 0) then

     mod_mean(:,:,:)=mod_mean(:,:,:)/float(jmod)
     mod_mean_ini(:,:,:)=mod_mean_ini(:,:,:)/float(jmod)
     if(jmod > 1) then
        mod_rms_ini(:,:,:)=sqrt(mod_rms_ini(:,:,:)-float(jmod)*mod_mean_ini(:,:,:)**2)/float(jmod-1)
     else
        mod_rms_ini(:,:,:)=0.
     endif
     
     mod_mean_ini(:,:,:)= (mod_mean(:,:,:)-mod_mean_ini(:,:,:))/mod_mean_ini(:,:,:)*100. !  ! *100000.    ! /mod_mean(:,:,:)*100.

     if(jmod > 1) then
        mod_rms(:,:,:)=sqrt(mod_rms(:,:,:)-float(jmod)*mod_mean(:,:,:)**2)/float(jmod-1)
        mod_rms(:,:,:)=mod_rms(:,:,:)/(mod_rms_ini(:,:,:)+0.001)
     else
        mod_rms(:,:,:)=0.
     endif
     if(iopt == 1) open(7,file='modelP.mean',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     if(iopt == 2) open(7,file='modelS.mean',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(7,rec=1) mod_mean
     close(7)

     if(iopt == 1) open(7,file='modelP_per.mean',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     if(iopt == 2) open(7,file='modelS_per.mean',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(7,rec=1) mod_mean_ini
     close(7)

     if(iopt == 1) open(7,file='modelP_ini.rms',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     if(iopt == 2) open(7,file='modelS_ini.rms',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(7,rec=1) mod_rms_ini
     close(7)

     if(iopt == 1) open(7,file='modelP_nor.rms',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     if(iopt == 2) open(7,file='modelS_nor.rms',access='direct',recl=4*n1inv*n2inv*n3inv,status='unknown')
     write(7,rec=1) mod_rms
     close(7)

     write(*,*) ' Min/max for model and rms'

     write(*,*) 'model (km/s) ',minval(mod_mean),maxval(mod_mean)
     write(*,*) 'model percentage',minval(mod_mean_ini),maxval(mod_mean_ini)

     write(*,*) 'model rms ini ',minval(mod_rms_ini),maxval(mod_rms_ini)
     write(*,*) 'model rms normalized ',minval(mod_rms),maxval(mod_rms)

  else

     write(*,*) ' no model found '

  endif

  deallocate(model,model_ini,mod_mean,mod_mean_ini,mod_rms,mod_rms_ini)

end program model_mean






