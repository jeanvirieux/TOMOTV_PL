!###############################################################################|
! Program fsphe_fmeca2focal      : fusion of two informations into a single file
!        
!---fsphe (binary rec=10*4)                                               
!---fmeca (binary rec=10.4)
!---focal (binary rec=10.4)       
!###############################################################################|

program fsphe_fmeca2focal
  implicit none


  integer(kind=4) :: id_dat0,lut_src0,lut_sta0,impulse0,iphase0
  integer(kind=4) :: id_dat,lut_src,lut_sta,impulse,iphase
  real(kind=4) :: px0,py0,pz0,azimuth0,dip0
  real(kind=4) :: px,py,pz,azimuth,dip

  integer(kind=4) :: nt0,ntp0,nts0,irec
  integer(kind=4) :: nt,ntp,nts,ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7

  open(10,file='fsphe',access='direct',recl=10*4)
  open(11,file='fmeca',access='direct',recl=10*4)
  open(12,file='focal',access='direct',recl=10*4)

  read(10,rec=1) nt0,ntp0,nts0,ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7
  read(11,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7

  if(nt0 /= nt) then
     write(*,*)  'error in the total number of picks ',nt0,nt
     write(*,*) ' please check fsph and fmeca files '
     stop
  endif

  if(ntp0 /= ntp) then
     write(*,*)  'error in the total number of P picks ',ntp0,ntp
     write(*,*) ' please check fsph and fmeca files '
     stop
  endif

  if(nts0 /= nts) then
     write(*,*)  'error in the total number of S picks ',nts0,nts
     write(*,*) ' please check fsph and fmeca files '
     stop
  endif

  write(12,rec=1) nt,ntp,nts,ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7

  do irec=1,nt
     read(10,rec=irec+1) id_dat0,lut_src0,lut_sta0,px0,py0,pz0,azimuth0,dip0,impulse0,iphase0
     read(11,rec=irec+1) id_dat,lut_src,lut_sta,px,py,pz,azimuth,dip,impulse,iphase

     if(id_dat /= id_dat0) then
        write(*,*) ' error in the data IDs: they are different',id_dat,id_dat0
        write(*,*) ' please check fsph and fmeca files '
        stop
     endif

     if(lut_src /= lut_src0) then
        write(*,*) ' error in the event IDs: they are different',lut_src,lut_src0
        write(*,*) ' please check fsph and fmeca files '
        stop
     endif

     if(lut_sta /= lut_sta0) then
        write(*,*) ' error in the station IDs: they are different',lut_sta,lut_sta0
        write(*,*) ' please check fsph and fmeca files '
        stop
     endif

     !======================== make the fusion between picking information and ray information
     !                                             impulse0,iphase0                  px,py,pz
     write(12,rec=irec+1) id_dat0,lut_src0,lut_sta0,px,py,pz,azimuth,dip,impulse0,iphase0

  enddo

  close(9)
  close(10)
  close(11)
  close(12)

  stop  
end program fsphe_fmeca2focal
