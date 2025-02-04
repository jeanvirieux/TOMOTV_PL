!###############################################################################|
! program src2fsrc    ASCII >>> fsrc                                            |
!                                          from NLLOC                           |
!                                                                               |
!   input (ascii file)                                                          |
!        src.txt                                                                |
!                                                                               |
!  output (binary file)                                                         |
!                                                                               |
!       fsrc        (sources)    8*4 = 32 bytes for each record                 |
! BINARY FORMAT     nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7                  |
!   record 1+irec   sr_x(irec),sr_y(irec),sr_z(irec),sr_to(irec),   &           |
!                   ki_pos(irec),ki_to(irec),ki_m(irec),ki_id(irec)             |
!                                                                               |
!###############################################################################|
!
!   src.txt is the summary file from Nlloc
!
!###############################################################################

program src2fsrc
  implicit none
  integer(kind=4) :: nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7,irec,inull
  integer(kind=4) :: ki_pos,ki_to,ki_m,ki_id,icount,it
  real(kind=4) :: sr_x,sr_y,sr_z,sr_to
  real(kind=4) :: x_orig,y_orig,scale
  integer(kind=4) :: iyear,imonth,iday,ihour,iminute
  character(len=23) :: date_string
  real(kind=8) :: utime

  character(len=4) :: name
  character(len=256) :: string,string1,string2,substring(14)
  character(len=256) :: event_string

  scale=1000.00 ! conversion from km to meters

  !================================================ we assume that they have been applied for NLLOC
  write(*,*) ' enter origin point (x0,y0) in meters'
  read(*,*) x_orig,y_orig

  write(*,*) ' WARNING: current origin is assumed to be the one used for location '

  !===================================== find the number of located events 
  open(11,file='src.txt',status='old')
  nsrc=0
  inull=0
100 continue
  read(11,'(a)',end=200) string
  call pick_line_space(string,substring)
  string2=substring(2)
  string=string2(25:len(trim(string2))-7)   ! just Location.20160810.203539
  string2=substring(3)
  string1=string2(2:len(trim(string2))-1)   ! just LOCATED
  !======================================== check if it is a located event (two formats expected)
  if(string(1:8) == 'Location' .and. trim(string1) == 'LOCATED') then
     nsrc=nsrc+1
     goto 100
  endif
  !======================================== check if it is a rejected event
  if(string(1:8) == 'Location' .and. trim(string1) == 'REJECTED') then
     inull=inull+1
     write(*,*) ' current number of rejected events ',inull
     goto 100
  endif
  goto 100
200 continue
  close(11)

  write(*,*) ' number of located events by NLLOC ',nsrc
  write(*,*) ' number of rejected events ',inull

  open(11,file='src.txt',status='old')
  open(10,file='fsrc',access='direct',recl=8*4)
  open(12,file='fdate.txt')

  write(*,*) ' file conversion from fsrc.dat into fsrc '
  write(*,*) ' writing in file "fdate" the event id and event date under the utime format in seconds '

  neqks=nsrc   ! we consider only quakes
  nblast=0
  nshot=0
  neqks_out=0
  i5=0;i6=0;i7=0

  write(10,rec=1) nsrc,neqks,nshot,nblast,neqks_out,i5,i6,i7

  write(*,*) '========================================== '
  write(*,*) 'number of sources ',nsrc
  write(*,*) ' number of quakes inside the domain ', neqks
  write(*,*) ' number of quakes outside the domain ', neqks_out
  write(*,*) ' number of quakes ',neqks+neqks_out
  write(*,*) ' number of shots ',nshot
  write(*,*) ' number of blasts ',nblast
  write(*,*) '========================================== '

  do irec=1,nsrc
600  continue
     read(11,'(a)',end=2000) string
     !=========================== detect an event
     call pick_line_space(string,substring)
     string2=substring(2)
     string=string2(25:len(trim(string2))-7)   ! just Location.20160810.203539
     string2=substring(3)
     string1=string2(2:len(trim(string2))-1)   ! just LOCATED
     event_string=string
     if(string(1:8) == 'Location' .and. trim(string1) == 'LOCATED') then
        !==================================== yes one need to read the information
700     continue
        if(string(1:9) /= 'END_NLLOC') then
           read(11,'(a)') string   ! read the string
           string1=string
           call pick_line_space(string1,substring)
           !================================================== detection of source parameters
           if(string(1:10) == 'HYPOCENTER') then
              ! x substring(3); y substring(5), z substring(7), OT substring(9)
              read(substring(3),*) sr_x
              read(substring(5),*) sr_y
              read(substring(7),*) sr_z
              read(substring(9),*) sr_to
              sr_x=sr_x*scale
              sr_y=sr_y*scale
              sr_z=sr_z*scale
           endif  ! hypocenter
           !================================================= detection of the date
           if(string(1:10) == 'GEOGRAPHIC') then
              ! iyear substring(3); imonth substring(4); iday substring(5); ihour substring(6); iminute substring(7);
              read(substring(3),*) iyear
              read(substring(4),*) imonth
              read(substring(5),*) iday
              read(substring(6),*) ihour
              read(substring(7),*) iminute
           endif   !  geographic
        else
           goto 701   ! save and move to next event
        endif   !  end_nnloc
        goto 700
701     continue
        ki_pos=irec   ! position inversion
        ki_to=irec    ! origin time inversion
        ki_m=irec     ! inside
        ki_id=irec    ! id of the event
        sr_to=0.      ! cancel out the origin time (substracted when building fobs - pick2fobs-fsphe.f90)
        write(10,rec=irec+1) sr_x,sr_y,sr_z,sr_to,ki_pos,ki_to,ki_m,ki_id
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  now fill in the file fdate.txt for relation between event and location file related to picking file
!       fdate.txt      event linked to Location.20160813.014438                through summary of Nlloc
!       pick_event.txt Location.20160813.014438 linked to Aug13_014428_nll.txt through output of Nlloc
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !========================== write the data string
        write(date_string,'(I4,".",I2.2,".",I2.2,"-",I2.2,".",I2.2,".",F6.3)') iyear,imonth,iday,ihour,iminute,sr_to
        !========== replace blanks by zeros
        do it=1,23
           if(date_string(it:it) == ' ') date_string(it:it)='0'
        enddo
        call str2utime1(date_string,utime)
        write(12,*) ki_id,int(utime),utime,date_string,'   ',trim(event_string),' '
     else
        goto 600 ! continue reading
     endif
  enddo

2000 continue
  Write(*,*) ' end of scanning the NNLOC event file summary: located events and rejected events ',nsrc,inull
  close(unit=10)
  close(unit=11)
  close(unit=12)
  stop
end program src2fsrc


