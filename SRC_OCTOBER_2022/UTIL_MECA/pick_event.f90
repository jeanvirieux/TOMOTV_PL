!################################################################################################|
! Program pick-fsrc.f90                                                                          |
!                                                                                                |
! reading the stdout of NLLOC for making the link between an event and the related picking file  |
! first string is the location name and second string is the related picking file                |
!
! Location.20160810.212111(.grid0)  picking_file_name
!
!================================================================================================|
program pick_event

  character(len=256) :: string,string1,substring(14)
  character(len=256) :: pick_string,event_string

  integer(kind=4) :: iflag

  open(7,file='NLLoc.out',status='old')
  open(8,file='pick_event.txt',status='unknown')

  iflag=0
  iassoc=0
100 continue
  read(7,'(a)',end=200) string

  call pick_line_space(string,substring)

  write(*,*) iflag
  write(*,*) trim(substring(3)),' ',trim(substring(4))
  write(*,*) trim(substring(8)),' ',trim(substring(10))

  string=substring(3)
  string1=substring(4)
  if(trim(string) == 'observation' .and. trim(string1) == 'file') then
     string=substring(5)
     pick_string=string(3:len(trim(string)))
     iflag=1
  endif


  string=substring(10)
  string1=substring(8)
  if(trim(string) == 'location' .and. trim(string1) == 'used') then
     string=substring(11)
     event_string=string(25:len(trim(string))-2)
     iflag=0
     iassoc=iassoc+1
     write(*,*) 'new association',iassoc
     write(8,*) trim(event_string),' ',trim(pick_string)
  endif
  goto 100
200 continue

  close(7);close(8)
  stop
end program pick_event
