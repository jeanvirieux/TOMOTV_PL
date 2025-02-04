!program lec

!character(len=256) :: string,substring(14)

!write(*,*) ' enter string '
!read(*,'(a)') string
!call pick_line_space(string,substring)
!do i=1,14
!write(*,*) trim(substring(i))
!enddo
!stop
!end program lec

subroutine pick_line_space(string,substring)

character(len=256) :: string,substring(14)
integer(kind=4) :: posBT(2),pos

!crop final whitespaces
string=adjustl(trim(string))

  ! Get first part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(1)=string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get second part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(2)=string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get third part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(3) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get fourth part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(4) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get five part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(5) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get six part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(6) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get seven part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(7) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get eight part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(8) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get nine part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(9) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get ten part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(10) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get eleven part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(11) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get twelve part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(12) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get thirteen part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(13) = string(1:pos)
  string = adjustl(string(pos+1:))

  ! Get fourteen part:
  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(32))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring(14) = string(1:pos)

return
end subroutine pick_line_space
