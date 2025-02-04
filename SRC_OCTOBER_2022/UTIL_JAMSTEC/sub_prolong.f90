!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! simple Full Multi Grid   FMG ...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!===============================================
! restriction operator
!    from fine to coarse grids
!
! dimensions    n=2**k+1    3,5,9,17,33 etc
! field_2h coarse grid
! field_h  fine grid
!===============================================
subroutine fine2coarse(field_h,field_2h,nx_h,ny_h,nx_2h,ny_2h)
implicit none

real(kind=4), dimension(nx_h,ny_h) :: field_h
real(kind=4), dimension(nx_2h,ny_2h) :: field_2h
integer(kind=4) :: nx_h,ny_h,nx_2h,ny_2h

integer(kind=4) :: ix,iy,ix2,iy2

do iy=2,ny_2h-1
do ix=2,nx_2h-1
ix2=(ix-1)*2+1
iy2=(iy-1)*2+1
field_2h(ix,iy)=0.25*field_h(ix2,iy2) &
              +0.125*(field_h(ix2+1,iy2)+field_h(ix2-1,iy2)+field_h(ix2,iy2-1)+field_h(ix2,iy2+1))  &
              +0.0625*(field_h(ix2+1,iy2+1)+field_h(ix2+1,iy2-1)+field_h(ix2-1,iy2+1)+field_h(ix2-1,iy2-1))
enddo
enddo
do ix=2,nx_2h-1
ix2=(ix-1)*2+1
field_2h(ix,1)=0.5*(field_h(ix2,1)+field_h(ix2+1,1))
field_2h(ix,ny_2h)=0.5*(field_h(ix2,ny_h)+field_h(ix2+1,ny_h))
enddo
do iy=2,ny_2h-1
iy2=(iy-1)*2+1
field_2h(1,iy)=0.5*(field_h(1,iy2)+field_h(1,iy2+1))
field_2h(nx_2h,iy)=0.5*(field_h(nx_h,iy2)+field_h(nx_h,iy2+1))
enddo
field_2h(1,1)=field_h(1,1)
field_2h(nx_2h,1)=field_h(nx_h,1)
field_2h(1,ny_2h)=field_h(1,ny_h)
field_2h(nx_2h,ny_2h)=field_h(nx_h,ny_h)
return
end subroutine fine2coarse


!===============================================
! prolongation operator
!    from coarse to fine grids
!
! dimensions    n=2**k+1    3,5,9,17,33 etc
!
! field_2h coarse grid
! field_2 fine grid
!===============================================
subroutine coarse2fine(field_2h,field_h,nx_2h,ny_2h,nx_h,ny_h)
implicit none

real(kind=4),dimension(nx_h,ny_h) :: field_h
real(kind=4),dimension(nx_2h,ny_2h) :: field_2h
integer(kind=4) :: nx_h,ny_h,nx_2h,ny_2h

integer(kind=4) :: ix,iy,ix2,iy2
!====================== injection
do iy=3,ny_h-1
do ix=3,nx_h-1
ix2=(ix-1)/2+1
iy2=(iy-1)/2+1
field_h(ix,iy)=field_2h(ix2,iy2)            ! take the value
field_h(ix+1,iy)=0.5*(field_2h(ix2+1,iy2)+field_2h(ix2-1,iy2))   ! average over x
field_h(ix,iy+1)=0.5*(field_2h(ix2,iy2+1)+field_2h(ix2,iy2-1))   ! average over y
field_h(ix+1,iy+1)=0.25*(field_2h(ix2+1,iy2+1)+field_2h(ix2-1,iy2+1)+field_2h(ix2-1,iy2-1)+field_2h(ix2+1,iy2-1)) ! corner average
enddo
enddo
! left  right
do ix=1,nx_h-1
ix2=(ix-1)/2+1
field_h(ix,1)=field_2h(ix2,1)
field_h(ix+1,1)=0.5*(field_2h(ix2+1,1)+field_2h(ix2,1))
field_h(ix,2)=field_2h(ix2,1)
field_h(ix+1,2)=0.5*(field_2h(ix2+1,1)+field_2h(ix2,1))
field_h(ix,ny_h)=field_2h(ix2,ny_2h)
field_h(ix+1,ny_h)=0.5*(field_2h(ix2+1,ny_2h)+field_2h(ix2,ny_2h))
enddo
!bottom up
do iy=1,ny_h-1
iy2=(iy-1)/2+1
field_h(1,iy)=field_2h(1,iy2)
field_h(1,iy+1)=0.5*(field_2h(1,iy2+1)+field_2h(1,iy2))
field_h(2,iy)=field_2h(1,iy2)
field_h(2,iy+1)=0.5*(field_2h(1,iy2+1)+field_2h(1,iy2))
field_h(nx_h,iy)=field_2h(nx_2h,iy2)
field_h(nx_h,iy+1)=0.5*(field_2h(nx_2h,iy2)+field_2h(nx_2h,iy2))
enddo
! corners
field_h(1,1)=0.5*field_2h(1,1)+0.25*(field_2h(2,1)+field_2h(1,2))
field_h(nx_h,1)=0.5*field_2h(nx_2h,1)+0.25*(field_2h(nx_2h-1,1)+field_2h(nx_2h,2))
field_h(nx_h,ny_h)=0.5*field_2h(nx_2h,ny_2h)+0.25*(field_2h(nx_2h-1,ny_2h)+field_2h(nx_2h,ny_2h-1))
field_h(1,ny_h)=0.5*field_2h(1,ny_2h)+0.25*(field_2h(2,ny_2h)+field_2h(1,ny_2h-1))
return
end subroutine coarse2fine
