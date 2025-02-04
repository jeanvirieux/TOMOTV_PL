!==================================================================
! very simple test for checking weights through homogeneous values
! on a FMG strategy
!==================================================================
program test_fmg

real(kind=4),dimension(:,:),allocatable :: coarse_grid,fine_grid
integer(kind=4) :: ix,iy,nx_fine,ny_fine,nx_coarse,ny_coarse,kx,ky


write(*,*) ' enter the level kx,ky '
read(*,*) kx,ky
nx_coarse=2**kx+1
nx_fine=2**(kx+1)+1
ny_coarse=2**ky+1
ny_fine=2**(ky+1)+1

write(*,*) ' coarse grid dimensions ',nx_coarse,ny_coarse
write(*,*) ' fine grid dimensions ',nx_fine,ny_fine

allocate(coarse_grid(nx_coarse,ny_coarse))
allocate(fine_grid(nx_fine,ny_fine))

coarse_grid(:,:)=0.
fine_grid(:,:)=0.

!================================= from fine to coarse grids

write(*,*) ' test from fine to coarse '
do iy=1,ny_fine
do ix=1,nx_fine
fine_grid(ix,iy)=1000. ! ix+iy     !1000
enddo
enddo
call fine2coarse(fine_grid,coarse_grid,nx_fine,ny_fine,nx_coarse,ny_coarse)
open(8,file='coarsegrid',access='direct',recl=4*nx_coarse*ny_coarse)
write(8,rec=1) coarse_grid

do iy=1,ny_coarse
do ix=1,nx_coarse
!if(abs(coarse_grid(ix,iy)-float(ix+iy)) < eps) write(*,*) ix,iy,coarse_grid(ix,iy)
if(abs(coarse_grid(ix,iy)-1000.) < eps) write(*,*) ix,iy,coarse_grid(ix,iy)
enddo
enddo
close(8)

!================================= from coarse to fine grids 

write(*,*) ' test from coarse to fine '

do iy=1,ny_coarse
do ix=1,nx_coarse
coarse_grid(ix,iy)=1000. !  ix+iy     ! 1000 
enddo
enddo
call coarse2fine(coarse_grid,fine_grid,nx_coarse,ny_coarse,nx_fine,ny_fine)
write(*,*) ' OK  '
open(8,file='finegrid',access='direct',recl=4*nx_fine*ny_fine)
write(8,rec=1) fine_grid
do iy=1,ny_fine
do ix=1,nx_fine
!if(abs(fine_grid(ix,iy)-float(ix+iy)) < eps) write(*,*) ix,iy,fine_grid(ix,iy)
if(abs(fine_grid(ix,iy)-1000.) < eps) write(*,*) ix,iy,fine_grid(ix,iy)
enddo
enddo
close(8)


deallocate(coarse_grid)
deallocate(fine_grid)

stop
end program test_fmg

