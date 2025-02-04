!==================================================================
! provide from coarse to fine sampling  <<<<
! for a FMG strategy
! for Jamstec, the sampling is 2**10+1 (1025), 2*8+1 (257) (x,z)
!
!==================================================================
program sampling_fine

real(kind=4),dimension(:,:),allocatable :: coarse_grid,fine_grid
integer(kind=4) :: ix,iz,nx_fine,nz_fine,nx_coarse,nz_coarse,kx,kz

real(kind=4) :: dx,dz,xtot,ztot,dx2,dz2

character(len=132) :: name_coarse,name_fine

write(*,*) ' enter the file name of the coarse grid '
read(*,'(a)') name_coarse
write(*,*) ' enter the file name of the fine grid '
read(*,'(a)') name_fine
write(*,*) ' enter the level kx,kz of the initial coarse grid '
write(*,*) ' for Jamstec, it is 10 and 8 '
read(*,*) kx,kz
write(*,*) ' enter the step lengths dx,dz of the initial coarse grid '
read(*,*) dx,dz

!========================= coarse

nx_coarse=2**kx+1
nz_coarse=2**kz+1
xtot=dx*float(nx_coarse-1)
ztot=dz*float(nz_coarse-1)

write(*,*) ' coarse grid dimensions ',nx_coarse,nz_coarse
write(*,*) ' total size ',xtot,ztot
write(*,*) ' sampling ',dx,dz

nx_fine=2**(kx+1)+1
nz_fine=2**(kz+1)+1
dx2=xtot/float(nx_fine-1)
dz2=ztot/float(nz_fine-1)


write(*,*) ' fine grid dimensions ',nx_fine,nz_fine
write(*,*) ' total size ',xtot,ztot
write(*,*) ' sampling ',dx2,dz2


open(7,file=name_coarse,access='direct',recl=4*nx_coarse*nz_coarse)
open(8,file=name_fine,access='direct',recl=4*nx_fine*nz_fine)

allocate(coarse_grid(nx_coarse,nz_coarse))
allocate(fine_grid(nx_fine,nz_fine))

fine_grid(:,:)=0.

!=============================== read it normally
read(7,rec=1) coarse_grid

call coarse2fine(coarse_grid,fine_grid,nx_coarse,nz_coarse,nx_fine,nz_fine)

vmin=1.e+29
vmax=-1.e+29

do iz=1,nz_fine
do ix=1,nx_fine
if(fine_grid(ix,iz) > vmax) vmax=fine_grid(ix,iz)
if(fine_grid(ix,iz) < vmin) vmin=fine_grid(ix,iz)
if(abs(fine_grid(ix,iz)) < 0.001) write(*,*) ix,iz,fine_grid(ix,iz)
enddo
enddo

write(*,*) ' vmin, vmax ',vmin,vmax

write(8,rec=1) fine_grid

close(8)
close(7)

deallocate(coarse_grid)
deallocate(fine_grid)

stop
end program sampling_fine

