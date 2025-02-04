!==================================================================
! provide the different sampling in 2**k +1 from an initial sampling
! for a FMG strategy
! for Jamstec, the sampling is 2**10+1 (1025), 2*8+1 (257) (x,z)
!==================================================================
program sampling_coarse

real(kind=4),dimension(:,:),allocatable :: coarse_grid,fine_grid
integer(kind=4) :: ix,iz,nx_fine,nz_fine,nx_coarse,nz_coarse,kx,kz

real(kind=4) :: dx,dz,xtot,ztot,dx2,dz2

character(len=132) :: name_coarse,name_fine

write(*,*) ' enter the file name of the fine grid '
read(*,'(a)') name_fine
write(*,*) ' enter the file name of the coarse grid '
read(*,'(a)') name_coarse
write(*,*) ' enter the level kx,kz of the initial grid'
write(*,*) ' for Jamstec, it is 10 and 8 '
read(*,*) kx,kz
write(*,*) ' enter the step lengths dx,dz of the initial grid '
read(*,*) dx,dz


nx_fine=2**kx+1
nz_fine=2**kz+1
xtot=dx*float(nx_fine-1)
ztot=dz*float(nz_fine-1)

write(*,*) ' fine grid dimensions ',nx_fine,nz_fine
write(*,*) ' total size ',xtot,ztot
write(*,*) ' sampling ',dx,dz

nx_coarse=2**(kx-1)+1
nz_coarse=2**(kz-1)+1
dx2=xtot/float(nx_coarse-1)
dz2=ztot/float(nz_coarse-1)


write(*,*) ' coarse grid dimensions ',nx_coarse,nz_coarse
write(*,*) ' total size ',xtot,ztot
write(*,*) ' sampling ',dx2,dz2

open(7,file=name_fine,access='direct',recl=4*nx_fine*nz_fine)
open(8,file=name_coarse,access='direct',recl=4*nx_coarse*nz_coarse)


allocate(coarse_grid(nx_coarse,nz_coarse))
allocate(fine_grid(nx_fine,nz_fine))

coarse_grid(:,:)=0.

read(7,rec=1) fine_grid

call fine2coarse(fine_grid,coarse_grid,nx_fine,nz_fine,nx_coarse,nz_coarse)

write(8,rec=1) coarse_grid

close(8)
close(7)

deallocate(coarse_grid)
deallocate(fine_grid)

stop
end program sampling_coarse

