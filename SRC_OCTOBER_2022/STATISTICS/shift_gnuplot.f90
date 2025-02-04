!###############################################################################|
! PROGRAM FOR STATISTICS ANALYSIS BETWEEN TWO SETS OF SOURCES                   | 
!                                                                               |
!  COMPUTATION OF RESIDUALS DISTRIBUTION FOR EACH STATION                       |
!     input                                                                     |
!                                                                               |
!  TWO fsrc files   fsrc.one and fsrc.two                                       |
!                                                                               |
!                                                                               |
!###############################################################################|
program source_analysis
use s_read_src


implicit none

real(kind=4), allocatable, dimension(:) :: sr1_x,sr1_y,sr1_z,sr1_to
real(kind=4), allocatable, dimension(:) :: sr2_x,sr2_y,sr2_z,sr2_to
integer(kind=4), allocatable, dimension(:) :: sr1_id,sr2_id
integer(kind=4), allocatable, dimension(:) :: ki_pos,ki_to,ki1_m,ki2_m,sr2_twin

real(kind=4), allocatable, dimension(:) :: sr_dx,sr_dy,sr_dz,sr_dt
real(kind=4), allocatable, dimension(:) :: dx_step,dy_step,dz_step,dt_step
integer(kind=4), allocatable, dimension(:) :: nx_step,ny_step,nz_step,nt_step
real(kind=4) :: dx_hist,dy_hist,dz_hist,dt_hist

real(kind=4) :: dxmin,dxmax,dymin,dymax,dzmin,dzmax,dtmin,dtmax
real(kind=4) :: dxmax0,dymax0,dzmax0,dtmax0
integer(kind=4) :: nsrc,neqks,nblast,nshot,nsrc1,nsrc2,isrc,isrc1,isrc2
integer(kind=4) :: n_hist,nx_hist,ny_hist,nz_hist,nt_hist,it,it0

integer(kind=4) :: iunit,flog,i4,i5,i6,i7

iunit=10     ! funit starts at 10
flog=49
open(flog,file='flog.shift_source')
 
write(flog,*) ' enter the number of boxes '
write(*,*) ' enter the number of boxes '
read(*,*) n_hist
write(flog,*) ' n_hist',n_hist

! =================================
!  reading sources from the first set
! =================================
!
! REMINBER    neqks then nshots then nblasts
!

write(flog,*) ' working on sources the first set     ----------------'
open(iunit,file='fsrc.one',access='direct',recl=8*4)  ! eight datas of 4 bytes 
read(iunit,rec=1) nsrc,neqks,nshot,nblast,i4,i5,i6,i7 ! reading numbers of events

allocate(sr1_x(nsrc))      ! x coordinate of source general frame
allocate(sr1_y(nsrc))      ! y
allocate(sr1_z(nsrc))      ! z
allocate(sr1_to(nsrc))     ! t0 origin time
allocate(ki_pos(nsrc))     ! flag on position (if # 0, we invert position)
allocate(ki_to(nsrc))      ! flag on origin time (if #0, we invert origin time)
allocate(ki1_m(nsrc))       ! flag on model   1 inside 0 outside the box
allocate(sr1_id(nsrc))     ! id of the source

call read_fsrc(iunit,sr1_x,sr1_y,sr1_z,sr1_to,ki_pos,ki_to,ki1_m,sr1_id,nsrc)
close(iunit)

deallocate(ki_pos);deallocate(ki_to)   ! free memory

nsrc1=neqks    ! all earthquakes    even when out  (ki1_m will trace that out)

! =================================
!  reading sources from the first set
! =================================

write(flog,*) ' working on sources the first set     ----------------'
open(iunit,file='fsrc.two',access='direct',recl=8*4)  ! eight datas of 4 bytes 
read(iunit,rec=1) nsrc,neqks,nshot,nblast,i4,i5,i6,i7 ! reading numbers of events

allocate(sr2_x(nsrc))      ! x coordinate of source general frame
allocate(sr2_y(nsrc))      ! y
allocate(sr2_z(nsrc))      ! z
allocate(sr2_to(nsrc))     ! t0 origin time
allocate(ki_pos(nsrc))     ! flag on position (if # 0, we invert position)
allocate(ki_to(nsrc))      ! flag on origin time (if #0, we invert origin time)
allocate(ki2_m(nsrc))       ! flag on model   1 inside 0 outside the box
allocate(sr2_id(nsrc))     ! id of the source

call read_fsrc(iunit,sr2_x,sr2_y,sr2_z,sr2_to,ki_pos,ki_to,ki2_m,sr2_id,nsrc)
close(iunit)

deallocate(ki_pos);deallocate(ki_to)   ! free memory
nsrc2=neqks  ! all earthquakes    even when out  (ki1_m will trace that out)

!=====================================
!  end of reading of the two sets
!=====================================
write(*,*) ' total of events: quakes  in and out '
write(*,*) ' number of quakes in the first set ',nsrc1
write(*,*) ' number of quakes in the second set',nsrc2

allocate(sr2_twin(nsrc2))
sr2_twin(:)=999               ! tracking untwin events on the second set

!=====================================
!  shifts in space and in time
!=====================================
isrc=0
dxmin=1.e+29; dxmax=-1.e+29
dymin=1.e+29; dymax=-1.e+29
dzmin=1.e+29; dzmax=-1.e+29
dtmin=1.e+29; dtmax=-1.e+29

allocate(sr_dx(nsrc1)); allocate(sr_dy(nsrc1)); allocate(sr_dz(nsrc1)); allocate(sr_dt(nsrc1))

do isrc1=1,nsrc1
if(ki1_m(isrc1) /= 0) then   ! do not consider events out of the box
  do isrc2=1,nsrc2
    if(sr1_id(isrc1) == sr2_id(isrc2)) then
      sr2_twin(isrc2)=1
      goto 1000
    endif
  enddo ! isrc2
  write(*,*) ' event id first set not in the second set ',sr1_id(isrc1)
  goto 2000
1000 continue   ! we have a twin event : compute differences
  isrc=isrc+1
!============= x
  sr_dx(isrc)=sr1_x (isrc)-sr2_x (isrc)
  if(sr_dx(isrc) < dxmin) dxmin=sr_dx(isrc)
  if(sr_dx(isrc) > dxmax) dxmax=sr_dx(isrc)
!============= y
  sr_dy(isrc)=sr1_y (isrc)-sr2_y (isrc)
  if(sr_dy(isrc) < dymin) dymin=sr_dy(isrc)
  if(sr_dy(isrc) > dymax) dymax=sr_dy(isrc)
!============= z  depth
  sr_dz(isrc)=sr1_z (isrc)-sr2_z (isrc)
  if(sr_dz(isrc) < dzmin) dzmin=sr_dz(isrc)
  if(sr_dz(isrc) > dzmax) dzmax=sr_dz(isrc)
!============= origin time
  sr_dt(isrc)=sr1_to(isrc)-sr2_to(isrc)
  if(sr_dt(isrc) < dtmin) dtmin=sr_dt(isrc)
  if(sr_dt(isrc) > dtmax) dtmax=sr_dt(isrc)
2000 continue
endif
enddo ! isrc1

do isrc2=1,nsrc2
if(ki2_m(isrc2) /= 0) then     ! do not consider event out of the box
  if(sr2_twin(isrc2) /=1) write(*,*) 'event id second set not in the first set',sr2_id(isrc2)
endif
enddo

deallocate(sr1_x);deallocate(sr1_y);deallocate(sr1_z);deallocate(sr1_to);deallocate(sr1_id)
deallocate(sr2_x);deallocate(sr2_y);deallocate(sr2_z);deallocate(sr2_to);deallocate(sr2_id)
deallocate(sr2_twin);deallocate(ki1_m);deallocate(ki2_m)

nsrc=isrc
write(*,*) ' total of twin events ',nsrc
write(*,*) ' compute statistics on these events '

dxmax0=dxmax
if(abs(dxmin) > dxmax0) dxmax0=abs(dxmin)
dymax0=dymax
if(abs(dymin) > dymax0) dymax0=abs(dymin)
dzmax0=dzmax
if(abs(dzmin) > dzmax0) dzmax0=abs(dzmin)
dtmax0=dtmax
if(abs(dtmin) > dtmax0) dtmax0=abs(dtmin)
!=======================================
!  compute structure of histograms
!=======================================
!************************************************************************ in x
nx_hist=n_hist ! resolution histogram  should be odd
dx_hist=dxmax0/(nx_hist/2)
allocate(dx_step(nx_hist)); allocate(nx_step(nx_hist))
dx_step(1)=-dxmax0
do it=2,nx_hist
  dx_step(it)=dx_step(it-1)+dx_hist
enddo ! it
write(*,*) ' dxmax, dx_step_max ',dxmax0,dx_step(nx_hist),dx_hist
!************************************************************************ in y
ny_hist=n_hist ! resolution histogram  should be odd
dy_hist=dymax0/(ny_hist/2)
allocate(dy_step(ny_hist)); allocate(ny_step(ny_hist))
dy_step(1)=-dymax0
do it=2,ny_hist
dy_step(it)=dy_step(it-1)+dy_hist
enddo 
write(*,*) ' dymax, dy_step_max ',dymax0,dy_step(ny_hist),dy_hist
!************************************************************************ in z
nz_hist=n_hist ! resolution histogram  should be odd
dz_hist=dzmax0/(nz_hist/2)
allocate(dz_step(nz_hist)); allocate(nz_step(nz_hist))
dz_step(1)=-dzmax0
do it=2,nz_hist
dz_step(it)=dz_step(it-1)+dz_hist
enddo 
write(*,*) ' dzmax, dz_step_max ',dzmax0,dz_step(nz_hist),dz_hist
!************************************************************************ in to
nt_hist=n_hist ! resolution histogram  should be odd
dt_hist=dtmax0/(nt_hist/2)
allocate(dt_step(nt_hist)); allocate(nt_step(nt_hist))
dt_step(1)=-dtmax0
do it=2,nt_hist
dt_step(it)=dt_step(it-1)+dt_hist
enddo 
write(*,*) ' dtmax, dt_step_max ',dtmax0,dt_step(nt_hist),dt_hist

nx_step(:)=0; ny_step(:)=0; nz_step(:)=0; nt_step(:)=0
!==================
do isrc=1,nsrc
!================== do histogram for dx
it0=0
!if(sr_dx(isrc) >= -dxmax0 .and. sr_dx(isrc) <= dxmax0) then
if(abs(sr_dx(isrc)) <= dxmax0) then
  it0=1
  do it=2,nx_hist
    if(sr_dx(isrc) > dx_step(it)) then
      it0=it
    endif
  enddo   ! it  
  nx_step(it0)=nx_step(it0)+1
endif ! inside the histogram in dx
!================== do histogram for dy
it0=0
!if(sr_dy(isrc) >= -dymax0 .and. sr_dy(isrc) <= dymax0) then
if(abs(sr_dy(isrc)) <= dymax0) then
  it0=1
  do it=2,ny_hist
    if(sr_dy(isrc) > dy_step(it)) then
      it0=it
    endif
  enddo   ! it  
  ny_step(it0)=ny_step(it0)+1
endif ! inside the histogram in dy
!================== do histogram for dz
it0=0
!if(sr_dz(isrc) >= -dzmax0 .and. sr_dz(isrc) <= dzmax0) then
if(abs(sr_dz(isrc)) <= dzmax0) then
  it0=1
  do it=2,nz_hist
    if(sr_dz(isrc) > dz_step(it)) then
      it0=it
    endif
  enddo   ! it  
  nz_step(it0)=nz_step(it0)+1
endif ! inside the histogram in dz
!================== do histogram for dt
it0=0
!if(sr_dt(isrc) >= -dtmax0 .and. sr_dt(isrc) <= dtmax0) then
if(abs(sr_dt(isrc)) <= dtmax0) then
  it0=1
  do it=2,nt_hist
    if(sr_dt(isrc) > dt_step(it)) then
      it0=it
    endif
  enddo   ! it  
  nt_step(it0)=nt_step(it0)+1
endif ! inside the histogram in dt
!==========================
enddo ! nsrc   twin events
!==========================
! open the shell file for gnuplot
open(9,file='SRC/run_source.sh',access='append')

open(iunit,file='SRC/histo_dx')
do it=1,nx_hist
write(iunit,*) dx_step(it)+0.5*dx_hist,nx_step(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "histo_dx.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"histo_dx"'//' linecolor 1 with boxes'
write(9,'(a)') 'EOD'
open(iunit,file='SRC/histo_dy')
do it=1,ny_hist
write(iunit,*) dy_step(it)+0.5*dy_hist,ny_step(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "histo_dy.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"histo_dy"'//' linecolor 1 with boxes'
write(9,'(a)') 'EOD'
open(iunit,file='SRC/histo_dz')
do it=1,nz_hist
write(iunit,*) dz_step(it)+0.5*dz_hist,nz_step(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "histo_dz.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"histo_dz"'//' linecolor 1 with boxes'
write(9,'(a)') 'EOD'
open(iunit,file='SRC/histo_dt')
do it=1,nt_hist
write(iunit,*) dt_step(it)+0.5*dt_hist,nt_step(it)
enddo
close(iunit)
write(9,'(a)') 'gnuplot <<EOD'
write(9,'(a)') 'set term jpeg'
write(9,'(a)') 'set output "histo_dt.jpg"'
write(9,'(a)') 'set style fill solid border -1'
write(9,'(a)') 'plot '//'"histo_dt"'//' linecolor 1 with boxes'
write(9,'(a)') 'EOD'

close(9)   ! shell file

deallocate(sr_dx); deallocate(sr_dy); deallocate(sr_dz); deallocate(sr_dt)
deallocate(dx_step); deallocate(nx_step)
deallocate(dy_step); deallocate(ny_step)
deallocate(dz_step); deallocate(nz_step)
deallocate(dt_step); deallocate(nt_step)

close(flog) ! log file
stop
end program source_analysis

