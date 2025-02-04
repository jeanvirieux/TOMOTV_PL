MODULE s_matrix
contains

!
!
!###################################################################################
!
!            routine permettant de lire une matrice sparse
!
!      ENTREES : iunit : free file unit (output the next free file unit) 
!                flog : log file
!                 
!      SORTIES : 
!                nnz
!                nl
!                mat
!                ic
!                id
!
!
subroutine read_matrix(fmat,fic,fid,mat,ic,id,nnz,nl,nc,iunit,flog)
implicit none
integer(kind=4) :: nnz,nl,nc,iunit,flog
integer(kind=4) :: umat,uic,uid
character(len=*) :: fmat,fic,fid
real(kind=4), dimension(:) :: mat(:)
integer(kind=4), dimension(:) :: ic(:),id(:)

! logical units
umat=iunit
iunit=iunit+1
uic=iunit
iunit=iunit+1
uid=iunit
iunit=iunit+1

write(flog,'(20x,a35,5x,a15)') 'elements non nuls :',fmat
write(flog,'(20x,a35,5x,a15)') 'index des colonnes :',fic
write(flog,'(20x,a35,5x,a15)') 'pointeur de lignes :',fid
write(flog,'(20x,a35,5x,i12)') 'nombre d''elements non nuls :',nnz
write(flog,'(20x,a35,5x,i12)') 'nombre de lignes :',nl
write(flog,'(20x,a35,5x,i12)') 'nombre de colonnes :',nc

!
! opening files
!
write(flog,*) '---------------------------------------------------------------' 
write(flog,'(35x,a)') ' FICHIERS BINAIRES' 

write(flog,'(20x,a16,a10)')'-> open/read of ',fmat
open(umat,file=fmat,access='direct',recl=4*nnz)    ! file where are values of the matrix X
read(umat,rec=1) mat
close(umat)

write(flog,'(20x,a16,a10)')'-> open/read of ',fic
open(uic,file=fic,access='direct',recl=4*nnz)      ! file where are row index where are non-zero elements
read(uic,rec=1) ic
close(uic)

write(flog,'(20x,a16,a10)')'-> open/read of ',fid
open(uid,file=fid,access='direct',recl=4*nl+4)   ! file where are the column index where non-zero elements
read(uid,rec=1) id   
close(uid)

write(flog,'(35x,a)') ' SUCCESS IN READING INVERSION BINARY FILES'

return
end subroutine read_matrix

END MODULE s_matrix
