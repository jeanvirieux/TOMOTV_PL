MODULE s_read_matrix
contains

!
!
!###################################################################################
!
!            routine permettant de lire une matrice sparse
!
!      ENTREES : umat,uic,uid : file units
!                flog : log file
!                nnz
!                nl
!                 
!      SORTIES : 
!                mat    : non-zero elements of sparse matrix A
!                ic     : column index
!                id     : row index
!===================================================================================
!  -fmat contient les elements non nuls de la matrice A. 
!  -fic contient les numeros de colonne de chaque elements de fmat:
!                l'elemement mat(j) est dans la colonne ic(j) de A.
!  -fid est un pointeur servant a extraire une ligne i de la matrice A:
!                les elements non nuls de la ligne i sont:
!
!       [mat(id(i)),mat(id(i)+1),...,mat(id(i)+k),...,mat(id(i+1)-1)]
!
!  Ces trois fichiers suffisent pour reconstituer la matrice:
!                l'element mat(id(i)+k) se situe sur la ligne i de A et
!                se situe sur la colonne ic(id(i)+k) de A.
!                    
!                A( i, ic(id(i)+k) ) = mat(id(i)+k).
!
!----------------------------------------------------------------------------------
!
subroutine read_matrix(mat,ic,id,nnz,nl,nc,umat,uic,uid,fmat,fic,fid,flog)
implicit none
integer(kind=4) :: nnz,nl,nc,flog
integer(kind=4) :: umat,uic,uid
character(len=*) :: fmat,fic,fid
real(kind=4), dimension(:) :: mat(:)
integer(kind=4), dimension(:) :: ic(:),id(:)

write(flog,'(20x,a35,5x,a15)') 'elements non nuls :',fmat
write(flog,'(20x,a35,5x,a15)') 'index des colonnes :',fic
write(flog,'(20x,a35,5x,a15)') 'pointeur de lignes :',fid
write(flog,'(20x,a35,5x,i10)') 'nombre d''elements non nuls :',nnz
write(flog,'(20x,a35,5x,i10)') 'nombre de lignes :',nl
write(flog,'(20x,a35,5x,i10)') 'nombre de colonnes :',nc

!
! opening files
!
write(flog,*) '---------------------------------------------------------------' 
write(flog,'(35x,a)') ' BINARY FILES ' 

write(flog,'(20x,a16,a10)')'-> open/read of ',fmat
open(umat,file=fmat,access='direct',recl=4*nnz)    ! file : non-zero elements of A
read(umat,rec=1) mat
close(umat)

write(flog,'(20x,a16,a10)')'-> open/read of ',fic
open(uic,file=fic,access='direct',recl=4*nnz)  ! file : row index of non-zero elements
read(uic,rec=1) ic
close(uic)

write(flog,'(20x,a16,a10)')'-> open/read of ',fid  ! longueur nl+1   attention au +1
open(uid,file=fid,access='direct',recl=4*nl+4) ! file : column index of non-zero elements
read(uid,rec=1) id   
close(uid)

write(flog,'(35x,a)') ' SUCCESS IN READING INVERSION BINARY FILES'

return
end subroutine read_matrix

!
! une version pour lire seulement les tailles demandees avec des offsets pour le stockage
!
!
subroutine read_matrix_check(mat,ic,id,nnz,nl,nc,nnz_off,nl_off,umat,uic,uid,fmat,fic,fid,flog)
implicit none
integer(kind=4) :: nnz,nl,nc,nnz_off,nl_off
integer(kind=4) :: flog,i
integer(kind=4) :: umat,uic,uid
character(len=*) :: fmat,fic,fid
real(kind=4), dimension(:) :: mat(:)
integer(kind=4), dimension(:) :: ic(:),id(:)

write(flog,'(20x,a35,5x,a15)') 'elements non nuls :',fmat
write(flog,'(20x,a35,5x,a15)') 'index des colonnes :',fic
write(flog,'(20x,a35,5x,a15)') 'pointeur de lignes :',fid
write(flog,'(20x,a35,5x,i10)') 'nombre d''elements non nuls :',nnz
write(flog,'(20x,a35,5x,i10)') 'nombre de lignes :',nl
write(flog,'(20x,a35,5x,i10)') 'nombre de colonnes :',nc

!
! opening files
!
write(flog,*) '---------------------------------------------------------------' 
write(flog,'(35x,a)') ' FICHIERS BINAIRES' 

write(flog,'(20x,a16,a10)')'-> open/read of ',fmat
open(umat,file=fmat,access='direct',recl=4*nnz)    ! file : non-zero elements of A
read(umat,rec=1) (mat(i),i=nnz_off+1,nnz_off+nnz)
close(umat)

write(flog,'(20x,a16,a10)')'-> open/read of ',fic
open(uic,file=fic,access='direct',recl=4*nnz)  ! file : row index of non-zero elements
read(uic,rec=1) (ic(i),i=nnz_off+1,nnz_off+nnz)
close(uic)

write(flog,'(20x,a16,a10)')'-> open/read of ',fid  ! longueur nl+1   attention au +1
open(uid,file=fid,access='direct',recl=4*nl+4) ! file : column index of non-zero elements
read(uid,rec=1) (id(i),i=nl_off+1,nl_off+nl+1)   
close(uid)

write(flog,'(35x,a)') ' SUCCESS IN READING INVERSION BINARY FILES: CHECK OPTION'

return
end subroutine read_matrix_check

END MODULE s_read_matrix
