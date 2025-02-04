MODULE s_write_matrix
contains

!-------------------------------------------------------------------
!
!     ecriture des fichiers constituant une matrice sparse
!
!-------------------------------------------------------------------

subroutine write_matrix(mat,ic,id,nnz,nl,umat,uic,uid,fmat,fic,fid,flog)
implicit none

integer(kind=4) :: umat,uic,uid,flog
integer(kind=4) :: nnz,nl
integer(kind=4), dimension(:) :: ic(:),id(:)   !ic(nnz),id(nl+1)
real(kind=4), dimension(:) :: mat(:)           ! mat(nnz)

character(len=*) fmat,fic,fid

write(flog,'(35x,a)') ' ECRITURE DES FICHIERS BINAIRES PONDERES'
write(flog,'(20x,a)') '-> ecriture de la matrice'
open(umat,file=fmat,access='direct',recl=4*nnz)
write(umat,rec=1) mat
close(umat)

write(flog,'(20x,a)') '-> ecriture de l''index des colonnes'
open(uic,file=fic,access='direct',recl=4*nnz)
write(uic,rec=1) ic
close(uic)

write(flog,'(20x,a)') '-> ecriture du pointeur de ligne'
open(uid,file=fid,access='direct',recl=4*nl+4)
write(uid,rec=1) id
close(uid)

write(flog,'(35x,a)') 'FIN D''ECRITURE'

return
end subroutine write_matrix

!===========
!   on ecrit moins dans ce cas
!===========
subroutine write_matrix_check(mat,ic,id,nnz,nl,umat,uic,uid,fmat,fic,fid,flog)
implicit none

integer(kind=4) :: umat,uic,uid,flog
integer(kind=4) :: nnz,nl,i
integer(kind=4), dimension(:) :: ic(:),id(:)   !ic(nnz),id(nl+1)
real(kind=4), dimension(:) :: mat(:)           ! mat(nnz)

character(len=*) fmat,fic,fid

write(flog,'(35x,a)') ' ECRITURE DES FICHIERS BINAIRES PONDERES'
write(flog,'(20x,a)') '-> ecriture de la matrice'
open(umat,file=fmat,access='direct',recl=4*nnz)
write(umat,rec=1) (mat(i),i=1,nnz)
close(umat)

write(flog,'(20x,a)') '-> ecriture de l''index des colonnes'
open(uic,file=fic,access='direct',recl=4*nnz)
write(uic,rec=1) (ic(i),i=1,nnz)
close(uic)

write(flog,'(20x,a)') '-> ecriture du pointeur de ligne'
open(uid,file=fid,access='direct',recl=4*nl+4)
write(uid,rec=1) (id(i),i=1,nl+1)
close(uid)

write(flog,'(35x,a)') 'FIN D''ECRITURE'

return
end subroutine write_matrix_check
END MODULE s_write_matrix

