gfortran -c read_fobs.f90
gfortran -c read_input.f90
gfortran -c read_matrix.f90
gfortran -c sub_lap.f90
gfortran -c sub_precond.f90
gfortran -c write_matrix.f90
gfortran -c write_rms.f90
gfortran -static -o precond.exe precond.tomo.f90 read_fobs.o read_input.o read_matrix.o sub_lap.o sub_precond.o write_matrix.o write_rms.o
del ..\..\BIN\precond.exe

copy precond.exe ..\..\BIN\precond.exe

del *.o

del *.mod
del *.exe

