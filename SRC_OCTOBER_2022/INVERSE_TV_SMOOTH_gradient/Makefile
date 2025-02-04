include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : inverse_TV_SMOOTH.tomo 

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(F90_OPTION) $(INC) -c $*.f90

.SUFFIXES: .f .o
.f.o:
	$(F77) $(F77_OPTION) $(INC) -c $*.f

.SUFFIXES: .c .o
.c.o:
	$(CC) $(C_OPTION) $(INC) -c $*.c

clean:
	$(RM) *.o *~ *.mod \#*

#

SUB_3D = aprod.o lsqrblas.o lsqr.o \
      new_parameters.o \
      read_matrix.o sub_read.o sub_smoothing.o sub_nesterov_3D.o


inverse_TV_SMOOTH.tomo : $(SUB_3D) inverse.tomo.o
	$(F90) $(F90_OPTION) inverse.tomo.o $(SUB_3D) $(LIBS) -o $(EXEC)/inverse_TV_SMOOTH.tomo

