include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : inverse.tomo 

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
#SUB = aprod.o lsqr_sub.o lsqr_m.o 

SUB = aprod.o lsqrblas.o lsqr.o \
      new_parameters.o \
      read_matrix.o sub_read.o
#

inverse.tomo : $(SUB) inverse.tomo.o
	$(F90) $(F90_OPTION) inverse.tomo.o $(SUB) $(LIBS) -o $(EXEC)/inverse.tomo

