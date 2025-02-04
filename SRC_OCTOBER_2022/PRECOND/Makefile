include	../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : precond.tomo

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
SUB = read_fobs.o read_input.o read_matrix.o sub_lap.o sub_precond.o write_matrix.o write_rms.o

#

precond.tomo : $(SUB) precond.tomo.o
	$(F90) $(F90_OPTION) precond.tomo.o $(SUB) $(LIBS) -o $(EXEC)/precond.tomo

