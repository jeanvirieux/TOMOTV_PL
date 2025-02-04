include	../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : detection_dd

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

detection_dd : detection_dd.o
	$(F90) $(F90_OPTION) detection_dd.o $(LIBS) -o $(EXEC)/detection_dd

