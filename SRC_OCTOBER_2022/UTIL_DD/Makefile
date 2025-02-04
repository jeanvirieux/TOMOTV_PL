include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=
SUB= sub_utime.o utmsubs.o read_fobs.o 

all : fsrc2dd fsta2dd


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

############   fsrc2dd
fsrc2dd :  fsrc2dd.o
	$(F90) $(F90_OPTION) fsrc2dd.o -o $(EXEC)/fsrc2dd

############   fsta2dd
fsta2dd :  fsta2dd.o
	$(F90) $(F90_OPTION) fsta2dd.o -o $(EXEC)/fsta2dd
