include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : smoothgauss3dfast smooth2d

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
SUB = 
#

smoothgauss3dfast : $(SUB) smoothgauss3dfast.o
	$(F90) $(F90_OPTION) smoothgauss3dfast.o $(SUB) $(LIBS) -o $(EXEC)/smoothgauss3dfast

smooth2d : $(SUB) smooth2d.o
	$(F90) $(F90_OPTION) smooth2d.o $(SUB) $(LIBS) -o $(EXEC)/smooth2d
