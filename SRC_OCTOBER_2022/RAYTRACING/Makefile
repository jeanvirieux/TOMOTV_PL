include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : raytrace.tomo

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

SUB = sub_hs.o submain.o read_par.o read_fsrc.o read_fsta.o read_write_fobs.o sub_inter.o time_3d.o

raytrace.tomo : $(SUB) raytrace.tomo.o
	$(F90) $(F90_OPTION) raytrace.tomo.o $(SUB) $(LIBS) -o $(EXEC)/raytrace.tomo

