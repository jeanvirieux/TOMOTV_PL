include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : station_res residue_gnuplot shift_gnuplot

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
SUB = read_fobs.o read_fsrc.o 

#

station_res : $(SUB) station_res.o
	$(F90) $(F90_OPTION) station_res.o $(SUB) $(LIBS) -o $(EXEC)/station_res

residue_gnuplot : $(SUB) residue_gnuplot.o
	$(F90) $(F90_OPTION) residue_gnuplot.o $(SUB) $(LIBS) -o $(EXEC)/residue_gnuplot

shift_gnuplot : $(SUB) shift_gnuplot.o
	$(F90) $(F90_OPTION) shift_gnuplot.o $(SUB) $(LIBS) -o $(EXEC)/shift_gnuplot

