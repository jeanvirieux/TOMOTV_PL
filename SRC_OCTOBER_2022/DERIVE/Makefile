include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : derive_SLOWNESS.tomo derive_SLOWNESS.tomo_20s

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
SUB = sub_inter.o read_fobs.o read_fray.o read_fsrc.o read_par.o \
      store_line.o submain.o valu_lent.o block_data.o slowness2angle.o

SUB1 = read_fobs.o read_fsrc.o read_matrix.o store_ddf_line.o 

#

derive_SLOWNESS.tomo : $(SUB) derive_SLOWNESS.tomo.o
	$(F90) $(F90_OPTION) derive_SLOWNESS.tomo.o $(SUB) $(LIBS) -o $(EXEC)/derive.tomo

derive_SLOWNESS.tomo_20s : $(SUB) derive_SLOWNESS.tomo_20s.o
	$(F90) $(F90_OPTION) derive_SLOWNESS.tomo_20s.o $(SUB) $(LIBS) -o $(EXEC)/derive.tomo_20s
