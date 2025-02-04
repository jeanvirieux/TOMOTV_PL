include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : facqui2tv model2power2 power2_sampling test_fmg model3D_setting power2_finer model2D_strip

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
SUB = sub_interpol.o sub_prolong.o
#

facqui2tv : facqui2tv.o
	$(F90) $(F90_OPTION) facqui2tv.o $(LIBS) -o $(EXEC)/facqui2tv

power2_sampling : power2_sampling.o $(SUB)
	$(F90) $(F90_OPTION) power2_sampling.o $(SUB) $(LIBS) -o $(EXEC)/power2_sampling

power2_finer : power2_finer.o $(SUB)
	$(F90) $(F90_OPTION) power2_finer.o $(SUB) $(LIBS) -o $(EXEC)/power2_finer

model2power2 : $(SUB) model2power2.o 
	$(F90) $(F90_OPTION) model2power2.o $(SUB) $(LIBS) -o $(EXEC)/model2power2

test_fmg : test_fmg.o $(SUB)
	$(F90) $(F90_OPTION) test_fmg.o $(SUB) $(LIBS) -o $(EXEC)/test_fmg

model3D_setting : model3D_setting.o 
	$(F90) $(F90_OPTION) model3D_setting.o $(LIBS) -o $(EXEC)/model3D_setting

model2D_strip : model2D_strip.o 
	$(F90) $(F90_OPTION) model2D_strip.o $(LIBS) -o $(EXEC)/model2D_strip
