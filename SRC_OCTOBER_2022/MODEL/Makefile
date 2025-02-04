include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : model.tomo model.tomo_old model.nll model.interpol model.coupe model.coupebis model.cross modelP2_5D modelS2_5D model_synt model_layer model_grid model_smooth model_spike model_spike_pattern model_spike_iter model_diff model_diff_sum model2Dto3D model2Dto3D_dummy  model2Dto3D_prior_x

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
SUB = sub_interpol.o 
#

model.tomo : $(SUB) model.tomo.o
	$(F90) $(F90_OPTION) model.tomo.o $(SUB) $(LIBS) -o $(EXEC)/model.tomo

model.tomo_old : $(SUB) model.tomo_old.o
	$(F90) $(F90_OPTION) model.tomo_old.o $(SUB) $(LIBS) -o $(EXEC)/model.tomo_old

model.nll : $(SUB) model.nll.o
	$(F90) $(F90_OPTION) model.nll.o $(SUB) $(LIBS) -o $(EXEC)/model.nll

model.interpol : $(SUB) model.interpol.o
	$(F90) $(F90_OPTION) model.interpol.o $(SUB) $(LIBS) -o $(EXEC)/model.interpol

model.coupe : $(SUB) model.coupe.o
	$(F90) $(F90_OPTION) model.coupe.o $(SUB) $(LIBS) -o $(EXEC)/model.coupe

model.coupebis : $(SUB) model.coupebis.o
	$(F90) $(F90_OPTION) model.coupebis.o $(SUB) $(LIBS) -o $(EXEC)/model.coupebis

model.cross : $(SUB) model.cross.o
	$(F90) $(F90_OPTION) model.cross.o $(SUB) $(LIBS) -o $(EXEC)/model.cross

modelP2_5D : $(SUB) modelP2_5D.o
	$(F90) $(F90_OPTION) modelP2_5D.o $(SUB) $(LIBS) -o $(EXEC)/modelP2_5D

modelS2_5D : $(SUB) modelS2_5D.o
	$(F90) $(F90_OPTION) modelS2_5D.o $(SUB) $(LIBS) -o $(EXEC)/modelS2_5D

model_synt : $(SUB) model_synt.o
	$(F90) $(F90_OPTION) model_synt.o $(SUB) $(LIBS) -o $(EXEC)/model_synt

model_layer : $(SUB) model_layer.o
	$(F90) $(F90_OPTION) model_layer.o $(SUB) $(LIBS) -o $(EXEC)/model_layer

model_grid : $(SUB) model_grid.o
	$(F90) $(F90_OPTION) model_grid.o $(SUB) $(LIBS) -o $(EXEC)/model_grid

model_smooth : model_smooth.o
	$(F90) $(F90_OPTION) model_smooth.o $(LIBS) -o $(EXEC)/model_smooth

model_diff : model_diff.o
	$(F90) $(F90_OPTION) model_diff.o $(LIBS) -o $(EXEC)/model_diff

model_diff_sum : model_diff_sum.o
	$(F90) $(F90_OPTION) model_diff_sum.o $(LIBS) -o $(EXEC)/model_diff_sum

model_spike : model_spike.o
	$(F90) $(F90_OPTION) model_spike.o $(LIBS) -o $(EXEC)/model_spike

model_spike_pattern : model_spike_pattern.o
	$(F90) $(F90_OPTION) model_spike_pattern.o $(LIBS) -o $(EXEC)/model_spike_pattern

model_spike_iter : model_spike_iter.o
	$(F90) $(F90_OPTION) model_spike_iter.o $(LIBS) -o $(EXEC)/model_spike_iter

model2Dto3D : model2Dto3D.o
	$(F90) $(F90_OPTION) model2Dto3D.o -o $(EXEC)/model2Dto3D

model2Dto3D_dummy : model2Dto3D_dummy.o
	$(F90) $(F90_OPTION) model2Dto3D_dummy.o -o $(EXEC)/model2Dto3D_dummy

model2Dto3D_prior_x : model2Dto3D_prior_x.o
	$(F90) $(F90_OPTION) model2Dto3D_prior_x.o -o $(EXEC)/model2Dto3D_prior_x
