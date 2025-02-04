include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : fsrc2a fsrc2in a2fsrc a2fsrc_out fsta2a a2fsta fobs2a fobs2tt a2fobs a2fobs_clean fcal2a a2model h2model setfwei proxy_VPVS  model_stat


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

############   fsrc2a
fsrc2a :  fsrc2a.o
	$(F90) $(F90_OPTION) fsrc2a.o -o $(EXEC)/fsrc2a

############   fsrc2a
fsrc2in :  fsrc2in.o
	$(F90) $(F90_OPTION) fsrc2in.o -o $(EXEC)/fsrc2in

############   a2fsrc
a2fsrc :  a2fsrc.o
	$(F90) $(F90_OPTION) a2fsrc.o -o $(EXEC)/a2fsrc

############   a2fsrc_out
a2fsrc_out :  a2fsrc_out.o
	$(F90) $(F90_OPTION) a2fsrc_out.o -o $(EXEC)/a2fsrc_out

############   fsta2a
fsta2a :  fsta2a.o
	$(F90) $(F90_OPTION) fsta2a.o -o $(EXEC)/fsta2a

############   a2fsta
a2fsta :  a2fsta.o
	$(F90) $(F90_OPTION) a2fsta.o -o $(EXEC)/a2fsta

############   fcal2a
fcal2a :  fcal2a.o
	$(F90) $(F90_OPTION) fcal2a.o -o $(EXEC)/fcal2a

############   fobs2a
fobs2a :  fobs2a.o
	$(F90) $(F90_OPTION) fobs2a.o -o $(EXEC)/fobs2a

############   fobs2tt
fobs2tt :  fobs2tt.o
	$(F90) $(F90_OPTION) fobs2tt.o -o $(EXEC)/fobs2tt

############   a2fobs
a2fobs :  a2fobs.o
	$(F90) $(F90_OPTION) a2fobs.o -o $(EXEC)/a2fobs

############   a2fobs_clean
a2fobs_clean :  a2fobs_clean.o
	$(F90) $(F90_OPTION) a2fobs_clean.o -o $(EXEC)/a2fobs_clean

############   a2model
a2model :  a2model.o
	$(F90) $(F90_OPTION) a2model.o -o $(EXEC)/a2model

############   a2model
h2model :  h2model.o
	$(F90) $(F90_OPTION) h2model.o -o $(EXEC)/h2model

############   proxy_VPVS
proxy_VPVS : proxy_VPVS.o
	$(F90) $(F90_OPTION) proxy_VPVS.o -o $(EXEC)/proxy_VPVS

############   setfwei
setfwei :  setfwei.o
	$(F90) $(F90_OPTION) setfwei.o -o $(EXEC)/setfwei

############   model_stat
model_stat :  model_stat.o
	$(F90) $(F90_OPTION) model_stat.o -o $(EXEC)/model_stat
