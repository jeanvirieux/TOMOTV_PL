include ../../make.inc
#
EXEC=../../BIN_LINUX_PL
LIBS=
INC=

all : hypo2fobsfsrc hypo2fobsfsrcfwei hypo2fsta


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

SUB = utmsubs.o sub_utime.o

############   hypo2fobsfsrc
hypo2fobsfsrc :  hypo2fobsfsrc.o $(SUB)
	$(F90) $(F90_OPTION) hypo2fobsfsrc.o $(SUB) -o $(EXEC)/hypo2fobsfsrc

############   hypo2fobsfsrcfwei
hypo2fobsfsrcfwei :  hypo2fobsfsrcfwei.o $(SUB)
	$(F90) $(F90_OPTION) hypo2fobsfsrcfwei.o $(SUB) -o $(EXEC)/hypo2fobsfsrcfwei

############   hypo2fsta
hypo2fsta :  hypo2fsta.o $(SUB)
	$(F90) $(F90_OPTION) hypo2fsta.o $(SUB) -o $(EXEC)/hypo2fsta
