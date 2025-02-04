####################### compilateur gnu ################ DEGUG - SLOW ##############
F77=gfortran
F90=gfortran
CC=gcc
F90_OPTION=-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicitnone -ftraceback -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
F90_OPTION=-cpp -fno-second-underscore -ffree-line-length-255 -D -g -fbounds-check -Wall -fbacktrace -finit-real=nan -fopenmp
F90_OPTION=-cpp	-O3 -fopenmp -Wunused -Wuninitialized -fcheck=all -fopenmp
F77_OPTION=-O3 
C_OPTION= -O3 -fopenmp
CC_OPTION= -O3 -fopenmp

# options pour le gnu g77 :
# -fno-second-underscore : elimine le second underscore lors de l'interfacage avec le c. 
# -w : elimine les warning
# -Wunused : cherche les variables inutilisees
# -Wuninitialized : regarde les variables non initialisees (s'utilise avec -O)
# -Wno-globals : enleve les warning lors des passages du buffer en parametre.
#
# OPTIMISATION 
#
# -O -O1 -O2 -O3 : niveau d'optimisation
#
# option Albanne sur intel
# -warn all -check all -O0 -fno-alias -fno-fnalias -g -debug -traceback -implicitnone -fpe0 -fp-stack-check -ftrapuv -heap-arrays -gen-interface -warn interface
#

