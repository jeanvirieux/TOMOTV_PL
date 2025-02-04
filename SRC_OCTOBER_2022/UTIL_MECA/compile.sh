# subroutine de conversion angles ... validation of slowness2angle.f90 (not used here ... integrated into tomoTV)
gfortran -c slowness2angle.f90
# conversion from shell format UTM into local XYZ coordinates in meters
gfortran -o ../bin/sta2fsta sta2fsta.f90
# conversion from NLLOC location output into fsrc in meters and seconds
gfortran -c sub_utime.f
gfortran -c pick_line_space.f90 
gfortran -o ../bin/src2fsrc src2fsrc.f90 sub_utime.o pick_line_space.o
# conversion from utime to date (to test file fdate.txt) not used for tomoTV
gfortran -o ../bin/utime2date utime2date.f90
# analysis of the output of NLLOC  ... linking event and picking file 
gfortran -o ../bin/pick_event pick_event.f90  pick_line_space.o
# construction of fobs and fsph
gfortran -c pick_line_tab.f90 
gfortran -o ../bin/pick2fobs-fsphe pick2fobs-fsphe.f90 sub_utime.f pick_line_space.o pick_line_tab.o
# fusion des infos pointage et rayon dans focal
gfortran -o../bin/fsphe_fmeca2focal fsphe_fmeca2focal.f90
# conversion vers ascii
gfortran -o ../bin/fmeca2a fmeca2a.f90
# conversion vers ascii
gfortran -o ../bin/fsphe2a fsphe2a.f90
# conversion vers ascii
gfortran -o ../bin/focal2a focal2a.f90
# conversion vers bin
gfortran -o ../bin/a2fmeca a2fmeca.f90
# conversion vers bin
gfortran -o ../bin/a2fsphe a2fsphe.f90
# conversion vers bin
gfortran -o ../bin/a2focal a2focal.f90

######################################
#
# small programs to test slowness2angle and to compare with Isabelle matlab code
#
######################################
#gfortran -o ../bin/test_angles test_angles.f90 slowness2angle.o
# comparison from Isabelle output and current output. not used for tomoTV
#gfortran -o ../bin/compare compare.f90 slowness2angle.o
rm *.o

