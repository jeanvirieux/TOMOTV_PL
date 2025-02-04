# TOMOTV_PL
First Arrival Delayed Time Tomography
tomoTV
By Jean Virieux
Version 0.00
Winter 2015

Contribution of Ortensia Amoroso, Stéphanie Gauthier, Jean-Luc Got, Jenny Jacques, Diana
Latorre, Vadim Monteiller, Stéphane Operto, Céline Ravaut, Tiziana Vanorio,.
- CAUTION A: files are different from the previous TLR3 software (conversion program is
provided as it is for helping you transforming TLR3 input files into tomoTV files).
- CAUTION B: this software has not yet the double difference feature developed in the original
tomoTV code of Monteiller …

# General purpose

We propose a modular hypocenter-velocity tomography computer code written in
FORTRAN90 (only one subroutine is written in C as the eikonal solver of Podvin & Lecomte).
When one wishes locating hypocenters and inverting for velocity structures simultaneously, the
different steps may require various computer resources and one may repeat a new step without
performing previous ones for which information has already been obtained.
For theoretical approaches, the reading of books and articles are recommended. Books we may
recommend are
- D. Zhao, 2015, Multiscale Seismic Tomography, Springer Geophysics, Springer
- K. Aki, 1993, Overview in Seismic tomography: theory and practice, edited by H.M.
Iyer and K. Hirahara, pp 1.-8., Chapman & Hall, London.
- J. Bee Bednar, 1989, Geophysical Inversion, SIAM
- G. Nolet, 1987, Seismic tomography with applications in global seismology and
exploration geophysics, Kluwer Academic Pub., 386 pp.
- G. Nolet, 2008, A breviary of seismic tomography, Imaging the interior of the Earth and
Sun, Cambridge University Press, 339pp.
Many papers are available for seismic travel time tomography: here are those where the
approach we are going to describe has been used
Le Meur, H., J. Virieux and P. Podvin, 1997, Seismic tomography of the Gulf of Corinth: a
comparison of methods, Annali di Geofisica, XL, 1, pp1-25
- Ghose,S., M.W. Hamburger & J. Virieux, 1998. Three-dimensional velocity structure and earthquake
locations beneath the northern Tien Shan of Kyrgyzstan, central Asia, J. Geophys. Res., 103, 2725-2748.
- Latorre, D., Virieux, J., Monfret, T., Monteiller, V., Vanorio, T., Got, J.-L. & H Lyon-Caen, 2004. A
new seismic tomography of Aigion area (Gulf of Corinth, Greece) from the 1991 data set, Geophys. J.
Int., 159, 1013-1031.- Monteiller V., Got J.-L., J. Virieux & P. Okubo, 2005, An efficient algorithm
for double-difference tomography and location in heterogeneous media, with an application to Kilauea
volcano, J. Geophys. Res., 110, B12306, doi:10.1029/2004JB003466.
- Vanorio T., J. Virieux, P. Capuano & G. Russo, 2005, Three-dimensional seismic tomography from
P wave and S wave microearthquake travel times and rock physics characterization of the Campi
Flegrei Caldera, J. Geophys. Res., 110, B03201, doi:10.1029/2004JB003102.
- Gautier S., Latorre D., Virieux J., Deschamps A., Skarpelos C., Sotiriou A., Serpetsidaki A. &
Tselentis A., 2006, A New Passive Tomography of the Aigion Area (Gulf Of Corinth, Greece) from
the 2002 Data set, Pure appl. Geophys, 163, 431-453.

# Root directory architecture

Inside this directory you could install anywhere, you will find the following structure

BIN_LINUX_PL: where should be binaries of computer codes
DOC_TOMOTV: where you will find this document
SRC_OCTOBER_2022: where you will find computer source codes
MAKE.DIR: where you will find typical files you should have in the root directory as make.inc
to be used by various files Makefile.This is the minimal content of the root directory.
You will find as well many directories for building your acquisition configuration and for
running few examples.

Although seismic tomography is based on pretty images, no graphic tools are provided and you
may rely on those you are familiar with. The graphical environment that is often used comes
from CWP project through SEISMIC UNIX (SU) software as well as from GNU. We use
essentially ximage for quick check-up of velocities and GNUPLOT for quick check-up of
hypcenters. Other familiar tools is the GMT tool.

# The directory SRC_OCTOBER_2022 contains :

Makefile: the make file for compiling of the various executable programs
- Directory MODEL: where are source codes for reading the information on the velocity structure
- Directory RAYTRACING: where are source codes for doing two-points ray tracing
- Directory DERIVE: where are source codes for computing synthetic times and
sensitivity matrix
- Directory PRECOND: where are source codes for weighting data and parameters as
well as smoothing techniques
- Directory INVERSE: where are source codes for the inversion and the updating of
parameters, velocity and hypocenter/blast parameters.
- Directory STATISTICS : where are source codes for doing statistics on stats, srcs and
velocities
- Directory UTIL: where are various useful programs
- Directory THIRD_PARTY: where you will find different programs which may not be
distributed (see first lines of each program)

# Compiling and executing

Inside the subdirectory SRC_OCTOBER_2022, one should perform the command make. The argument should be
all if you want to compile the different modules. Therefore, you should type “make all”. The
execution of the Makefile will proceed through dedicated directories of the various modules
and codes should be compiled. We have tested the code using the gfortran 4.8. of GNU
distribution as well as the ifort 14 of INTEL. You may go down inside one module and type
make if you want to recompile only this one.

# Data description

Aside the velocity description, one has to deal with seismic stations information, sources
information and travel-times information. As various dimensions and configurations could be
met, after different trial experiences, we have decided to consider only binary files. This
may be a difficulty but analysing unexpected outputs is also time-consuming.

Four files should be constructed.

# File fsta

Station information: the file is a binary file fsta with equal length of 28 octets.
Utilities are provided a2fsta and fsta2a for performing conversion from ascii to binary and from
binary to ascii.
The first record has only one useful information: the number of stations nsta with six other
dummy 4-octets variables.
The number of other records is, therefore, nsta, leading to a total of (nsta+1) records of 28
octets.
Each record following the first one describes a station with the following information
(1) id_sta: integer(kind=4) a unique id for the station (different between stations)
(2) x_sta: real(kind=4) the x position of the station in the reference frame
(3) y_sta: real(kind=4) the y position of the station in the reference frame
(4) z_sta: real(kind=4) the z position of the station in the reference frame
(5) dtp_sta: real(kind=4) the static delay time for P waves at this station
(6) dts_sta: real(kind=4) the static delay time for S waves at this station
(7) name: character(len=4) the label of the station in four letters (quite common for geographical
identification … not used in computer codes). Having four letters is mandatory (no empty
space).
The easiest way for constructing the fsta file is through a raw text file where each line describes
one record (and, therefore, one station). By using the a2fsta program, one can convert this text
file into the binary file used for the tomography. The reading of this program might be useful
for tuned conversion of existing data file.

# File fsrc

Source information: the file is a binary file fsrc with an equal length of 32 octets
Utilities are provided a2fsrc and fsrc2a for performing conversion from ascii to binary and from
binary to ascii.
The first record has the following information: the number of sources nsrc, the number of
earthquakes neqks, the number of shots nshot, the number of blasts nblast, the number of out-
of-the-box earthquakes neqks_out which may increase during iterations. Once a hypocenter
moves out, there is no attempt to put it back without modifying fsrc (or starting with a new
velocity structure). The total of earthquakes (neqks+neqks_out) stays constant. Three dummy
variables of 4 octets pad the end of the record.
The number of other records is, therefore, nsrc, leading to a total of (nsrc+1) records of 32
octets.
Each record following the first one describes a source with the following information
(1) x_src: real(kind=4) the x position of the source in the reference frame
(2) y_src: real(kind=4) the y position of the source in the reference frame
(3) z_src: real(kind=4) the z position of the source in the reference frame
(4) t0_src: real(kind=4) the origin time of the source
(5) kp_src: integer(kind=4) position flag for quakes (continuous increment for easier
debugging) (0 for others). If the flag is non-zero, we consider this event for position inversion.
(6) kt_src: integer(kind=4) origin time flag for quakes and blasts (continuous increment for
easier debugging) (0 for others). If the flag is non-zero, we consider this event for origin time.
(7) in_src: integer(kind=4) status of the source (in=1 means in the box, out=0 means out). Once
this number is set to zero, only specific editing manipulation will put back the event into the
tomography optimization
(8) id_src: integer(kind=4) unique id of the source (should be unique). Although this
identification number is arbitrary, keep it small will simplify the tracking during the
tomography. An incremental strategy is advisable.
The easiest way for constructing the fsrc file is through a raw text file where each line describes
one record. By using the a2fsrc program, one can convert this text file into the binary file used
for the tomography.
BE AWARE THAT THIS FILE WILL BE MODIFIED THROUGH THE INVERSION
PROCEDURE.

# File fobs

Time observation: the file is a binary file fobs with an equal length of 24 octets
Utilities are provided a2fobs and fobs2a for performing conversion from ascii to binary
(fobs.asc  fobs) and from binary to ascii (fobs  fobs.asc).
The first record has the following information: the total number of picked times nt, the number
of P wave picked times ntp and the number of S wave picked times nts. Three dummy variables
of four octets pad the end of the record.
The number of other records is, therefore, nt, leading to a total of (nt+1) records of 24 octets.
Records up to ntp+1 are related to P wave times and records from ntp+2 to ntp+2+nts related
to S wave times.
Each record following the first one describes a data picked time with the following information
(1) id_obs: integer(kind=4) a unique integer for defining this data point
(2) t_obs: real(kind=4) arrival time of the picked wave (P or S)
(3) dt_obs: real(kind=4) precision in the picking
(4) id_src: integer(kind=4) id of the associated source
(5) id_sta: integer(kind=4) id of the associated station
(6) id_ray: integer(kind=4) id of the associated rayThe last information is built during the inversion and is not required as input: you could set it
to zero. This is a way of tracking ray information while scanning the data in a loop from 1 to nt
(total number of observations).
The easiest way for constructing the fobs file is through a raw text file where each line describes
one record. By using the a2fobs program, one can convert this text file into the binary file used
for the tomography.
This file will be modified through the inversion procedure for the id_ray parameter ONLY.
PLEASE consider that the ID of the data, id_dat, should stay the same during the investigation
with one dataset.
Computer codes for the conversion of ascii files into binary files might be useful for tuning data
in another format. They are in the directory UTIL under the directory SRC where you will find
programs for conversion.

# Definition of the velocity models

%File modelP
The P wave velocity structure, as well as the S wave velocity structure, is described in a binary
file (direct-access file using the terminology of FORTRAN) with the following order of storage:
x is the fastest, y the second and z the third. Only float values of the velocity are stored.

%File modelS
The S wave velocity is described as well following the same grid description.

% File model.head
Information regarding the velocity model is provided inside the file model.head: we have
always one line of comments and then input values in the second line and so on.
This rectangular grid is attached to a reference frame with an absolute origin (xinv,yinv,zinv)
related to the tomographic box. Numbers of points along the three directions, nxinv,nyinv,nzinv
are defined as well as the grid stepping in the three direction hxinv,hyinv,hzinv. These
quantities are provided through the input file “modelP.head” when scanned by the first
computer code “model.tomo” which is dedicated to the construction of the fine eikonal cubic
grid used by the eikonal solver (Podvin & Lecomte, 1991).
We assume that the origin of the forward cubic grid is identical to the origin of the tomographic
box. Numbers of points along the three directions (nx,ny,nz) as well as the unique grid stepping
h which is the same on both directions for cubic sampling. This stepping should be small enough
for making accurate travel-time computations.
We consider a sampling dray along the ray which could be similar or smaller than the eikonal
sampling h: you have often to tune it when difficulties are met for performing the two-points
ray tracing.Finally, one may define which nodes should be inverted. For the present time, we limit
ourselves to minimal indexes (ixo,iyo,izo) and maximal indexes (nix,niy,niz) which means that
we consider a smaller box inside the tomographic box for the inversion. Uptonow, we have
used full grid strategy and, for potential use, the origin of the eikonal grid should be different
from the inversion grid.
The basic input could be seen through a shell file for running this velocity model transformation.

#FILE model.head
% P velocity only (1) or P & S velocities (2)
1
% origin of the grid
-1000.00000
-375.000000
-1000.00000
% dimensions nx,ny,nz
65
4
25
strategy for 2D application
% samplings along x,y,z
500.000000
250.000000
500.000000
% forward problem sampling (cube)
250.000000
% ray sampling
10.0000000
% starting node
1 1 1
% ending node
65 4 25
END OF THE FILE

There are files you need to proceed sequentially through different steps defined by isolated
programs : they have files for communication between them.

Complete information: inside DOC_TOMOTV, look at the file

# USER_MANUAL_2015.pdf

# Conclusion

These codes are those designed around 1995 using first F77 writing. It has been converted
around 2002 using simple F90 writing.
Better understanding of these tomographic problems would certainly lead to better writing of
these codes as they are combined inputs from many students.
We may hope that more efficient algorithms could be available including topography,
anisotropy, attenuation … and so on.
Good luck …