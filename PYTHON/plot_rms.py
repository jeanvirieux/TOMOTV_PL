##############################################################################
# cmd: python ../PYTHON/plot_rms.py repertoire
# cmd: python ../PYTHON/plot_rms.py repertoire bottom top
# cmd: python ../PYTHON/plot_rms.py repertoire xmin xmax ymin ymax step
##############################################################################
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import glob

print('nbre of arguments',len(sys.argv))
##############################
# Model and dimension inputs #
##############################
if len(sys.argv) == 2 :
  repert = str(sys.argv[1])
  step = 1
  bottom = -1.e+29
  top=1.e+29
  xmin=1.e+29
  xmax=-1.e+29
  ymin=1.e+29
  ymax=-1.e+29
elif len(sys.argv) == 4 :
  repert = str(sys.argv[1])
  step = 1
  bottom = float(sys.argv[2])
  top = float(sys.argv[3])
  xmin=1.e+29
  xmax=-1.e+29
  ymin=1.e+29
  ymax=-1.e+29
elif len(sys.argv) == 7 :
  repert = str(sys.argv[1])
  xmin = float(sys.argv[2])
  xmax = float(sys.argv[3])
  ymin = float(sys.argv[4])
  ymax = float(sys.argv[5])
  step = int(sys.argv[6])
  bottom = ymin
  top = ymax
else:
  print(len(sys.argv))
  print(' either 1 or 3 or 6 arguments ')
  sys.exit()

############################
# plotting source and receiver
#############################
icolor = -1
for file in glob.glob(repert+"/rms.*"):
   print('file ',file)
   icolor=icolor+1
   X_pos=[]
   Y_pos=[]
   infile = open(file,'r')
############ premiere lecture pour initialiser line
   line=infile.readline()
   coords = line.split()
   if float(coords[1]) > bottom and float(coords[1]) < top:
     X_pos.append(float(coords[0]))
     Y_pos.append(float(coords[1]))
############ boucle sur le fichier
     ilist=0
     while line:
       line=infile.readline()
       ilist=ilist+1
       if line:
          if ilist % step == 0:
            coords = line.split()
            if float(coords[1]) > bottom and float(coords[1]) < top:
              X_pos.append(float(coords[0]))
              Y_pos.append(float(coords[1]))
              
   infile.close()
   if icolor == 0 :
     plt.plot(X_pos,Y_pos,'b')
   if icolor == 1 :
     plt.plot(X_pos,Y_pos,'r')
   if icolor == 2 :
     plt.plot(X_pos,Y_pos,'g')

############ il faut calculer les extremums
   if len(sys.argv) == 2 or len(sys.argv) == 4 :
     xmin=min(X_pos)
     xmax=max(X_pos)
     ymin=min(Y_pos)
     ymax=max(Y_pos)
     print('box',xmin,xmax,ymin,ymax)

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

plt.xlabel('Iteration')
plt.ylabel('Weighted data misfit')
plt.title('Misfit Evolution')

plt.savefig("Misfit.pdf",dpi=600,format='pdf')
#plt.show()
