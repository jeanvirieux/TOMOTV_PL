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
  xmin_t=1.e+29
  xmax_t=-1.e+29
  ymin_t=1.e+29
  ymax_t=-1.e+29
else:
  print(len(sys.argv))
  print(' 1 argument ')
  sys.exit()

############################
# plotting source and receiver
#############################
icolor = -1
for file in glob.glob(repert+".*"):
   print('file ',file)
   icolor=icolor+1
   X_pos=[]
   Y_pos=[]
   infile = open(file,'r')
   iter=0
############ premiere lecture pour initialiser line
   line=infile.readline()
   coords = line.split()
   if float(coords[0]) > bottom and float(coords[0]) < top:
     Y_pos.append(float(coords[0]))
     X_pos.append(float(iter))
############ boucle sur le fichier
     ilist=0
     while line:
       line=infile.readline()
       ilist=ilist+1
       if line:
          if ilist % step == 0:
            coords = line.split()
            if float(coords[0]) > bottom and float(coords[0]) < top:
              iter=iter+1
              Y_pos.append(float(coords[0]))
              X_pos.append(float(iter))
              
   infile.close()
   plt.plot(X_pos,Y_pos,"o-",label=file)

############ il faut calculer les extremums
   if len(sys.argv) == 2 :
     xmin=min(X_pos)
     xmax=max(X_pos)
     ymin=min(Y_pos)
     ymax=max(Y_pos)
     print('box',xmin,xmax,ymin,ymax)

     xmin_t=min(xmin,xmin_t)
     xmax_t=max(xmax,xmax_t)
     ymin_t=min(ymin,ymin_t)
     ymax_t=max(ymax,ymax_t)

     ymin_t=0.52
     ymax_t=0.79
     
plt.xlim(xmin_t,xmax_t)
plt.ylim(ymin_t,ymax_t)

plt.legend()
plt.xlabel('Iteration')
plt.ylabel('Weighted data misfit')
plt.title('Misfit Evolution')

plt.savefig("Misfit.pdf",dpi=600,format='pdf')
#plt.show()
