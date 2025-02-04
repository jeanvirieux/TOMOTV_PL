##############################################################################
# cmd: python ../PYTHON/plot_dir.py dir   all files "mod.*"
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
  file_xy = str(sys.argv[1])
  step = 1
  bottom = -1.e+29
  top=1.e+29
else:
  print(len(sys.argv))
  print(' either 1 ')
  sys.exit()
  
############################
# plotting source and receiver
#############################
iplot=0
for file in glob.glob(file_xy+"/rms.*"):
   print('file ',file)
   iplot=iplot+1
   X_pos=[]
   Y_pos=[]
   infile = open(file,'r')
############ premiere lecture pour initialiser line
   iter=0
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
   if iplot == 1:
     plt.plot(X_pos,Y_pos,"o-",label="TV")
   if iplot == 2:
     plt.plot(X_pos,Y_pos,"x:",label="Smooth+TV2")
   if iplot == 3:
     plt.plot(X_pos,Y_pos,">-.",label="Smooth")
   if iplot == 4:
      plt.plot(X_pos,Y_pos,"<-",label="Penalty")
      
xmin=0
xmax=22
ymin=0.7
ymax=1.1

plt.legend()
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel('Iteration')
plt.ylabel('RMS (sec)')
plt.title('RMS Evolution')
plt.savefig(file_xy+"/total_rms.pdf",dpi=600,format='pdf')
plt.show()
