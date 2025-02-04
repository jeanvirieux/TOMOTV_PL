##############################################################################
# cmd: python ../PYTHON/plot_scatter.py file_ini file_final
#               python plot_misfit.py misfit.ini misfit.final
##############################################################################
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

print('nbre of arguments',len(sys.argv))
##############################
# Model and dimension inputs #
##############################
if len(sys.argv) == 3 :
  file_ini = str(sys.argv[1])
  file_final = str(sys.argv[2])
  step = 1
  bottom = -1.e+29
  top=1.e+29
elif len(sys.argv) == 5 :
  file_ini = str(sys.argv[1])
  file_final = str(sys.argv[2])
  step = 1
  bottom = float(sys.argv[3])
  top = float(sys.argv[4])
elif len(sys.argv) == 8:
  file_ini = str(sys.argv[1])
  file_final = str(sys.argv[2])
  xmin = float(sys.argv[3])
  xmax = float(sys.argv[4])
  ymin = float(sys.argv[5])
  ymax = float(sys.argv[6])
  step = int(sys.argv[7])
  bottom = ymin
  top = ymax
else:
  print(len(sys.argv))
  print(' either 1 or 3 or 6 arguments ')
  sys.exit()

X_pos=[]
Y_pos=[]
 
############################
# plotting source and receiver
#############################
infile = open(file_ini,'r')

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
############ il faut calculer les extremums
if len(sys.argv) == 3 or len(sys.argv) == 5 :
 xmin=min(X_pos)
 xmax=max(X_pos)
 ymin=min(Y_pos)
             
 ymin=0.45
             
 ymax=max(Y_pos)
 print('box',xmin,xmax,ymin,ymax)


plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

plt.xlabel('Number of models')
plt.ylabel('Weighted Data Misfit (L2-RMS)')
plt.title('Initial and Final Misfits')

plt.scatter(X_pos,Y_pos,s=1,c='r')

X_pos=[]
Y_pos=[]
 
############################
# plotting source and receiver
#############################
infile = open(file_final,'r')

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
             
plt.scatter(X_pos,Y_pos,s=1,c='b')

 
plt.savefig("misfit_models.pdf",dpi=600,format='pdf')
#plt.show()
