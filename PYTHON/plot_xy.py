##############################################################################
# cmd: python ../PYTHON/plot_xy.py file_xy
# cmd: python ../PYTHON/plot_xy.py file_xy bottom top
# cmd: python ../PYTHON/plot_xy.py file_xy xmin xmax ymin ymax step
##############################################################################
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

print('nbre of arguments',len(sys.argv))
##############################
# Model and dimension inputs #
##############################
if len(sys.argv) == 2 :
  file_xy = str(sys.argv[1])
  step = 1
  bottom = -1.e+29
  top=1.e+29
elif len(sys.argv) == 4 :
  file_xy = str(sys.argv[1])
  step = 1
  bottom = float(sys.argv[2])
  top = float(sys.argv[3])
elif len(sys.argv) == 7 :
  file_xy = str(sys.argv[1])
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

X_pos=[]
Y_pos=[]
 
############################
# plotting source and receiver
#############################
infile = open(file_xy,'r')

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
if len(sys.argv) == 2 or len(sys.argv) == 4 :
 xmin=min(X_pos)
 xmax=max(X_pos)
 ymin=min(Y_pos)
 ymax=max(Y_pos)
 print('box',xmin,xmax,ymin,ymax)

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Y=f(X)')

plt.plot(X_pos,Y_pos,'b')
 
plt.savefig("XY_2D.pdf",dpi=600,format='pdf')
#plt.show()
