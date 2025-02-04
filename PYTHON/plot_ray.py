#############################################
# cmd: python ../PYTHON/plot_ray.py fsrc.dat fsta.dat cyl_bottom.dat cyl_top.dat 5
##############################################
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import glob
import re
#import os
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

scale=1.   # units

fig = plt.figure()
ax = plt.axes(projection="3d")

xmin=0.
xmax=40.
ymin=0.
ymax=40.
zmin=0.
zmax=120. 

##############################
# Model and dimension inputs #
##############################
evt = sys.argv[1]
sta = sys.argv[2]
bas = sys.argv[3]
top = sys.argv[4]

step = sys.argv[5]

xpos = []

#############################
# plotting bas
#############################
infile = open(str(bas),"r")
######### read the first line
line=infile.readline()
n_bas=re.findall(r"[-+]?\d+", line) # extract real numbers from the string
print('number of events:',n_bas[0])
######### read other lines (positions: x,y,z)
ilist=0
while True:
 line=infile.readline()
 if not line:
   break
 ilist=ilist+1
 if ilist % 1 == 0:  # modulo 10
   result=re.findall(r"[-+]?\d*\.\d+|\d+", line)  # extract real numbers from the string
#   print(' coords',float(result[0]),float(result[1]),float(result[2]))
   ax.scatter3D(float(result[0])/scale,float(result[1])/scale,float(result[2])/scale,c='green',marker='.',s=2,edgecolors='none')
   xpos=np.append(xpos,float(result[0])/scale)
   xpos=np.append(xpos,float(result[1])/scale)
   xpos=np.append(xpos,float(result[2])/scale)

nt=int(len(xpos)/3)        # number of points
xpos_3=xpos.reshape(nt,3)  # reshape it
verts=[xpos_3]             # define vertices of a polygone

#################### transparency working around
pc=Poly3DCollection(verts, linewidths=1, edgecolors='r', alpha=.5)
pc.set_facecolor('C0')
ax.add_collection3d(pc)

#############################
# plotting top
#############################
infile = open(str(top),"r")
######### read the first line
line=infile.readline()
n_top=re.findall(r"[-+]?\d+", line) # extract real numbers from the string
print('number of events:',n_top[0])
######### read other lines (positions: x,y,z)
ilist=0
while True:
 line=infile.readline()
 if not line:
   break
 ilist=ilist+1
 if ilist % 1 == 0:  # modulo 10
   result=re.findall(r"[-+]?\d*\.\d+|\d+", line)  # extract real numbers from the string
#   print(' coords',float(result[0]),float(result[1]),float(result[2]))
   ax.scatter3D(float(result[0])/scale,float(result[1])/scale,float(result[2])/scale,c='green',marker='.',s=2,edgecolors='none')

#############################
# plotting evts
#############################
infile = open(str(evt),"r")
######### read the first line
line=infile.readline()
n_evt=re.findall(r"[-+]?\d+", line) # extract real numbers from the string
print('number of events:',n_evt[0])
######### read other lines (positions: x,z)
ilist=0
while True:
 line=infile.readline()
 if not line:
   break
 ilist=ilist+1
 if ilist % 1 == 0:  # modulo 10
   result=re.findall(r"[-+]?\d*\.\d+|\d+", line)  # extract real numbers from the string
#   print(' coords',float(result[0]),float(result[1]),float(result[2]))
   ax.scatter3D(float(result[0])/scale,float(result[1])/scale,float(result[2])/scale,c='blue',marker='x',s=20,edgecolors='none')

#############################
# plotting stations
#############################
infile = open(str(sta),"r")
######### read the first line
n_sta=infile.readline()
print('number of stations:',n_sta[0])
######### read other lines (positions: x,z)
ilist=0
while True:
 line=infile.readline()
 if not line:
   break
 ilist=ilist+1
 if ilist % 1 == 0:   # modulo 1
   result=re.findall(r"[-+]?\d*\.\d+|\d+", line) # extract real numbers from the string
#   print(' coords',float(result[0]),float(result[1]),float(result[2]))
   ax.scatter3D(float(result[0])/scale,float(result[1])/scale,float(result[2])/scale,c='red',marker='v',s=20,edgecolors='none')

#############################
# plotting rays from subdirectory TEMP_RAY
#############################

for file in glob.glob("PRAY"+"/SRC*"+".txt"):
#      print("event",chaine_evt)
    ilist=ilist+1
    if ilist % int(step) == 0:   # modulo step
      print(' fichier ',file)
      with open(file,'r') as f:
         x=[]    # declaration d'un vecteur
         y=[]
         z=[]
         ilist=0
         for line in f:
            result=re.findall(r"[-+]?\d*\.\d+|\d+", line)
#           plt.scatter(float(result[0])/scale,float(result[1])/scale,c='green',marker='.',s=5,edgecolors='none')
            x.append(float(result[0])/scale)    # ajout une valeur au vecteur x
            y.append(float(result[1])/scale)    # idem
            z.append(float(result[2])/scale)    # idem

         ax.plot3D(x,y,z,linewidth=0.5,linestyle='--',color='g',label='ray')

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_zlim(zmin,zmax)

plt.show()
