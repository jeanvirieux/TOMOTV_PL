#############################################
# cmd: python ../PYTHON/plot_vel_ray.py 41 41 121 1 1 1 modelP.15 3 10 [2000 2500]
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

##############################
# Model and dimension inputs #
##############################
if len(sys.argv) == 10 :
  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  nz = int(sys.argv[3])
  dx = float(sys.argv[4])
  dy = float(sys.argv[5])
  dz = float(sys.argv[6])
  model_vel = sys.argv[7] 
  n_section = int(sys.argv[8])
  step = int(sys.argv[9])
  speed = np.fromfile(str(model_vel), "float32")
  speed = speed.reshape(nz,ny,nx)
  vmin=np.amin(speed)
  vmax=np.amax(speed)
elif len(sys.argv) == 12 :
  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  nz = int(sys.argv[3])
  dx = float(sys.argv[4])
  dy = float(sys.argv[5])
  dz = float(sys.argv[6])
  model_vel = sys.argv[7] 
  n_section = int(sys.argv[8])
  step = int(sys.argv[9])
  vmin = float(sys.argv[10])
  vmax = float(sys.argv[11])
  speed = np.fromfile(str(model_vel), "float32")
  speed = speed.reshape(nz,ny,nx)
else:
  print(len(sys.argv))
  print(' either 9 or 11 arguments ')
  sys.exit()

speed = np.fromfile(str(model_vel), "float32")
speed = speed.reshape(nz,ny,nx)

xmin=0.
ymin=0.
zmin=0.
xmax=(nx-1)*dx
ymax=(ny-1)*dy
zmax=(nz-1)*dz

fig = plt.figure()
ax = plt.axes(projection="3d")

xx, yy = np.meshgrid(np.linspace(xmin,xmax,num=nx), np.linspace(ymin,ymax,num=ny))

print('velocity max',vmax)
print('velocity min',vmin)

value = []
for iz in range(0,100):
  value=np.append(value,vmin + (vmax-vmin)*iz/100)

#print(value)

dsec=0.5*(zmax-zmin)/float(n_section)
iz=int(nz/2)
zlevel=zmin+(iz-1)*dz
cset = ax.contourf(xx, yy, speed[iz,:,:], value, zdir='z', offset=zlevel, cmap=cm.seismic, vmin=vmin,vmax=vmax, alpha=0.2, extend='neither')

if n_section == 3:
 iz=int(nz/2+nz/4)
 zlevel=zmin+(iz-1)*dz
 cset = ax.contourf(xx, yy, speed[iz,:,:], value, zdir='z', offset=zlevel, cmap=cm.seismic, vmin=vmin,vmax=vmax, alpha=0.1, extend='neither')
 iz=int(nz/2-nz/4)
 zlevel=zmin+(iz-1)*dz
 cset = ax.contourf(xx, yy, speed[iz,:,:], value, zdir='z', offset=zlevel, cmap=cm.seismic, vmin=vmin,vmax=vmax, alpha=0.1, extend='neither')
else:
 print('only one section')

plt.colorbar(cset)

xpos = []

#############################
# plotting bas
#############################
infile = open("cyl_bottom.dat","r")
######### read the first line
line=infile.readline()
n_bas=re.findall(r"[-+]?\d+", line) # extract real numbers from the string
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
#   ax.scatter3D(float(result[0])/scale,float(result[1])/scale,float(result[2])/scale,c='green',marker='.',s=2,edgecolors='none')
   xpos=np.append(xpos,float(result[0])/scale)
   xpos=np.append(xpos,float(result[1])/scale)
   xpos=np.append(xpos,float(result[2])/scale)

nt=int(len(xpos)/3)        # number of points
xpos_3=xpos.reshape(nt,3)  # reshape it
verts=[xpos_3]             # define vertices of a polygone

#################### transparency working around
pc=Poly3DCollection(verts, linewidths=0, edgecolors='r', alpha=.25)
pc.set_facecolor('C0')
ax.add_collection3d(pc)

#############################
# plotting top
#############################
infile = open("cyl_top.dat","r")
######### read the first line
line=infile.readline()
n_top=re.findall(r"[-+]?\d+", line) # extract real numbers from the string
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
infile = open('fsrc.dat','r')
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
infile = open('fsta.dat','r')
######### read the first line
line=infile.readline()
n_sta=re.findall(r"[-+]?\d+", line) # extract real numbers from the string
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
#      print(' fichier ',file)
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
