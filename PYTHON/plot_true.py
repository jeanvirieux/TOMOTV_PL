#############################################
# cmd: python ../PYTHON/plot_vel.py 41 41 121 1 1 1 modelP.15 3 [2000 2500]
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
if len(sys.argv) == 9 :
  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  nz = int(sys.argv[3])
  dx = float(sys.argv[4])
  dy = float(sys.argv[5])
  dz = float(sys.argv[6])
  model_vel = sys.argv[7] 
  n_section = int(sys.argv[8])
  speed = np.fromfile(str(model_vel), "float32")
  speed = speed.reshape(nx,ny,nz)
  vmax=np.amax(speed)
  vmin=np.amin(speed)
elif len(sys.argv) == 11 :
  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  nz = int(sys.argv[3])
  dx = float(sys.argv[4])
  dy = float(sys.argv[5])
  dz = float(sys.argv[6])
  model_vel = sys.argv[7] 
  n_section = int(sys.argv[8])
  vmin = float(sys.argv[9])
  vmax = float(sys.argv[10])
  speed = np.fromfile(str(model_vel), "float32")
  speed = speed.reshape(nx,ny,nz)
else:
  print(len(sys.argv))
  print(' either 8 or 10 arguments ')
  sys.exit()

print('velocity max',vmax)
print('velocity min',vmin)

xmin=-14.
ymin=-14.
zmin=-1.
xmax=xmin+(nx-1)*dx
ymax=ymin+(ny-1)*dy
zmax=zmin+(nz-1)*dz

print('box x',xmin,xmax)
print('box y',ymin,ymax)
print('box z',zmin,zmax)

fig = plt.figure()
ax = plt.axes(projection="3d")

xx, yy = np.meshgrid(np.linspace(xmin,xmax,num=nx), np.linspace(ymin,ymax,num=ny))

value = []
for iz in range(0,100):
  value=np.append(value,vmin + (vmax-vmin)*iz/100)

#print(value)

dsec=0.5*(zmax-zmin)/float(n_section)

if n_section == 3:
 iz=int(nz/2+nz/4)
 zlevel=zmin+(iz-1)*dz
 cset = ax.contourf(xx, yy, speed[:,:,iz], levels=value, zdir='z', offset=zlevel, cmap=cm.seismic, vmin=vmin,vmax=vmax, alpha=0.1, extend='neither')
 iz=int(nz/2)
 zlevel=zmin+(iz-1)*dz
 cset = ax.contourf(xx, yy, speed[:,:,iz], levels=value, zdir='z', offset=zlevel, cmap=cm.seismic, vmin=vmin,vmax=vmax, alpha=0.5, extend='neither')
 iz=int(nz/2-nz/4)
 zlevel=zmin+(iz-1)*dz
 cset = ax.contourf(xx, yy, speed[:,:,iz], levels=value, zdir='z', offset=zlevel, cmap=cm.seismic, vmin=vmin,vmax=vmax, alpha=0.1, extend='neither')
else:
 print('only one section')
 iz=int(nz/2)
 zlevel=zmin+(iz-1)*dz
 cset = ax.contourf(xx, yy, speed[iz,:,:], levels=value, zdir='z', offset=zlevel, cmap=cm.seismic, vmin=vmin,vmax=vmax, alpha=0.8, extend='neither')

ax.plot_wireframe(xx, yy, speed[:,:,29], color='b')
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
   print(' coords',float(result[0]),float(result[1]),float(result[2]))
   ax.scatter3D(float(result[0])/scale,float(result[1])/scale,float(result[2])/scale,c='green',marker='.',s=2,edgecolors='none')


#############################
# plotting shots (or true quakes)
#############################
infile = open('fshots.dat','r')
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
   print(' coords',float(result[0]),float(result[1]),float(result[2]))
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
   print(' coords',float(result[0]),float(result[1]),float(result[2]))
   ax.scatter3D(float(result[0])/scale,float(result[1])/scale,float(result[2])/scale,c='red',marker='v',s=30,edgecolors='none')


ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_zlim(zmin,zmax)
ax.view_init(elev=10., azim=145.)

#plt.savefig("tomo_events.pdf",dpi=600,format='pdf')

plt.show()
