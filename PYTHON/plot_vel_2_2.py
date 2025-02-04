#############################################
#
# python ../PYTHON/plot_vel_2_2.py 41 41 121 1 1 1 modelP.15 [2000 2200]
#
# model on binary format
# fevt and fsta in ascii format
##############################################

import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import glob
import re

##############################
# Model and dimension inputs #
##############################
if len(sys.argv) == 8 :
  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  nz = int(sys.argv[3])
  dx = float(sys.argv[4])
  dy = float(sys.argv[5])
  dz = float(sys.argv[6])
  model_vel = sys.argv[7] 
  speed = np.fromfile(str(model_vel), "float32")
  speed = speed.reshape(nz,ny,nx)
  vmin=np.amin(speed)
  vmax=np.amax(speed)
elif len(sys.argv) == 10 :
  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  nz = int(sys.argv[3])
  dx = float(sys.argv[4])
  dy = float(sys.argv[5])
  dz = float(sys.argv[6])
  model_vel = sys.argv[7] 
  vmin = float(sys.argv[8])
  vmax = float(sys.argv[9])
  speed = np.fromfile(str(model_vel), "float32")
  speed = speed.reshape(nz,ny,nx)
else:
  print(len(sys.argv))
  print(' either 7 or 9 arguments ')
  sys.exit()

xmin=0.
ymin=0.
zmin=0.
xmax=(nx-1)*dx
ymax=(ny-1)*dy
zmax=(nz-1)*dz

################ automatic level definition (index)
istep=int(nz/5)
iz1=istep
iz2=iz1+istep  
iz3=iz2+istep
iz4=iz3+istep
################ tunable level definition (index)
iz1=50
iz2=55
iz3=65
iz4=70

print('Four depth levels')
print('level1 (mm)',float(iz1-1)*dz)
print('level4 (mm)',float(iz4-1)*dz)

print('velocity max',vmax)
print('velocity min',vmin)

fig = plt.figure()

if len(sys.argv) == 8 :

 ax=fig.add_subplot(2,2,1)
 im=ax.imshow(speed[iz1,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz1-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(2,2,2)
 im=ax.imshow(speed[iz2,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz2-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(2,2,3)
 im=ax.imshow(speed[iz3,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.set_xlabel('in-line (mm)')
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz3-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(2,2,4)
 im=ax.imshow(speed[iz4,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.set_xlabel('in-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz4-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

elif len(sys.argv) == 10 :

 ax=fig.add_subplot(2,2,1)
 im=ax.imshow(speed[iz1,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz1-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(2,2,2)
 im=ax.imshow(speed[iz2,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz2-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(2,2,3)
 im=ax.imshow(speed[iz3,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.set_xlabel('in-line (mm)')
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz3-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(2,2,4)
 im=ax.imshow(speed[iz4,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.set_xlabel('in-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz4-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
cbar_ax.set_ylabel('velocity (mm/ms)')

fig.suptitle('Four Depth Sections')
plt.savefig("tomo_2_2.pdf",dpi=1200,format='pdf')
plt.show()
