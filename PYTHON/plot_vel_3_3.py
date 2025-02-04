#############################################
#
# python ../PYTHON/plot_vel_3_3.py 41 41 121 1 1 1 modelP.15 [2000 2200]
#
# model on binary format
##############################################

import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
#import matplotlib.pylab as plt
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

istep=int(nz/9)
iz1=istep
iz2=iz1+istep  
iz3=iz2+istep
iz4=iz3+istep
iz5=iz4+istep
iz6=iz5+istep
iz7=iz6+istep
iz8=iz7+istep
iz9=iz8+istep

################### tuned levels 
istep=5
iz5=53
iz6=iz5+istep
iz7=iz6+istep
iz8=iz7+istep
iz9=iz8+istep
iz4=iz5-istep
iz3=iz4-istep
iz2=iz3-istep
iz1=iz2-istep

print('Nine vertical levels')
print('level1 (mm)',float(iz1-1)*dz)
print('level9 (mm)',float(iz9-1)*dz)

print('velocity max',vmax)
print('velocity min',vmin)

bbox_kwargs = {'fc': 'w', 'alpha': .75, 'boxstyle': "round4"}
ann_kwargs = {'xycoords': 'axes fraction','textcoords': 'offset points','bbox': bbox_kwargs}

fig = plt.figure(figsize=(10,10))
#fig = plt.figure()

if  len(sys.argv) == 8 :

 ax=fig.add_subplot(3,3,1)
 im=ax.imshow(speed[iz1,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz1-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,2)
 im=ax.imshow(speed[iz2,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.tick_params(axis='y',right=False,left=False,labelleft=False,labelright=False)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz2-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,3)
 im=ax.imshow(speed[iz3,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz3-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,4)
 im=ax.imshow(speed[iz4,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=False,labelbottom=False,labeltop=False)
 ax.tick_params(axis='y',right=False,left=True,labelleft=True,labelright=False)
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz4-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,5)
 im=ax.imshow(speed[iz5,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='y',right=False,left=False,labelleft=False,labelright=False)
 ax.tick_params(axis='x',top=False,bottom=False,labeltop=False,labelbottom=False)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz5-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,6)
 im=ax.imshow(speed[iz6,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=False,labelbottom=False,labeltop=False)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz6-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,7)
 im=ax.imshow(speed[iz7,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=False,left=True,labelleft=True,labelright=False)
 ax.set_xlabel('in-line (mm)')
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz7-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,8)
 im=ax.imshow(speed[iz8,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=False,left=False,labelleft=False,labelright=False)
 ax.set_xlabel('in-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz8-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,9)
 im=ax.imshow(speed[iz9,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.set_xlabel('in-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz9-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

elif  len(sys.argv) == 10 :
 ax=fig.add_subplot(3,3,1)
 im=ax.imshow(speed[iz1,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz1-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,2)
 im=ax.imshow(speed[iz2,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.tick_params(axis='y',right=False,left=False,labelleft=False,labelright=False)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz2-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,3)
 im=ax.imshow(speed[iz3,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz3-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,4)
 im=ax.imshow(speed[iz4,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=False,labelbottom=False,labeltop=False)
 ax.tick_params(axis='y',right=False,left=True,labelleft=True,labelright=False)
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz4-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,5)
 im=ax.imshow(speed[iz5,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='y',right=False,left=False,labelleft=False,labelright=False)
 ax.tick_params(axis='x',top=False,bottom=False,labeltop=False,labelbottom=False)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz5-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,6)
 im=ax.imshow(speed[iz6,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=False,labelbottom=False,labeltop=False)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz6-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,7)
 im=ax.imshow(speed[iz7,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=False,left=True,labelleft=True,labelright=False)
 ax.set_xlabel('in-line (mm)')
 ax.set_ylabel('x-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz7-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,8)
 im=ax.imshow(speed[iz8,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=False,left=False,labelleft=False,labelright=False)
 ax.set_xlabel('in-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz8-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

 ax=fig.add_subplot(3,3,9)
 im=ax.imshow(speed[iz9,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',vmin=vmin,vmax=vmax,extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
 ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
 ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
 ax.set_xlabel('in-line (mm)')
 ax.text(xmin+0.04*(xmax-xmin),ymin+0.94*(ymax-ymin), str(int(float(iz9-1)*dz))+"mm", size=8,ha="left",va="top",bbox=dict(boxstyle="round",ec=(0.0, 1.0, 0.0),fc=(0.0,1.0, 0.0)))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
cbar_ax.set_ylabel('velocity (mm/ms)')

fig.suptitle('Nine Depth Sections')
plt.savefig("tomo_3_3.pdf",dpi=600,format='pdf')
plt.show()
