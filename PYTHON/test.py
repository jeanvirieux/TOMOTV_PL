#############################################
#
# python voir.py 51 200 model.bin fsta.asc fevt.asc fevt_ref.asc
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
if len(sys.argv) == 9 :
  nx = int(sys.argv[1])
  ny = int(sys.argv[2])
  nz = int(sys.argv[3])
  dx = float(sys.argv[4])
  dy = float(sys.argv[5])
  dz = float(sys.argv[6])
  model_vel = sys.argv[7] 
  n_section = int(sys.argv[8])
  vmin=2000.
  vmax=3000.
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
else:
  print(len(sys.argv))
  print(' either 8 or 10 arguments ')
  sys.exit()

xmin=0.
ymin=0.
zmin=0.
xmax=(nx-1)*dx
ymax=(ny-1)*dy
zmax=(nz-1)*dz


speed = np.fromfile(str(model_vel), "float32")
speed = speed.reshape(nz,ny,nx)

fig = plt.figure()

ax=fig.add_subplot(2,2,1)
im=ax.imshow(speed[30,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
ax.set_ylabel('Offset (m)')

ax=fig.add_subplot(2,2,2)
im=ax.imshow(speed[30,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
ax.tick_params(axis='x',top=True,bottom=False,labelbottom=False,labeltop=True)
ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)

ax=fig.add_subplot(2,2,3)
im=ax.imshow(speed[30,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
#ax.xlabel('Offset (m)')
#ax.ylabel('Depth (m)')
ax.set_xlabel('Depth (m)')
ax.set_ylabel('Offset (m)')

ax=fig.add_subplot(2,2,4)
im=ax.imshow(speed[30,:,:],origin='upper',aspect='equal',interpolation='bilinear',resample='true',extent=[xmin,xmax,ymin,ymax],cmap=cm.seismic,clip_on=True)
#ax.xlabel('Offset (m)')
ax.tick_params(axis='x',top=False,bottom=True,labelbottom=True,labeltop=False)
ax.tick_params(axis='y',right=True,left=False,labelleft=False,labelright=True)
ax.set_xlabel('Depth (m)')

#plt.subplots_adjust(bottom=0.1, right=2.9, top=1.9)
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

#fig.colorbar(im)
fig.suptitle('Toy Example')
plt.show()
