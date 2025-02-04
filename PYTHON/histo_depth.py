##############################################################################
# cmd: python ../PYTHON/histo.py file_x
# cmd: python ../PYTHON/histo.py file_x bottom top
##############################################################################
import sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

print('nbre of arguments',len(sys.argv))

if len(sys.argv) == 2 :
  file_x = str(sys.argv[1])
  step = 1
  bottom = -1.e+29
  top=1.e+29
elif len(sys.argv) == 4 :
  file_x = str(sys.argv[1])
  bottom = float(sys.argv[2])
  top = float(sys.argv[3])
  step = 1
else:
  print(len(sys.argv))
  print(' either 1 or 3 arguments ')
  sys.exit()

step=1
xval = []

#############################
infile = open(file_x,'r')

############ premiere lecture pour initialiser line
line=infile.readline()
coords = line.split()
if float(coords[0]) > bottom and float(coords[0]) < top:
   xval.append(float(coords[0]))
ilist=0

while line:
    line=infile.readline()
    ilist=ilist+1
    if line:
       if ilist % step == 0:
         coords = line.split()
         if float(coords[0]) > bottom and float(coords[0]) < top:
            xval.append(float(coords[0]))
infile.close()

xmin=min(xval)
xmax=max(xval)
print('min/max',xmin,xmax)

#x1 = [1, 2, 2, 3, 4, 4, 4, 4, 4, 5, 5, 5.2]

#inter = [-3.00, -2.75, -2.5, -2.25, -2.00, -1.75, -1.5, -1.25, -1.00, -0.75,
#            -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.00, 1.25, 1.5, 1.75, 2.00, 2.25, 2.5, 2.75, 3.00]

debut=int(min(xval))
fin=int(max(xval))
fine=50
            
inter = np.arange(debut,fin,(fin-debut)/fine)

plt.hist(xval, bins = inter, color = 'red', alpha=0.5,
            edgecolor = 'red', label = 'residues') # hatch = '\\',
#plt.hist(xval, bins = inter, color = 'blue', alpha=0.5,
#            edgecolor = 'red', hatch = '/', label = 'residues')

plt.ylabel('number of events')
plt.xlabel('Depth (km) ')
plt.title('Histogramme')

plt.show()

