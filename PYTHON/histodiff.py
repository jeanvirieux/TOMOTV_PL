##############################################################################
# cmd: python ../PYTHON/histodiff.py file1_dat file2_dat 
# cmd: python ../PYTHON/histodiff.py file1_dat file2_dat bottom top 
# cmd: python ../PYTHON/histodiff.py file1_dat file2_dat bottom top step
# avoid to count values outside [xmin,xmax]
# cmd: python ../PYTHON/histodiff.py file1_dat file2_dat bottom top step fine xmin xmax
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
  file1_dat = str(sys.argv[1])
  file2_dat = str(sys.argv[2])
  step = 1
  auto = 1
if len(sys.argv) == 5 :
  file1_dat = str(sys.argv[1])
  file2_dat = str(sys.argv[2])
  bottom = float(sys.argv[3])
  top = float(sys.argv[4])
  step = 1
  auto = 2
elif len(sys.argv) == 6 :
  file1_dat = str(sys.argv[1])
  file2_dat = str(sys.argv[2])
  bottom = float(sys.argv[3])
  top = float(sys.argv[4])
  step = int(sys.argv[5])
  auto = 3
elif len(sys.argv) == 9 :
  file1_dat = str(sys.argv[1])
  file2_dat = str(sys.argv[2])
  bottom = float(sys.argv[3])
  top = float(sys.argv[4])
  step = int(sys.argv[5])
  fine = int(sys.argv[6])
  xmin = float(sys.argv[7])
  xmax = float(sys.argv[8])
  auto = 4
else:
  print(len(sys.argv))
  print(' either 2, 4, 5, 8 arguments ')
  sys.exit()

xval = []
yval = []

#############################
infile = open(file1_dat,'r')

############ premiere lecture pour initialiser line
line=infile.readline()
coords = line.split()

# get the second value only
if auto == 4:
   if float(coords[0]) > xmin and float(coords[0]) < xmax:
      xval.append(float(coords[0]))
      
else:
   xval.append(float(coords[0]))

ilist=0

while line:
    line=infile.readline()
    ilist=ilist+1
    if line:
       if ilist % step == 0:
         coords = line.split()
         if auto == 4:
           if float(coords[0]) > xmin and float(coords[0]) < xmax:
              xval.append(float(coords[0]))
              
         else:
           xval.append(float(coords[0]))
infile.close()

#############################
infile = open(file2_dat,'r')

############ premiere lecture pour initialiser line
line=infile.readline()
coords = line.split()

# get the second value only
if auto == 4:
  if float(coords[0]) > xmin and float(coords[0]) < xmax:
     yval.append(float(coords[0]))
     
else:
  yval.append(float(coords[0]))
            
ilist=0

while line:
    line=infile.readline()
    ilist=ilist+1
    if line:
       if ilist % step == 0:
         coords = line.split()
         if auto == 4:
           if float(coords[0]) > xmin and float(coords[0]) < xmax:
             yval.append(float(coords[0]))
             
         else:
           yval.append(float(coords[0]))
infile.close()

if auto == 1 or auto == 2 or auto == 3:
  fine = 50

if auto == 1:
  debut=int(min(xval)*10)/10
  fin=int(max(xval)*10)/10
  if debut > int(min(yval)*10)/10:
    debut=int(min(yval)*10)/10
  if fin < int(max(yval)*10)/10:
    fin=int(max(yval)*10/10)

if auto == 2 or auto == 3 or auto == 4:
  debut=int(bottom*10)/10
  fin=int(top*10)/10

print(' choix des bornes',debut,fin)

print(np.shape(xval))
print(np.shape(yval))

inter = np.arange(debut,fin,(fin-debut)/fine)

#print(inter)

# inter= inter / 1000

plt.hist(xval, bins = inter, color = 'red', alpha=0.5,
            edgecolor = 'red', label = 'before') # hatch = '\\',

plt.hist(yval, bins = inter, color = 'blue', alpha=0.5,
            edgecolor = 'blue', label = 'after') # hatch = '/',

plt.ylabel('number of events')
plt.xlabel('Residual Time (sec)')
plt.title('Histogram')
plt.savefig("histodiff.pdf",dpi=600,format='pdf')

#plt.show()

