#!/usr/bin/python

import os 
from pylab import *

# Sub-functions
def getNumData(fileName):
    fileR = open(fileName, 'r')
    data = fileR.read()
    fileR.close()
    data = data.split('\n')
    data.pop(0) # Header
    data.pop(0) # Header
    data.pop() # Last line
    data.pop() # Last line
    x = []; y = [];
    for line in data:
        # print fileName
        # print line
        line = line.split(' ') 
        x.append(float(line[0]))
        y.append(float(line[2]))
    return [x,y]
    
def getLabData(fileName):   
    fileR = open(fileName, 'r')
    data = fileR.read()
    fileR.close()
    data = data.split('\n')
    data.pop() # Last line
    x = []; y = [];
    for line in data:
        line = line.split(' ')
        x.append(float(line[0]))
        y.append(float(line[1]))
    return [x,y]
    
# Definitions
times = ['0.00', '0.35', '0.75', '1.15', '1.55', '1.95'];
pathNum = './freeSurface';
fileNum = 'alpha1_topFreeSurface.raw'
pathLab = './labData';

# Plot
figure(num=1)
subplots_adjust(hspace=0.6)

for i in range(len(times)):
    subplot(2,3,i+1)
    # Porous medium shading
    fill([0.3, 0.3, 0.59, 0.59],[0.0, 0.5, 0.5, 0.0], '0.5')
    # Numerical data
    if(i!=0):
        [x,y] = getNumData(os.path.join(pathNum,times[i],fileNum))
        plot(x,y,'r.')
    # Experimental data
    fileLab = 'b%i.dat' % (i+1)
    [x,y] = getLabData(os.path.join(pathLab,fileLab))
    if(i!=0):
        plot(x,y,'ko')
    else:
        plot(x,y,'b-',linewidth=4)
    
    # Cosmetics
    title('t = {0} s'.format(times[i]))
    xlabel('x (m)')
    ylabel('y (m)')
    xlim(0, 0.9)
    ylim(0.0, 0.5)
    if(i!=0):
        legend(('Num', 'Lab'),loc=1)

suptitle('IHFOAM validation. Lin (1998) dam break experiments: crushed rocks, h = 0.35 m.', fontsize=22)
show()