#!/usr/bin/python

import os 
from pylab import *

pathname = os.path.abspath('.')
readPath = os.path.join(pathname,'sensorsPres')
a = os.listdir(readPath)

# Sorting
remove = []
for i in range(len(a)):
    if (a[i].rfind('.')+1): # Includes point
        remove.append(i)
remove.reverse()
for i in remove:
    a.pop(i)
a.sort(lambda a,b: cmp(int(a.split('_')[1]), float(b.split('_')[1])))

# Plot  
index = 0
indexFig = 0
for gauge in a:
    index = index + 1
    if index >= 4 or index == 1:
        index = 1
        indexFig = indexFig + 1
        figure(num=indexFig)
        
    subplots_adjust(hspace=0.6)
    
    fileR = open(os.path.join(readPath,gauge), 'r')
    data = fileR.read()
    fileR.close()
    data = data.split('\n')
    x = []
    y = []
    for i in range(len(data)-1):
        line = data[i]
        line = line.split(' ') 
        x.append(float(line[0]))
        y.append(float(line[1]))
  
    subplot(3,1,index)
    plot(x,y)
    xlabel('t (s)')
    ylabel('P (Pa)')
    title(gauge)

show()
