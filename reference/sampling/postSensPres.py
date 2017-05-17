#!/usr/bin/python

import os 

pathname = os.path.abspath('.')
savePath = os.path.join(pathname,'sensorsPres')
if not os.path.isdir(savePath):
    os.makedirs(savePath)

postPath = os.path.join(pathname,'postProcessing')
if os.path.isdir(postPath):
    postPath = 'postProcessing/sets'
else:
    postPath = 'sets'

# List of time dirs in order
a = os.listdir('./'+postPath)
a.sort(lambda a,b: cmp(float(a), float(b)))

# Get number of sensors
dir1 = os.path.join(pathname,postPath,a[int(len(a)/2.0)])
b = os.listdir(dir1)
nSens = 0
index = []
for i in range(len(b)):
    test1 = b[i].find('P') + 1
    test2 = b[i].find('p') + 1
    if test1 and test2:
        index.append(i)
        nSens += 1

first = True

for i in range(nSens):
    fileBaseName = b[index[i]][0:b[index[i]].find('_')] + '_'
    print 'Sensor ' + '%i' % int(i+1) + ' of ' + '%i' % nSens + '.'
    
    # Read files time by time
    for j in range(len(a)):
        directory = os.path.join(pathname,postPath,a[j])
        try:
            fileR = open(os.path.join(directory,b[index[i]]), 'r')
        except:
            print 'WARNING - File not present: ' + os.path.join(directory,b[index[i]])
        else:
            data = fileR.read()
            fileR.close()
            data = data.split('\n')
            
            if first: # First time step
                coord = j
                first = False
    
            # x y z pres : saving
            for k in range(len(data)-1):
                line = data[k]
                line = line.split('\t') 
                # x = float(line[0]) # y = float(line[1]) 
                # z = float(line[2]) # pres = float(line[3])
                fileName = fileBaseName + '%i' % (k+1)
    
                if j == coord: # First time step
                    # Create coordinate files
                    fileWXYZ = open(os.path.join(savePath,fileName + '.xyz'), 'w')
                    fileWXYZ.write( line[0] + line[1] + line[2] )
                    fileWXYZ.close()
                    # Reset files to write
                    fileW = open(os.path.join(savePath,fileName), 'w')  
                    fileW.close()
                    
                # Save sensors
                time = a[j]
                fileW = open(os.path.join(savePath,fileName), 'a') 
                fileW.write(time + ' ' + line[3] + '\n')
                fileW.close()
            
print 'Done'
