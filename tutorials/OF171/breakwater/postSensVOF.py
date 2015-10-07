#!/usr/bin/python

import os 

pathname = os.path.abspath('.')
savePath = os.path.join(pathname,'gaugesVOF')
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
    test1 = b[i].find('VOF') + 1
    test2 = b[i].find('alpha') + 1
    if test1 and test2:
        index.append(i)
        nSens += 1

first = True

for i in range(nSens):
    # Create files to write
    fileName = b[index[i]][0:b[index[i]].find('_')]
    fileW = open(os.path.join(savePath,fileName), 'w')
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
            
            x = []
            y = []
            z = []
            alpha = []
    
            # x y z alpha1 calculation
            for k in range(len(data)-1):
                line = data[k]
                line = line.split('\t') 
                # x = float(line[0]) # y = float(line[1]) 
                # z = float(line[2]) # pres = float(line[3])
                
                z.append(float(line[2]))
                alpha.append(float(line[3]))
                
                if j == coord: # First time step
                    # Create coordinate files
                    fileWXYZ = open(os.path.join(savePath,fileName + '.xy'), 'w')
                    fileWXYZ.write( line[0] + line[1] )
                    fileWXYZ.close()
    
            # Integrate in Z
            wLevel = z[0]
            for k in range(len(z)-1):
                wLevel = wLevel + alpha[k]*(z[k+1]-z[k])
    
            # Write to file
            time = a[j]
            fileW.write(time + ' ' + '%.6f' % wLevel + '\n')

    fileW.close()

print 'Done'
