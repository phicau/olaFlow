#!/usr/bin/python

import os 
from pylab import *
import numpy as np
import Image

def waveNumHunt(T, h):
    if type(T) == list:
        T = np.array(T)
    if type(h) == list:
        h = np.array(h)       
    
    D = np.array([0.6666666666, 0.3555555555, 0.1608465608, 0.0632098765, 0.0217540484, 0.0065407983])
    
    L0 = 9.81 * T**2 / (2.*np.pi)
    k0 = 2. * np.pi / L0
    k0h = k0 * h
    
    aux = D[0] * k0h + D[1] * k0h**2 + D[2] * k0h**3 + D[3] * k0h**4 + D[4] * k0h**5 + D[5] * k0h**6
    
    kh = k0h * np.sqrt( 1. + ( k0h * ( 1. + aux ))**-1  )
    
    k = kh / h
    L = 2. * np.pi / k
    return L

#H = 0.01
#h = 0.4
#T = 1.0

H = float(raw_input("Please enter H: "))
h = float(raw_input("Please enter h: "))
T = float(raw_input("Please enter T: "))

# Calculations
L = waveNumHunt(T, h)
print 'L =',L,'m'

g = 9.81

a = h/(g*pow(T,2.0))
b = H/(g*pow(T,2.0))

xMin = 5.34*pow(10,-4.0)
xMinPix = 402
xMax = 0.2
xMaxPix = 1977

yMin = 0.00005
yMinPix = 1865
yMax = 0.05
yMaxPix = 28

if a > xMax:
    a = xMax
elif a < xMin:
    a = xMin

if b > yMax:
    b = yMax
elif b < yMin:
    b = yMin

auxX = 612*log10(xMax/a)
auxY = 612*log10(b/yMin)

coordX = xMaxPix - auxX
coordY = yMinPix - auxY

# Plot
#figure(num=1, figsize=(40, 50), dpi=80)

pathname = os.path.abspath('.')
image=Image.open('teoria.png')
arr=np.asarray(image)
imshow(arr)

plot(coordX,coordY,'ro',markersize=12)
plot(coordX,coordY,'k.')
plot(coordX,coordY,'kx',markersize=12)


title('H = '+"%.2f"%H+' m; T = '+"%.2f"%T+' s; h = '+"%.2f"%h+' m.');

xlim(-30,2030)
ylim(2200,-30)

show()

