# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 13:38:45 2018

@author: shihao
"""

"""
find the center burst of the QCL Channerl 2
"""

import numpy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os, os.path
from matplotlib import animation as animation
from pylab import *
import scipy.io

#specify all imaging parameters
data_dir = r'D:\irimages\irholography\spectrum'  #directory of all raw data
zstep = 0.25               #stage moving stepsize, delta_m, lightpath length difference in the hologram
num_frames = 15         #number of averaged frames for each image, to minimize laser fluctuation
positions = 1          #number of different positions of the stage. i.e. number of measurements
lp = 81                 #laser power
inte = 30               #integration time of FPA
wn_start = 1420         #starting wavenumber of the spectrum
wn_end = 1660           #ending wavenumber of the spectrum
wn_num = 240             #number of bands to be imaged
#badBands = [1505, 1532, 1557, 1558, 1559, 1633, 1634, 1645, 1650, 1651, 1652, 1661, 1667, 1673]
realBands = numpy.loadtxt(r'D:\irimages\irholography\spectrum\realBands.txt')
powerSeq = numpy.loadtxt(r'D:\irimages\irholography\spectrum\powerSeq.txt')
##load in holograms into dataset
def readFile(fileType):
    #fileType: specify the input hologram type, can be sample or bg(background) or ref(reference)
    fileSeq = numpy.zeros((128, 128, positions, wn_num))                    #initialize an array to store hologram
    for wnIndex in range(wn_num):                                           #for each wavenumber
        wn = int(wn_start + wnIndex * (wn_end - wn_start) / wn_num)         #calculate the current wavenumber1
        fileDir = data_dir + '\\' + str(fileType) + '\\' + str(wn)          #get the directory of the wavenumber folder
        fileList = os.listdir(fileDir)                                      #get the list of images of the hologram
        fileNum = len(fileList)                                             #get the number of images of the hologram
        for i in range(fileNum):                                            #for each image
            imgDir = fileDir + '\\' + fileList[i]                           #get the directory of the image
            fileData = numpy.fromfile(imgDir, dtype = numpy.uint16, count = -1, sep = '')       #read in the image
            fileData = numpy.reshape(fileData, (128, 128), order = 'C')     #reshape it to 128x128
            fileSeq[:, :, i, wnIndex] = fileData                            #store it in a sequence
        
    return fileSeq

dataSeq5 = readFile('powerissue5')

dataSeq7 = readFile('powerissue7')

wnX = numpy.arange(1420, 1660, 1)
subtractionBands = realBands - wnX
mean5 = numpy.zeros((wn_num))
mean7 = numpy.zeros((wn_num))

for i in range(wn_num):
    mean5[i] = numpy.mean(dataSeq5[:,:,0,i])
    mean7[i] = numpy.mean(dataSeq7[:,:,0,i])

_min, _max = numpy.amin(mean5), numpy.amax(mean5)  

plt.figure()
plt.subplot(311)    
plt.plot(wnX,mean5,'b-',label = 'Before (threshold 12000)')
plt.plot(wnX,mean7,'y-',label = 'After (threshold 11500)')
plt.xlabel('wavenumber')
plt.ylabel('Intensity')
plt.ylim(_min-100, _max)
plt.legend()

plt.subplot(312)
plt.plot(wnX, powerSeq, 'm-', label = 'Laser Power')
plt.xlabel('wavenumber')
plt.ylabel('Power / %')
plt.legend()

plt.subplot(313)
plt.plot(wnX, subtractionBands, 'r-', label = 'Wavenumber Shifts')
plt.xlabel('wavenumber')
plt.ylabel('wavenumber')
plt.legend()


centralPixel = dataSeq[64, 64, 0, :]
#sumWave = sum(centralPixel, axis = 1)

x = numpy.arange(0, 30, 1)

plt.figure()
plt.plot(x, centralPixel)

for i in range(10):

    plt.subplot(2,5,i+1)
    plt.plot(x, centralPixel[:,i])
    plt.title(str(i) + 'th wavenumber')


#get the maximal and minimal values of all images to turn off auto-scaling
_min, _max = numpy.amin(dataSeq7[:,:,0,:]), numpy.amax(dataSeq7[:,:,0,:])    


#plot the hologram and write out a .mp4 file
fig = plt.figure()
img = []
for i in range(240):
    img.append([plt.imshow(dataSeq7[:,:,0,i], vmin = _min, vmax = _max)])

ani = animation.ArtistAnimation(fig,img,interval=50)
writer = animation.writers['ffmpeg'](fps=20)

ani.save(data_dir + '\\afterCompen.mp4',writer=writer)
plt.figure() 
for j in range(10):

    plt.subplot(2,5,j+1)
    posIndex = j
    
    meanSpec = numpy.zeros((wn_num))
    for i in range(wn_num):
        meanSpec[i] = numpy.mean(imageSeqSub[:,:,posIndex,i])
        
    _min, _max = numpy.amin(meanSpec), numpy.amax(meanSpec)
    wnX = numpy.arange(1420, 1660, (wn_end - wn_start) / wn_num)
 
    plt.plot(wnX,meanSpec,'b-',label = 'Difference')
    plt.xlabel('wavenumber')
    plt.ylabel('Intensity')
    plt.ylim(_min-100, _max)
    plt.legend()
