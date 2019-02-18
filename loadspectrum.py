# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 11:04:27 2018

@author: shihao
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import os, os.path
from matplotlib import animation as animation
import subprocess
import scipy.io
from pylab import *


numIm = 300
fileSeq = np.zeros((128, 128, numIm), dtype = np.uint16)
fileDir = r'D:\irimages\irholography\purged\QCL3scan2kFrames'
fileList = os.listdir(fileDir)                                      #get the list of images of the hologram
fileNum = len(fileList)                                             #get the number of images of the hologram
for i in range(fileNum):                                            #for each image
    imgDir = fileDir + '\\' + fileList[i]                           #get the directory of the image
    fileData = np.fromfile(imgDir, dtype = np.uint16, count = -1, sep = '')       #read in the image
    fileData = np.reshape(fileData, (128, 128), order = 'C')     #reshape it to 128x128
    fileSeq[:, :, i] = fileData                            #store it in a sequence
    

scipy.io.savemat(fileDir +'\\'+'spec.mat', {'spec':(fileSeq)})

#NIDir = r'D:\irimages\irholography\oldQCL\SNR-test-ambient\cover.txt'
#
#pointsC = []
#with open(NIDir) as file:
#    for line in file:
#        line = line.strip()
#        pointsC.append(float(line))
#
##fft = np.fft.fft(points)
#
#
#plt.figure()
##plt.plot(pointsC, 'r-', label = 'With Cover')
#plt.plot(pointsN, 'b-', label = 'Without Cover')
#
#plt.grid('on')
#plt.ylabel('Pixel Value')
#plt.xlabel('Frame Number')
#plt.title('FPA Detector Noise')
#plt.legend()

_min, _max = np.amin(np.abs(fileSeq)), np.amax(np.abs(fileSeq))

#plot and save the animaiton
fig = plt.figure()

img = []
for i in range(numIm):
    img.append([plt.imshow(np.abs(fileSeq[:,:,i]), vmax = _max, vmin = _min)])

ani = animation.ArtistAnimation(fig,img,interval=1000/numIm*16)
writer = animation.writers['ffmpeg'](fps=numIm/16) 
ani.save(r'D:\irimages\irholography\purged\QCL3scan2kFrames\scan.mp4',writer=writer)