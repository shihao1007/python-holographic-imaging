# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 13:17:29 2018

@author: shihao
"""
import numpy
import os, os.path
from matplotlib import pyplot as plt
import subprocess

imfolder = r'C:\Users\shihao\Desktop\ir-images\ir-holography\new-inte-test'
ncap = str(100);                         #number of captured images passed into snap.exe 
inte_fold = 30;                         #effective integration time = 15 microsecond * inte_fold
pixel_mean = list()
pixel_std = list()
for i in range(1, inte_fold + 1):
    nframe = str(15 * i * 4)
    subprocess.call([r"qcl-snap.exe", imfolder, '--inte_time',nframe,'--ncap',ncap,'--WN','1220'])

    
    
for j in range (1, inte_fold):
    eff_inte = j * 15 * 4
    intefoldername = imfolder + '\\' + str(eff_inte)
    flst = os.listdir(intefoldername) # dir is your directory path
    number_files = len(flst)
    pixel = list()
    for i in range (0, number_files):
        imname = intefoldername + '\\' + flst[i]
        im = numpy.fromfile(imname, dtype=numpy.uint32,count=-1,sep='')
        im = numpy.reshape(im, (128, 128))
        pixel.append(im [60,60])
    
    
    pixel_mean.append(numpy.mean(pixel, axis = 0))
    pixel_std.append(numpy.std(pixel, axis = 0))    

plt.figure(1)
plt.plot(pixel_mean)
plt.xlabel('Integration Time')
plt.ylabel('Mean')
plt.title('Plot of Mean along Integration Time')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([15, eff_inte, numpy.min(pixel_mean),numpy.max(pixel_mean)])
plt.grid(True)
plt.show()

plt.figure(2)
plt.plot(pixel_std)
plt.xlabel('Integration Time')
plt.ylabel('Standard Deviation')
plt.title('Plot of Standard Deviation along Integration Time')
plt.grid(True)
plt.show()

plt.figure(3)
plt.imshow(im, interpolation='nearest')
plt.show()