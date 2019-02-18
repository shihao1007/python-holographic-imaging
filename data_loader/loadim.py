# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 10:29:34 2018

@author: shihao
"""

import numpy
from matplotlib import pyplot as plt
#import os, os.path


p_mean = list()
p_std = list()

imname = r'D:\irimages\irholography\15um-imstream-new\75000\sbf161_img_000001_1220'

im = numpy.fromfile(imname, dtype=numpy.uint16,count=-1,sep='')
im = numpy.reshape(im, (128, 128, 5000), order='F')

for inte_number in range (1, 25):
    s = (128,128,100)
    sub_im = numpy.zeros(s)
    pixel = numpy.zeros(99)
    for trial_index in range (0, 100):
        temp_im = numpy.zeros((128,128,trial_index+1))
        a = inte_number*trial_index
        b = (trial_index+1)*inte_number
        temp_im = im[:,:,a:b]
        sub_im[:,:,trial_index] = numpy.sum(temp_im, axis = 2)
    
    pixel = sub_im[62,62,:]
    p_mean.append(numpy.mean(pixel, axis = 0))
    p_std.append(numpy.std(pixel, axis = 0))

    
plt.figure(1)
plt.plot(p_mean)
plt.xlabel('Integration Time')
plt.ylabel('Mean')
plt.title('Plot of Mean along Integration Time')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([15, eff_inte, numpy.min(pixel_mean),numpy.max(pixel_mean)])
plt.grid(True)
plt.show()

plt.figure(2)
plt.plot(p_std)
plt.xlabel('Integration Time')
plt.ylabel('Standard Deviation')
plt.title('Plot of Standard Deviation along Integration Time')
plt.grid(True)
plt.show()


plt.figure(1)
plt.imshow(im[:,:,60], interpolation='nearest')
plt.show()

