# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 14:57:58 2018

@author: shihao
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 10:29:34 2018

@author: shihao
"""

import numpy
from matplotlib import pyplot as plt
import os, os.path
from matplotlib import cm
import subprocess

imfolder = r'C:\Users\shihao\Desktop\ir-images\ir-holography\3d-test-1\3'
ltnum = 10                               #number if laser tuning
inilp = 52                               #initial laser power
imnum = 20000                             #number of images taken for each laser power
ncap = str(imnum)
trialnum = 100                           #number of trials for each integraiton time
maxintetime = int(imnum / trialnum)          #maximal number of integrated images
#for i in range(1, ltnum):
#    lp = str(inilp + i * 2)
#    subprocess.call([r"qcl-snap-gns.exe", imfolder,'--laserpower',lp,'--ncap',ncap])

p = (maxintetime,ltnum)
p_mean_lp = numpy.zeros(p)
p_std_lp = numpy.zeros(p)
s = (128,128,imnum,ltnum)
im = numpy.zeros(s)

for j in range(1, ltnum + 1):
    lp = str(inilp + j * 2)
    folername = imfolder + '\\' + lp
    flst = os.listdir(folername)
    number_files = len(flst)
    for i in range (0, number_files):
        imname = folername + '\\' + flst[i]
        indi_m = numpy.fromfile(imname, dtype=numpy.uint16,count=-1,sep='')
        indi_m = numpy.reshape(indi_m, (128, 128), order='C')
        im[:,:,i,j-1] = indi_m
    

for j in range(0, ltnum):
    for inte_number in range (1, maxintetime + 1):
        s = (128,128,trialnum,ltnum)
        sub_im = numpy.zeros(s)
        pixel = numpy.zeros(trialnum)
        for trial_index in range (0, trialnum):
            temp_im = numpy.zeros((128,128,trial_index+1))
            a = inte_number*trial_index
            b = (trial_index+1)*inte_number
            temp_im = im[:,:,a:b,j]
            sub_im[:,:,trial_index,j] = numpy.sum(temp_im, axis = 2)
        
        pixel = sub_im[62,62,:,j]
        p_mean_lp[inte_number - 1, j] = numpy.mean(pixel, axis = 0)
        p_std_lp[inte_number - 1, j] = numpy.std(pixel, axis = 0)

    
#plt.figure(1)
#plt.plot(p_mean)
#plt.xlabel('Integration Time')
#plt.ylabel('Mean')
#plt.title('Plot of Mean along Integration Time')
#plt.axis([0.015, 4.5, 0, 3000])
##plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
##plt.axis([15, eff_inte, numpy.min(pixel_mean),numpy.max(pixel_mean)])
#plt.grid(True)
#plt.show()
#
#plt.figure(2)
#plt.plot(p_std)
#plt.xlabel('Integration Time')
#plt.ylabel('Standard Deviation')
#plt.title('Plot of Standard Deviation along Integration Time')
#plt.grid(True)
#plt.show()
#
#
#plt.figure(1)
#plt.imshow(im[:,:,60], interpolation='nearest')
#plt.show()
        
x = numpy.arange(0,ltnum - 4, 1)
y = numpy.arange(0, maxintetime, 1)
intetime = y * 15
laserpower = 52 + x * 2
xx,yy = numpy.meshgrid(laserpower,intetime)

fig1 = plt.figure(1)
amean = fig1.gca(projection='3d')
amean.set_xlabel('Laser Power/%')
amean.set_ylabel('Integration Time/us')
amean.set_zlabel('Intensity')
surf1 = amean.plot_surface(xx, yy, p_mean_lp_u, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig1.colorbar(surf1, shrink=0.5, aspect=5)
plt.show()

fig2 = plt.figure(2)
astd = fig2.gca(projection='3d')
astd.set_xlabel('Laser Power/%')
astd.set_ylabel('Integration Time/us')
astd.set_zlabel('Intensity')
surf2 = astd.plot_surface(xx, yy, p_std_lp_u, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig2.colorbar(surf2, shrink=0.5, aspect=5)
plt.show()

fig3 = plt.figure(3)
asnr = fig3.gca(projection='3d')
asnr.set_xlabel('Laser Power/%')
asnr.set_ylabel('Integration Time/us')
asnr.set_zlabel('SNR')
surf3 = asnr.plot_surface(xx, yy, p_snr_lp_u, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig3.colorbar(surf3, shrink=0.5, aspect=5)
plt.show()