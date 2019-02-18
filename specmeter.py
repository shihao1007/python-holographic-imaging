# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 15:10:57 2018
read in sample and bg images from the single beam spectrumeter
calculate the absorbance image of each bands and save as a mat file
@author: shihao
"""
import numpy
from matplotlib import pyplot as plt
import os, os.path
from matplotlib import animation as animation
import subprocess
import scipy.io
from pylab import *


#specify all imaging parameters
data_dir = r'D:\irimages\irholography\3rdQCL\firstVerification'  #directory of all raw data
zstep = 0.25               #stage moving stepsize, delta_m, lightpath length difference in the hologram
num_frames = 50         #number of averaged frames for each image, to minimize laser fluctuation
positions = 1          #number of different positions of the stage. i.e. number of measurements
ref_pos = 16
lp = 81                 #laser power
inte = 30               #integration time of FPA
wn_start = 912        #starting wavenumber of the spectrum
wn_end = 1910           #ending wavenumber of the spectrum
wn_num = 1070             #number of bands to be imaged
#realBands = numpy.loadtxt(r'D:\irimages\irholography\spectrum\realBands.txt')
#powerSeq = numpy.loadtxt(r'D:\irimages\irholography\spectrum\powerSeq.txt')

dnf_seq = numpy.zeros((128, 128))
ab= numpy.zeros((128, 128, wn_num))
#read in detector noise
dn_dir = r'D:\irimages\irholography\spectrum\DN'
dn_list = os.listdir(dn_dir)
dnf_dir = dn_dir + '\\' + dn_list[0]
dnf_data = numpy.fromfile(dnf_dir, dtype = numpy.uint16, count = -1, sep = '')
dnf_data = numpy.reshape(dnf_data, (128, 128), order = 'C')
dnf_seq[:, :] = dnf_data

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

#read in raw data and background data
imageSeqRaw = readFile('newQCL')
#bgSeq = readFile('bg')

imageSeqRaw = numpy.squeeze(imageSeqRaw)
bgSeq = numpy.squeeze(bgSeq)
#mean = numpy.zeros((wn_num))

#for i in range(wn_num):
##    imageSeqRaw[:,:,i] = imageSeqRaw[:,:,i] - dnf_seq
##    bgSeq[:,:,i] = bgSeq[:,:,i] - dnf_seq
#    ab[:,:,i] = -numpy.log10(imageSeqRaw[:,:,i]/bgSeq[:,:,i])
##    mean[i] = numpy.mean(bgSeq[:,:,i])
    
scipy.io.savemat(data_dir +'\\'+'ab.mat', {'ab':(ab)})
    
#D = scipy.io.loadmat(r'D:\irimages\ir-short-path\1st-test\noise\1502\sbf161_img_000_1600.mat')
#
#decN = D['s'] / 4800

#plt.figure()
#plt.imshow(imageSeqRaw[:,:,0,79])
#plt.plot(mean)
#plt.colorbar()