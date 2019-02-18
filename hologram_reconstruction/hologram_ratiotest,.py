# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 14:12:59 2018

@author: shihao
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 10:20:40 2018

@author: shihao
"""
"""
Created on Thu Aug  2 13:35:03 2018

@author: shihao
"""

import numpy
from matplotlib import pyplot as plt
import os, os.path
from matplotlib import animation as animation
from pylab import *
import scipy.io

#specify all imaging parameters
data_dir = r'D:\irimages\irholography\spectrum\auto1'  #directory of all raw data
zstep = 0.25               #stage moving stepsize, delta_m, lightpath length difference in the hologram
num_frames = 25         #number of averaged frames for each image, to minimize laser fluctuation
positions = 16          #number of different positions of the stage. i.e. number of measurements
ref_pos = 16
lp = 81                 #laser power
inte = 30               #integration time of FPA
wn_start = 1420         #starting wavenumber of the spectrum
wn_end = 1660           #ending wavenumber of the spectrum
wn_num = 120             #number of bands to be imaged
realBands = numpy.loadtxt(r'D:\irimages\irholography\spectrum\realBands.txt')
powerSeq = numpy.loadtxt(r'D:\irimages\irholography\spectrum\powerSeq.txt')

##calculate the field with phase information from the hologram
def cal_field(start, num, wn, image_seq, ref_sqrt):
    #num: number of positions in the hologram
    #wn: current wavenumber of the hologram
    #image_seq: dataset of the hologram, can be raw image, background-subed image
    #ref_sqrt: square-root of reference image of current wavenumber
    
    lambdA = 10000 / wn                  #wavelength of the laser
    k = 2 * numpy.pi / lambdA   #wave vector k
    
    #D matrix denoting the difference between consecutive measuremnts
    D = [(image_seq[:,:,start+ i + 1] - image_seq[:,:,start+i]) for i in range(num)]
    
    #P matrix denoting the phase shifting terms
    P = [[numpy.exp(-1j * k * zstep * 2 * i), numpy.exp(1j * k * zstep * 2 * i)] for i in range(num)]
    
    #convert D and P to matrices or numpy arrays from lists
    D = numpy.asarray(D)
    P = numpy.matrix(numpy.asarray(P))
    
    #calculate the Hermitian conjugate of P matrix and its inverse
    P_hermint = P.getH()
    PH_P_inv = numpy.linalg.inv(P_hermint * P)
    
    #linear P matrix term before D matrix when calculating interferometric cross terms u
    PH_P_Ph =  PH_P_inv * P_hermint
    
    #initialize interferometric cross terms u as a 2x1 matrix for each point on the image with itself and its conjugate
    u = numpy.zeros((2,1,128,128)) + 0j
    
    #transform P matrix term and duplicate itself into higher dimention for dot product with D matrix for the whole image
    PH_P_Ph = numpy.asarray(PH_P_Ph)
    PH_P_Ph = numpy.repeat((numpy.repeat(numpy.reshape(PH_P_Ph,(2,num,1,1)),128,axis = 2)),128,axis = 3)
    
    #for each pixel on the image, do a doc product to calculate interferometric cross terms u
    D = D.reshape(num,1,128,128)  + 0j
    for i in range(128):
        for j in range(128):
            u[:,:,i,j] = numpy.dot(PH_P_Ph[:,:,i,j], D[:,:,i,j])
    
    #initialize a field to store u from the cross terms
    f = numpy.zeros((1,1,128,128)) + 0j
    f[:] = u[0,:,:,:]
    f = f.reshape((128,128))
    
    
    #calculate the scattering field F from reference field and extra terms in the equation
    expterm = numpy.exp(-1j * k * zstep) - 1
    dominator = expterm * ref_sqrt
    F = f / dominator
    
    return F

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
imageSeqRaw = readFile('edge')
bgSeq = readFile('bg')

#calculate the subtraction data set
imageSeqSub = numpy.zeros((128, 128, positions, wn_num))
for i in range(positions):
    for j in range(wn_num):
        imageSeqSub[:,:,i,j] = imageSeqRaw[:,:,i,j] - bgSeq[:,:,i,j]

#read in reference data
refSeq = readFile('ref')
refIntensity = numpy.mean(refSeq, axis = 2)
refSqrt = numpy.sqrt(refIntensity)

#for each wavenumber calculate the field
FRaw = numpy.zeros((128, 128, wn_num)) + 0j
FBg = numpy.zeros((128, 128, wn_num)) + 0j
F = numpy.zeros((128, 128, wn_num)) + 0j
FSub = numpy.zeros((128, 128, wn_num)) + 0j
#for i in range(wn_num):
#    index = i * int((wn_end - wn_start) / wn_num)
#    wn = realBands[index]
#    FRaw[:,:,i] = cal_field(positions - 1, wn, imageSeqRaw[:,:,:,i], refSqrt[:,:,i])
#    FBg[:,:,i] = cal_field(positions - 1, wn, bgSeq[:,:,:,i], refSqrt[:,:,i])
#    FSub[:,:,i] = numpy.abs(FRaw[:,:,i])-numpy.abs(FBg[:,:,i])
##    F[:,:,i] = -numpy.log10(numpy.abs(FRaw[:,:,i])/numpy.abs(FBg[:,:,i]))
#save the fields into a mat file
#scipy.io.savemat(r'D:\irimages\irholography\spectrum\auto1\FSub_edge', {'FSub':numpy.abs(FSub)})



#p = numpy.loadtxt(r'D:\irimages\irholography\powerseq.txt')
##
#plt.figure()
#waveNumber = 0
##positions = 5
##plt.plot(numpy.abs(F[80,40,:]))
##   
##    get the maximal and minimal values of all images to turn off auto-scaling
##_min, _max = numpy.amin(imageSeqSub[:,:,:,waveNumber]), numpy.amax(imageSeqSub[:,:,:,waveNumber])    
##_min, _max = numpy.amin(imageSeqSub[:,:,positions,:]), numpy.amax(imageSeqSub[:,:,positions])  
#_min, _max = numpy.amin(numpy.abs(FSub[:,:,:])), numpy.amax(numpy.abs(FSub[:,:,:]))  
##plot the hologram and write out a .mp4 file
#fig = plt.figure()
##for i in range(positions):
##    plt.imshow(image_seq_raw[:, :, i], vmin = _min, vmax = _max)
##    plt.pause(0.05)
##plt.show()
#img = []
#for i in range(wn_num):
#    img.append([plt.imshow(numpy.abs(FSub[:,:,i]), vmin = _min, vmax = _max)])
#    
##for i in range(wn_num):
##    img.append([plt.imshow(imageSeqSub[:,:,positions,i], vmin = _min, vmax = _max)])
#
#ani = animation.ArtistAnimation(fig,img,interval=100)
#writer = animation.writers['ffmpeg'](fps=10)
##
#ani.save(data_dir + '\\spectrum.mp4',writer=writer)
    
def plot_b4_ratio(image_seq):
    plt.figure()
    
    
    for i in range(8):
        F = cal_field(0 ,15, realBands[i*10], image_seq[:,:,:,i*10], refSqrt[:,:,i*10])
#        F_bg = cal_field(num_pri, start, bg_seq)
        plt.subplot(2,4,i+1)
        plt.imshow(numpy.abs(F))
        plt.title(str(i*10) + 'th Band')
        plt.colorbar()
        
#    plt.suptitle('Field Ratio Before Reconstruction plots from ' + str(start) + 'th of the hologram', fontsize = 15)
    
#def plot_b4_sub(image_seq):
#    plt.figure()
#    
##    plt.subplot(231)
##    plt.imshow(image_seq[:,:,start])
##    plt.title('1st Raw Image')
##    plt.colorbar()
#    6
#    for i in range(7):
#        num_pri = num * (2 ** i)
#        F_sub = cal_field(num_pri, start, image_seq)
##        F_bg = cal_field(num_pri, start, bg_seq)
#        plt.subplot(2,4,i+1)
#        plt.imshow(numpy.abs(F_sub))
#        plt.title(str(num_pri) + ' measurements')
#        plt.colorbar()
#        
#    num_last = positions - 1
#    F_sub = cal_field(num_last, start, image_seq)
##        F_bg = cal_field(num_pri, start, bg_seq)
#    plt.subplot(2,4,8)
#    plt.imshow(numpy.abs(F_sub))
#    plt.title(str(num_last) + ' measurements')
#    plt.colorbar()
#        
#    plt.suptitle('Field Sub', fontsize = 15)
    


#F_bead = cal_field(num, start, image_seq_raw)
#F_bg = cal_field(num, start, bg_seq)
#F_ratio = cal_field(num, start, image_seq_ratio)
#F_sub = cal_field(num, start, image_seq_sub)
#plot_fields(image_seq_raw, F_bead)
#plot_fields(bg_seq,F_bg)
#plot_fields(image_seq_raw - bg_seq, F_bead-F_bg, 'Field Sub')
#plot_fields(image_seq_raw - bg_seq, F_ratio, 'Ratio')
#plot_fields(image_seq_raw - bg_seq, F_sub, 'Subtraction')
#plot_fields(numpy.log10(image_seq_raw[:,:,start] / bg_seq[:,:,start]), numpy.log10(F_bead/F_bg))
#plot_nums_sub(image_seq_raw, bg_seq)
#plot_nums_ratio(image_seq_raw, bg_seq)
plot_b4_ratio(imageSeqSub)
#plot_b4_sub(image_seq_raw)
#plot_holo(image_seq_sub)