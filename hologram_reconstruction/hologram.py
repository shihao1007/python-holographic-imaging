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

import numpy
from matplotlib import pyplot as plt
import os, os.path
from matplotlib import animation as animation
import subprocess
from pylab import *

#specify all imaging parameters
data_dir = r'D:\irimages\irholography\oldQCL\bead8\bead_holo'  #directory of all raw data
zstep = 1               #stage moving stepsize, delta_m, lightpath length difference in the hologram
num_frames = 50         #number of averaged frames for each image, to minimize laser fluctuation
positions = 240          #number of different positions of the stage. i.e. number of measurements
ref_pos = 20
lp = 60                 #laser power
inte = 33               #integration time of FPA

#initialize a array for storing all image data
image_seq = numpy.zeros((128, 128, positions))
image_seq_raw = numpy.zeros((128, 128, positions))
image_seq_ratio = numpy.zeros((128, 128, positions))
image_seq_sub = numpy.zeros((128, 128, positions))
ref_seq = numpy.zeros((128, 128, positions))
ref_sqrt = numpy.zeros((128, 128))
bg_seq = numpy.zeros((128, 128, positions))
bg_obj_seq = numpy.zeros((128, 128, ref_pos))
bead_obj_seq= numpy.zeros((128, 128, ref_pos))
dnf_seq = numpy.zeros((128, 128, 5))

#call .exe file to do imaging
#subprocess.call([r"qcl-holo-cap.exe", data_dir, '--zstep', str(zstep), '--positions', str(positions), '--frames', str(num_frames), '--laserpower', str(lp), '--inte', str(inte)])

#read in detector noise
dn_dir = r'D:\irimages\irholography\oldQCL\detector_noise'
dn_list = os.listdir(dn_dir)
number_dn = len(dn_list)
for i in range(number_dn):
    dnf_dir = dn_dir + '\\' + dn_list[i]
    dnf_data = numpy.fromfile(dnf_dir, dtype = numpy.uint16, count = -1, sep = '')
    dnf_data = numpy.reshape(dnf_data, (128, 128), order = 'C')
    dnf_seq[:, :, i] = dnf_data
    
dnf_intensity = sum(dnf_seq, axis = 2) / number_dn

#read in reference image
re_dir = r'D:\irimages\irholography\oldQCL\bead8\ref'
re_list = os.listdir(re_dir)
number_res = len(re_list)
if( number_res != ref_pos):
    print("Error: Reference Image is missing!")
else:
    for i in range(number_res):
        ref_dir = re_dir + '\\' + re_list[i]
        ref_data = numpy.fromfile(ref_dir, dtype = numpy.uint16, count = -1, sep = '')
        ref_data = numpy.reshape(ref_data, (128, 128), order = 'C')
        ref_seq[:, :, i] = ref_data
        
    ref_intensity = sum(ref_seq, axis = 2) / number_res
    ref_sqrt = numpy.sqrt(ref_intensity)
    
#read in background image
bgs_dir = r'D:\irimages\irholography\oldQCL\bead8\bg_holo'
bgs_list = os.listdir(bgs_dir)
number_bgs = len(bgs_list)
if( number_bgs != positions):
    print("Error: Background Image is missing!")
else:
    for i in range(number_bgs):
        bg_dir = bgs_dir + '\\' + bgs_list[i]
        bg_data = numpy.fromfile(bg_dir, dtype = numpy.uint16, count = -1, sep = '')
        bg_data = numpy.reshape(bg_data, (128, 128), order = 'C')
        bg_seq[:, :, i] = bg_data
        
    bg_average = sum(bg_seq, axis = 2) / positions
    bg_normalize = bg_average 

#reading in all images to the data array initialized before
image_list = os.listdir(data_dir)
number_files = len(image_list)
if( number_files != positions):
    print("Error: Image is missing!")
else:
    for i in range(number_files):
        image_dir = data_dir + '\\' + image_list[i]
        image_data = numpy.fromfile(image_dir, dtype = numpy.uint16, count = -1, sep = '')
        image_data = numpy.reshape(image_data, (128, 128), order = 'C')
        image_data = image_data - dnf_intensity
        image_seq_raw[:, :, i] = image_data
#        image_seq_ratio[:, :, i] = - numpy.log10(image_data / bg_seq[:, :, i])
        image_seq_sub[:, :, i] = image_data - bg_seq[:, :, i]

#raw = False
#ratio =False
#sub = True
#
#if(raw == True): image_seq = image_seq_raw
#elif(raw == False and ratio == True): image_seq == image_seq_ratio
#elif(raw == False and ratio == False and sub == True): image_seq == image_seq_sub


##get the maximal and minimal values of all images to turn off auto-scaling
#_min, _max = numpy.amin(image_seq_raw), numpy.amax(image_seq_raw)    
#
#
##plot the hologram and write out a .mp4 file
#fig = plt.figure()
##for i in range(positions):
##    plt.imshow(image_seq_raw[:, :, i], vmin = _min, vmax = _max)
##    plt.pause(0.05)
##plt.show()
#img = []
#for i in range(positions):
#    img.append([plt.imshow(image_seq_raw[:, :, i], vmin = _min, vmax = _max)])
#
#ani = animation.ArtistAnimation(fig,img,interval=200)
#writer = animation.writers['ffmpeg'](fps=5)

#ani.save(data_dir + '\\demo_' + str(num_frames) + '.mp4',writer=writer)


#specify which portion of hologram will be used to calculate field
num = 6                    #number of measurements used to do phase reconstruction. theoritically more than 3 is enough
start = 0                   #the starting point of sampling from the raw dataset
end = start + num           #the end point of sampling
lambdA = 8                  #wavelength of the laser
k = 2 * numpy.pi / lambdA   #wave vector k



#D matrix denoting the difference between consecutive measuremnts
D = [(image_seq_sub[:,:,start + i + 1] - image_seq_sub[:,:,start + i]) for i in range(num)]

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
#PH_P_Ph_dot = numpy.dot(PH_P_inv, P_hermint)



#initialize interferometric cross terms u as a 2x1 matrix for each point on the image with itself and its conjugate
u = numpy.zeros((2,1,128,128)) + 0j
#u = numpy.zeros((2,1)) + 0j


#transform P matrix term and duplicate itself into higher dimention for dot product with D matrix for the whole image
PH_P_Ph = numpy.asarray(PH_P_Ph)
PH_P_Ph = numpy.repeat((numpy.repeat(numpy.reshape(PH_P_Ph,(2,num,1,1)),128,axis = 2)),128,axis = 3)

#for each pixel on the image, do a doc product to calculate interferometric cross terms u
D = D.reshape(num,1,128,128)  + 0j
for i in range(128):
    for j in range(128):
        u[:,:,i,j] = numpy.dot(PH_P_Ph[:,:,i,j], D[:,:,i,j])

#u = numpy.dot(PH_P_Ph, D)

f = numpy.zeros((1,1,128,128)) + 0j
f[:] = u[0,:,:,:]
f = f.reshape((128,128))

#calculate the scattering field F from reference field and extra terms in the equation
expterm = numpy.exp(-1j * k * zstep) - 1
dominator = expterm * ref_sqrt
F = f / dominator

#plot out scattering field
plt.figure()

plt.subplot(221)
plt.imshow(numpy.abs(dominator))
plt.title('dominator')
plt.colorbar()

plt.subplot(222)
plt.imshow(numpy.abs(F))
plt.title('Absolute Value')
plt.colorbar()

plt.subplot(223)
plt.imshow(numpy.real(F))
plt.title('Real Part')
plt.colorbar()

plt.subplot(224)
plt.imshow(numpy.imag(F))
plt.title('Imaginary Part')
plt.colorbar()

plt.suptitle('Bead7')