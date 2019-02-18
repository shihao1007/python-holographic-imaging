# -*- coding: utf-8 -*-
"""
Created on Sat Jul 28 14:26:28 2018

@author: shihao
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 15:10:29 2018

@author: shihao
"""

##fit the image into a sin function

import numpy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os, os.path
from matplotlib import animation as animation
import subprocess
from pylab import *
import scipy.io

def func(x, A, theta):
    return numpy.real(A * numpy.exp(x * theta * 1j))

#specify all imaging parameters
data_dir = r'D:\irimages\irholography\oldQCL\bead7\bead_holo'  #directory of all raw data
zstep = 0.1               #stage moving stepsize, delta_m, lightpath length difference in the hologram
num_frames = 50         #number of averaged frames for each image, to minimize laser fluctuation
positions = 80          #number of different positions of the stage. i.e. number of measurements
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

#read in reference image
re_dir = r'D:\irimages\irholography\oldQCL\bead7\ref'
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
bgs_dir = r'D:\irimages\irholography\oldQCL\bead7\bg_holo'
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
        image_seq_raw[:, :, i] = image_data
        image_seq_sub[:, :, i] = image_data - bg_seq[:, :, i]
        
    
#specify which portion of hologram will be used to calculate field
num = positions - 1                    #number of measurements used to do phase reconstruction. theoritically more than 3 is enough
start = 0                   #the starting point of sampling from the raw dataset
end = start + num           #the end point of sampling
lambdA = 8                  #wavelength of the laser
k = 2 * numpy.pi / lambdA   #wave vector k
Dx_index = 56
Dy_index = 43


##calculate D and P term
#D matrix denoting the difference between consecutive measuremnts
D = [(image_seq_raw[:,:,start + i + 1] - image_seq_raw[:,:,start + i]) for i in range(num)]
D_sub = [(image_seq_sub[:,:,start + i + 1] - image_seq_sub[:,:,start + i]) for i in range(num)]
#P matrix denoting the phase shifting terms
P = [[numpy.exp(-1j * k * zstep * 2 * i), numpy.exp(1j * k * zstep  * 2  * i)] for i in range(num)]

#convert D and P to matrices or numpy arrays from lists
D = numpy.asarray(D)   #D here is 199x128x128
D_sub = numpy.asarray(D_sub)
P = numpy.matrix(numpy.asarray(P))

#calculate the Hermitian conjugate of P matrix and its inverse
P_hermint = P.getH()
PH_P_inv = numpy.linalg.inv(P_hermint * P)

#linear P matrix term before D matrix when calculating interferometric cross terms u
PH_P_Ph =  PH_P_inv * P_hermint

#initialize interferometric cross terms u as a 2x1 matrix for each point on the image with itself and its conjugate
u = numpy.zeros((2,1)) + 0j

#for each pixel on the image, do a doc product to calculate interferometric cross terms u
D_pixel_A = D[:, 56, 43]
D_pixel_B = D[:, 100, 100]
#D_sub_pixel = D_sub[:, Dx_index, Dy_index]
PH_P_Ph_a = PH_P_Ph[0,:]
PH_P_Ph_a = PH_P_Ph_a.reshape((num,1))

#scipy.io.savemat(r'D:\irimages\irholography\bead7\D_pixel.mat', {'D_pixel':D_pixel})


#u = numpy.dot(PH_P_Ph, D_pixel)
#
##initialize a field to store u from the cross terms
#f = numpy.zeros((1,1,128,128)) + 0j
#f[:] = u[0,:,:,:]
#f = f.reshape((128,128))
#
#
##calculate the scattering field F from reference field and extra terms in the equation
#expterm = numpy.exp(-1j * k * zstep) - 1
#dominator = expterm * ref_sqrt
#F = f / dominator
#



x = numpy.linspace(0,num - 1, num)
xx = numpy.linspace(0, (num - 1)/10, num)
def fit_(ydata):
    ##fit in D to a sin wave

    ini = [80, 0.16]
    
    popt, pcov = curve_fit(func, x, ydata, p0=ini)
    
    return popt


#Amp = numpy.zeros((128,128))
#for i in range(128):
#    for j in range(128):
#
#        ydata = image_seq_sub[i,j,:]
#        popt, pcov = curve_fit(func, xdata, ydata)
#        Amp[i,j] = popt[0]

popt_D_A = fit_(D_pixel_A)
popt_D_B = fit_(D_pixel_B)

plt.figure()
plt.subplot(211)
plt.plot(xx, D_pixel_A, 'bo',ms=3, label='Samples at (56, 43)')
plt.plot(xx, func(x, *popt_D_A), 'k-', label="Fit")
plt.ylabel('Pixel Intensity')
plt.legend()

plt.subplot(212)
plt.plot(xx, D_pixel_B, 'go',ms=3, label='Samples at (100, 100)')
plt.plot(xx, func(x, *popt_D_B), 'c-', label="Fit")
plt.xlabel('Mirror Positions (um)')
plt.ylabel('Pixel Intensity')
plt.legend()

plt.suptitle('Curve Fits', fontsize = 15)
plt.show()





