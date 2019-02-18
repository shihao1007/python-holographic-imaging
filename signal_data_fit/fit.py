import numpy
from matplotlib import pyplot as plt
import os, os.path
from matplotlib import animation as animation
from pylab import *
import scipy.io

#specify all imaging parameters
data_dir = r'D:\irimages\irholography\spectrum\HD2'  #directory of all raw data
zstep = 0.25               #stage moving stepsize, delta_m, lightpath length difference in the hologram
num_frames = 25         #number of averaged frames for each image, to minimize laser fluctuation
positions = 32          #number of different positions of the stage. i.e. number of measurements
ref_pos = 32
lp = 81                 #laser power
inte = 30               #integration time of FPA
wn_start = 1420         #starting wavenumber of the spectrum
wn_end = 1680           #ending wavenumber of the spectrum
wn_num = 130             #number of bands to be imaged
badBands = [1456, 1538, 1557, 1558, 1559, 1634, 1652]

smoothrange = numpy.arange(2, 128, 1)

def spectrumSmooth(spectrum):
    smoothed = spectrum
    for i in smoothrange:
        smoothed[i] = (spectrum[i - 2] + spectrum[i - 1] + spectrum[i] + spectrum[i + 1] + spectrum[i + 2]) / 5
    
    return smoothed

x = numpy.arange(0,130,1)
plt.figure()
#plt.imshow(numpy.abs(FSub[:,:,0]))

plt.subplot(311)
plt.plot(x, spectrumSmooth(numpy.abs(FSub[100,20,:])))
plt.subplot(312)
plt.plot(x, spectrumSmooth(numpy.abs(FSub[100,80,:])))
plt.subplot(313)
plt.plot(x, spectrumSmooth(numpy.abs(FSub[100,80,:])- numpy.abs(FSub[100,20,:])))

