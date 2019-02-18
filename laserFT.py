# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 10:15:00 2018

@author: shihao
"""

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

dataname = r'D:\irimages\irholography\oldQCL\laser-labview\test3\100k\laseroff1\laseroff1.txt'
a = np.zeros(250000)
a = np.loadtxt(dataname)
sampling_rate = 250000
b = np.arange(0,1,0.000004)
n = a.size
timestep = 0.000004
freq = np.fft.fftfreq(n, d=timestep)
sp = np.fft.fft(a)
T = 1/sampling_rate # inverse of the sampling rate
x = np.linspace(0.0, 1.0/(2.0*T), int(n/2))
y = 2/n * np.abs(sp[0:np.int(n/2)])

plt.figure()
plt.plot(b, -a)
plt.xlabel('Time (s)')
plt.ylabel('Amplitede')
plt.title('Time Domain (Laser Off)');
plt.grid(True)
plt.show()


spec = np.linspace(900, 1900, 500)
plt.figure()
plt.plot(x, y)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitede')
plt.title('Frequency Domain (Laser Off)');
plt.grid(True)
plt.show()


plt.figure()
plt.imshow(fileSeq[..., 100])
plt.colorbar()

spec = np.linspace(900, 1900, 500)
matplotlib.rcParams.update({'font.size': 18})
plt.figure()
plt.plot(spec, fileSeq[60, 60, :], 'r', label = 'Pixel A (60, 60)')
plt.plot(spec, fileSeq[100, 40, :], 'y', label = 'Pixel B (100, 40)')
plt.xlabel('Wavenumber')
plt.ylabel('Pixel Value')
plt.title('Laser Emission Spectrum')
plt.grid(True)
plt.legend()
plt.show()

dataname = r'D:\irimages\irholography\MIRcatFix_v4\90kHz\2sec.txt'
b = np.arange(0,2,0.000004)
samples = np.loadtxt(dataname)
plt.figure()
plt.plot(b, samples)
plt.xlabel('Time (s)')
plt.ylabel('Amplitede')
plt.title('Laser Pulse Rate 90kHz')