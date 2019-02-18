# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 11:51:52 2018

@author: shihao
"""
import numpy
from matplotlib import pyplot as plt
import os, os.path
from matplotlib import animation as animation
import subprocess
import scipy.io
from pylab import *

wn_num = 1000
dnf_seq = numpy.zeros((128, 128, wn_num))
dn_dir = r'D:\irimages\irholography\newqcl\spectrum2'
dn_list = os.listdir(dn_dir)
for i in range(len(dn_list)):
    dnf_dir = dn_dir + '\\' + dn_list[i]
    dnf_data = numpy.fromfile(dnf_dir, dtype = numpy.uint16, count = -1, sep = '')
    dnf_data = numpy.reshape(dnf_data, (128, 128), order = 'C')
    dnf_seq[:, :, i] = dnf_data

raw_spectrum = dnf_seq
scipy.io.savemat(dn_dir +'\\'+'raw_spectrum.mat', {'raw_spectrum':(raw_spectrum)})

#wn_x = numpy.arange(900, 1900, 1)
#
#plt.figure()
#plt.plot(wn_x, dnf_seq[64,64,:])


