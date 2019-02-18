# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 14:40:56 2018

@author: shihao
"""

from scipy import signal
import matplotlib.pyplot as plt
import numpy

t = numpy.linspace(0, 1, 1000, endpoint = False)
s = numpy.random.uniform(0,1,1000)
plt.plot(t, 7+signal.square(2 * numpy.pi * 3 * t), label = 'TTL Trigger')
plt.plot(t, 4+signal.square(2 * numpy.pi * 3 * t, duty = 0.1), label = 'Pulse Trigger')
plt.plot(t, 1 + s * signal.square(2 * numpy.pi * 100 * t), label = 'Laser Pulse')
plt.plot(t, -2+signal.square(2 * numpy.pi * 3 * t, duty = 0.04), label = 'DAQ')
plt.plot(t, -5+signal.square(2 * numpy.pi * 3 * (t - 0.02), duty = 0.025), label = 'FPA')
plt.xlabel('Time')
plt.ylim(-10,10)
plt.legend()