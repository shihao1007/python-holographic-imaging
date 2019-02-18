# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import subprocess

ncap = str(50);                         #number of captured images passed into snap.exe 
inte_fold = 10;                         #effective integration time = 15 microsecond * inte_fold
for i in range(1, inte_fold):
    nframe = str(15 * i)
    subprocess.call([r"qcl-snap.exe",r"C:\Users\shihao\Desktop\ir-images\ir-holography\exe-test",
                     '--inte_time',nframe,'--ncap',ncap,'--WN','1230'])