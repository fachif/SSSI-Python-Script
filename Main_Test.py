'''
Test script for signal filter!
'''

import numpy as np
import scipy.io
import math
import matplotlib.pyplot as plt

from test_filter import test_filter
from sssi import sssi

mat = scipy.io.loadmat('7bus_fault.mat')

time = mat['ts']
signal = mat['P_line']
step_min = 0.04
var_min = 0.1
f = [0.1,2.5]
d = [0,5,10]
Nm = 10

sss_smi, sss_ami, sss_gmi, out_t, out_y, model_Poles, model_Res, model_K, model_that, model_yhat, detail_d, detail_f, detail_l = sssi(signal, time, step_min, var_min, f, d, Nm)

print("Small Signal Stability SMI index: {}" .format(sss_smi))
print("Small Signal Stability AMI index: {}" .format(sss_ami))
print("Small Signal Stability GMI index: {}" .format(sss_gmi))
