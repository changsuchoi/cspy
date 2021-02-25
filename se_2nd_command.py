# 2st photometry
# final catalog *.dat

import os
import glob
import astropy.io.fits as fits
import numpy as np
import subprocess
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clipping
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Process,Pool


# input files, config and params
seconfigdir ='/data7/cschoi/code/cspy/sex.config/'
seconfig    ='se1.sex'
separam     ='se1.param'
separam_noPSF = 'se1_noPSF.param'
seconv      ='default.conv'
sennw       ='default.nnw'
DETECT_MINAREA = str(3)
DETECT_THRESH  = str(1.5)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)
