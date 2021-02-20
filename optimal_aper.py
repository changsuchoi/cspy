#growth curve fit to estimate optimal aperture size

import os
import glob
from sys import flags
import astropy.io.fits as fits
import numpy as np
import subprocess
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clipping
import matplotlib.pyplot as plt
import seaborn as sns

def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'

def pixelscale(i):
	cd11 = fits.getheader(i)['CD1_1']
	cd12 = fits.getheader(i)['CD1_2']
	cd21 = fits.getheader(i)['CD2_1']
	cd22 = fits.getheader(i)['CD2_2']
	pixscale=round(np.sqrt(cd11**2 + cd21**2) *3600 ,4)
	puthdr(i,'PSCALE',pixscale)
	print('Pixel scale =', pixscale,'\"')
	return pixscale

# input files, config and params
seconfigdir ='/data7/cschoi/code/cspy/sex.config/'
seconfig    ='se1.sex'
separam     ='se1.param'
growthparam     ='growth.param'
separam_noPSF = 'se1_noPSF.param'
seconv      ='default.conv'
sennw       ='default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH  = str(5)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)

def secom(im,psf=False):
    #PSCALE=fits.getheader(i)['PSCALE']
	PSCALE=pixelscale(im)
	#aper_list=[3,5,7]
	#aper_input = ''
	#for i in aper_list: aper_input += '{},'.format(round(i/pixelscale(im),1))
	#aper_input = aper_input[:-1]
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.growth'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b=' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2b=' -PARAMETERS_NAME '+seconfigdir+growthparam

	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt8 + im
	print(secommand)
	sexout = subprocess.getoutput(secommand)
	line = [s for s in sexout.split('\n') if 'RMS' in s]
	skymed, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	return skymed, skysig
	os.system(secommand)

growthcat=ascii.read(fn+'.growth')    
aper = np.linspace(1,60,60)
snrwin=growthcat['SNR_WIN']
# catalog cleaning
flag
snr_win
magautoerr

# signal to noise ratio

flux/fluxerr

# mean values

# curve fitting

# first maximum 

