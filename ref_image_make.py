# reference image making
# zp, fwhm, skyval, skysig, num_star, UL5
# based on aperture F15

import astropy.io.fits as fits
import astropy.io.ascii as ascii
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys

# after se1st()
# ZP, fwhm_pix, skyval,
os.system('rm *ap.fits *seg.fits')
os.system('gethead saCalib*com.fits DATE-OBS ZP_F15 UL5_F15 FWHM_PIX SKYVAL > refcand.info')

refcand=ascii.read('refcand.info',names=['FILE','DATE','ZP','UL','FWHM','SKYVAL'])

refcand.sort('FILE')

def ds9list(imlist):
	ll=''
	for l in imlist: ll+=l +' '
	os.system('ds9 '+ ll+' &')

def swarp_list(ii):
	''' input list'''
	if len(ii)==1 :
		print("single image, change name only")
		if os.path.isfile(centertimeheader(ii[0])[1]) : pass
		else:
			os.system('cp '+ii[0]+'s '+centertimeheader(ii[0])[1])
			#puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])
	else:
		if os.path.isfile(centertimeheader(ii)[1]) : pass
		else:
			print(len(ii),'images will be combined', centertimeheader(ii)[1])
			swarpcom(ii)
			puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])
