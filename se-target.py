#source extractor on subtracted images

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
import astropy.io.ascii as ascii
from scipy.interpolate import UnivariateSpline
from multiprocessing import Process,Pool


seconfigdir ='/data7/cschoi/code/cspy/sex.config/'
seconfig    ='se1.sex'
separam     ='se1.param'
separam_noPSF = 'se1_noPSF.param'
growthparam = 'growth.param'
seconv      ='default.conv'
sennw       ='default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH  = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)
lowmag=13
highmag=19
filname,filerr='R','Rerr'
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
magtype=magtypes[0]


# isolated field object, no need to subtraciotn
psf=False
def se2com(im,psf=psf):
	PSCALE=fits.getheader(im)['PSCALE']
	fwhm_pix=fits.getheader(im)['FWHM_PIX']
	opt_ap=fits.getheader(im)['OPT_AP']
	aper_list=[3,5,7]
	aper_list_fwhm=[1.0,1.5,2.0,2.5,3.0]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	#aper_input = aper_input[:-1]
	for i in aper_list_fwhm: aper_input += '{},'.format(round(i*fwhm_pix,1))
	aper_input = aper_input+str(opt_ap)
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.sef'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	opt5a=' -CHECKIMAGE_TYPE NONE '
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	opt9=' -SEEING_FWHM '+str(round(fwhm_pix*PSCALE))+ ' '
	if psf==True:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 + im
	else:
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 + im
	os.system(secommand)


# for subtracted images
psf=True
# se_sub(im,psf=psf)
def se_sub(im,psf=psf,dual=False,det=False):
	PSCALE=fits.getheader(im)['PSCALE']
	fwhm_pix=fits.getheader(im)['FWHM_PIX']
	opt_ap=fits.getheader(im)['OPT_AP']
	aper_list=[3,5,7]
	aper_list_fwhm=[1.0,1.5,2.0,2.5,3.0]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	#aper_input = aper_input[:-1]
	for i in aper_list_fwhm: aper_input += '{},'.format(round(i*fwhm_pix,1))
	aper_input = aper_input+str(opt_ap)
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.sub'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	opt5a=' -CHECKIMAGE_TYPE NONE '
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	opt9=' -SEEING_FWHM '+str(round(fwhm_pix*PSCALE))+ ' '
	if psf==True & dual==True:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 +det+' '+ im
	elif psf==True & dual==False:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 + im
	elif psf==False & dual==True:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt8+opt9 +det+' '+ im
	else:
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 + im
	os.system(secommand)


tra,tdec=[]
def targetfind(tra, tdec, refra, refdec, sep=2.0):
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	targ_coord	= SkyCoord(tra, tdec, unit=(u.deg, u.deg))
	phot_coord	= SkyCoord(refra, refdec, unit=(u.deg, u.deg))
	indx, d2d, d3d	= targ_coord.match_to_catalog_sky(phot_coord)
	return indx.item(), round(d2d.arcsec.item(),2)
