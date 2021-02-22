# growth curve fit to estimate optimal aperture size

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
from scipy.interpolate import UnivariateSpline

def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr = fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'


def pixelscale(i):
	cd11 = fits.getheader(i)['CD1_1']
	cd12 = fits.getheader(i)['CD1_2']
	cd21 = fits.getheader(i)['CD2_1']
	cd22 = fits.getheader(i)['CD2_2']
	pixscale = round(np.sqrt(cd11**2 + cd21**2) * 3600, 4)
	puthdr(i, 'PSCALE', round(pixscale, 3))
	# print('Pixel scale =', pixscale,'\"')
	return pixscale


# input files, config and params
seconfigdir = '/data7/cschoi/code/cspy/sex.config/'
seconfig = 'se1.sex'
separam = 'se1.param'
growthparam = 'growth.param'
separam_noPSF = 'se1_noPSF.param'
seconv = 'default.conv'
sennw = 'default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)


def segrowthcom(im, psf=False):
	PSCALE = pixelscale(im)
	fn = os.path.splitext(im)[0]
	aper_list=[s for s in range(1,51)]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(i,1)
	aper_input = aper_input[:-1]
	opt1 = seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME ' + fn+'.growth'
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	#opt2a = ' -PARAMETERS_NAME '+seconfigdir+separam
	#opt2b = ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2b = ' -PARAMETERS_NAME '+seconfigdir+growthparam
	opt2 = ' -FILTER_NAME '+seconfigdir+seconv + ' -STARNNW_NAME '+seconfigdir+sennw
	opt3 = ' -DETECT_MINAREA ' + DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4 = ' -DEBLEND_NTHRESH ' + DEBLEND_NTHRESH + \
	    ' -DEBLEND_MINCONT ' + DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE NONE '
	opt6 = ' -PHOT_APERTURES '+aper_input+' '
	opt7 = ' -PSF_NAME '+fn+'.psf '
	opt8 = ' -PIXEL_SCALE '+str(PSCALE)+' '
	secommand = 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8 + im
	print(secommand)
	sexout = subprocess.getoutput(secommand)
	line = [s for s in sexout.split('\n') if 'RMS' in s]
	skymed = float(line[0].split('Background:')[1].split('RMS:')[0])
	skysig= float(line[0].split('RMS:')[1].split('/')[0])
	os.system(secommand)
	return skymed, skysig

def fwhm_img(im,mtbl1):
	fwhm_img=sigma_clipped_stats(mtbl1['FWHM_IMAGE'], sigma=3, maxiters=10)
	filtered_data=sigma_clip(mtbl1['FWHM_IMAGE'],sigma=3,maxiters=10)
	selected, nonselected= ~filtered_data.mask, filtered_data.mask
	print('FWHM_IMAGE','{}'.format(round(fwhm_img[0],3)),
		len(mtbl1[selected]),'stars from',len(mtbl1))
	puthdr(im, 'FWHM_PIX', round(fwhm_img[0],3), hdrcomment='FWHM PIXEL VALUE')
	return round(fwhm_img[0],3)

def find_opt_aper(im):
	skymed, skysig=segrowthcom(im)
	fn=os.path.splitext(im)[0]
	growthcat=ascii.read(fn+'.growth')
	aper = np.linspace(1,50,50)
	# catalog cleaning
	idx= np.where( (growthcat['FLAGS']==0) &
            (growthcat['SNR_WIN']>100) &
            (growthcat['MAGERR_AUTO']<0.1)
        	)
	mcat=growthcat[idx]
	# signal to noise ratio
	fwhm=fwhm_img(im,mcat)
	plt.plot(np.ones(len(mcat))*aper[0],mcat['FLUX_APER']/mcat['FLUXERR_APER'],'k.')
	for q in aper[:-1]:
		p=int(q)
		plt.plot(np.ones(len(mcat))*aper[p],mcat['FLUX_APER_'+str(p)]/mcat['FLUXERR_APER_'+str(p)],'k.')
	# mean values
	plt.plot(aper[0],np.mean(mcat['FLUX_APER']/mcat['FLUXERR_APER']),'bo')
	for q in aper[:-1]:
		p=int(q)
		plt.plot(aper[p],np.mean(mcat['FLUX_APER_'+str(p)]/mcat['FLUXERR_APER_'+str(p)]),'bo')
	# plt.plot(np.ones(len(mcat))*0,mcat['FLUX_PSF']/mcat['FLUXERR_PSF'],'ro')
	# mcat['FLUX_AUTO']/mcat['FLUXERR_AUTO']
	plt.vlines(fwhm*1.5, 0,10000)
	# curve fitting
	snrval=[0]
	snrval[0]=np.mean(mcat['FLUX_APER']/mcat['FLUXERR_APER'])
	for q in aper[:-1]:
		p=int(q)
		snrval.append(np.mean(mcat['FLUX_APER_'+str(p)]/mcat['FLUXERR_APER_'+str(p)]))
	sarr=sarr=np.asarray(snrval)
	sarr1=sarr[~np.isnan(sarr)]
	aper1=aper[~np.isnan(sarr)]
	s = UnivariateSpline(aper1, sarr1, s=1)
	x=np.linspace(1,50,500)
	ys = s(x)
	idx_opt=np.where(ys==np.max(ys))
	plt.plot(x,ys)
	plt.vlines(x[idx_opt],0,10000)
	plt.title(fn+' '+'optimal ap size')
	plt.xlabel('aperture size (pixel)')
	plt.ylabel('SNR (FLUX / FLUX err)')
	plt.savefig(fn+'_growth.png')
	plt.close()

	print('Optimal aperture',round(x[idx_opt][0],2))
	return skymed, skysig, fwhm, round(x[idx_opt][0],2)
# first maximum
def opt_ap_fwhm(im):
	skyval, skysig,fwhm,opt_ap=find_opt_aper(im)
	puthdr(im, 'SKYVAL', skymed,
		hdrcomment='sky median value form sextractor')
	puthdr(im, 'SKYSIG', skysig,
		hdrcomment='sky sigma value form sextractor')
	puthdr(im, 'FWHM_PIX', fwhm, hdrcomment='FWHM PIXEL VALUE')
	puthdr(im, 'OPT_AP', opt_ap, hdrcomment='Optimal aperture size')
	return skyval, skysig,fwhm,opt_ap
