from photutils.segmentation import make_source_mask
import astropy.io.fits as fits
from astropy.stats import sigma_clipped_stats
import numpy as np
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground
import os
import parmap


sigma_clip = SigmaClip(sigma=3.0)
bkg = SExtractorBackground(sigma_clip)

lines= oswalkfunc()
lines.sort()

fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']

mao=[i for i in fitslist if 'MAO' in i]
snucam=[i for i in mao if 'SNUCAM' in i]
snucamI=[i for i in snucam if '-I-' in i]
snucamI60=[i for i in snucamI if '-I-60' in i]

def ds9list(imlist):
	ll=''
	for l in imlist: ll+=l +' '
	os.system('ds9 '+ ll+' &')

for i in snucamI60 :
     os.system('rename MAIDANAK MAO_SNUCAM '+i)
     print(i)
for i in Ilist:
     os.system('cp '+i+' .')

# remove bad images from list : Ilist (510 ->359)
# MAO_I_fringe dircetory
# /data7/cschoi/IMSNGgalaxies/MAO_I_fringe
# find background value, mask
# normalize masked image

import glob
Ilist=glob.glob('Calib*.fits')
Ilist.sort()

def mask_bg_norm(im):
	data=fits.getdata(im)
	fn=os.path.splitext(im)[0]
	bkg_value = bkg.calc_background(data)
	normdata=data/bkg_value
	mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=15)
	#mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
	#mean, median, std = sigma_clipped_stats(data, sigma=3.0)
	masked_data = ~mask * normdata 
	fits.writeto( fn+'_mask.fits', masked_data,overwrite=True)

for n,im in enumerate(Ilist):
	print(n+1, len(Ilist), im)
	mask_bg_norm(im)

Ilist=[i.split('/')[-1] for i in Ilist]

cpunum=5
result=parmap.map(mask_bg_norm, Ilist, pm_pbar=True, pm_processes=cpunum)

masklist=glob.glob('Calib*mask.fits')

from pyraf import iraf

def imcombine(group,output):
	if type(group) ==list : group=(",".join(group)) # if group datatype is list
	#else : print(group.split(','))
	iraf.imcombine(group,output=output,combine='average',
					project='no',reject="none",scale="none",zero="mode")

imcombine(masklist,'MAO_I_master_norm_fringe,fits')



