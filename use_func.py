)import os,sys
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import numpy as np


# useful functions
# import and declare these fuctions first


def pixelscale(i):
	cd11 = fits.getheader(i)['CD1_1']
	cd12 = fits.getheader(i)['CD1_2']
	cd21 = fits.getheader(i)['CD2_1']
	cd22 = fits.getheader(i)['CD2_2']
	pixscale=round(np.sqrt(cd11**2 + cd21**2) *3600 ,5)
	puthdr(i,'PSCALE',round(pixscale,5))
	#print('Pixel scale =', pixscale,'\"')
	return pixscale

def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'

def limitmag(N, zp, aper, skysigma):			# 3? 5?, zp, diameter [pixel], skysigma
	import numpy as np
	R           = float(aper)/2.				# to radius
	braket      = N*skysigma*np.sqrt(np.pi*(R**2))
	upperlimit  = float(zp)-2.5*np.log10(braket)
	return round(upperlimit, 3)

def radec_center(im):
	from astropy.wcs import WCS
	from astropy.coordinates import SkyCoord
	from astropy.coordinates import ICRS, Galactic, FK4, FK5
	from astropy.coordinates import Angle, Latitude, Longitude
	from astropy.io import fits
	import astropy.units as u
	import astropy.coordinates as coord
	import numpy as np
	hdr = fits.getheader(im)
	#	RA, Dec center for reference catalog query
	xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
	w = WCS(im)
	racent, deccent = w.all_pix2world(xcent, ycent, 1)
	c=SkyCoord(racent,deccent,unit="deg")
	rastr=c.ra.to_string(unit=u.hourangle,sep=':')
	decstr=c.dec.to_string(unit=u.deg,sep=':')
	racent, deccent = racent.item(), deccent.item()
	return rastr,decstr,racent,deccent

def imcenter_offset(im, ra=161.645641, dec=13.750859,seplimit=10):
	import numpy as np
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	a,b,c,d=radec_center(im)
	c1=SkyCoord(c*u.deg,d*u.deg)
	c2=SkyCoord(ra*u.deg, dec*u.deg)
	sep = c1.separation(c2)
	#print(sep.arcmin, 'away from given position to image center')
	if sep.arcmin > seplimit :
		print(im,'Too far away',sep.arcmin)
	return sep



def targetfind(tra, tdec, refra, refdec, sep=2.0):
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	targ_coord	= SkyCoord(tra, tdec, unit=(u.deg, u.deg))
	phot_coord	= SkyCoord(refra, refdec, unit=(u.deg, u.deg))
	indx, d2d, d3d	= targ_coord.match_to_catalog_sky(phot_coord)
	return indx.item(), round(d2d.arcsec.item(),2)
	#return indx, d2d, d3d

def trim(inim, position, size, outim='trim.fits'):
	'''
	position=(ra,dec) in deg unit
	size=(px,py) in pixel unit
	'''
	from astropy.nddata import Cutout2D
	from astropy import units as u
	import numpy as np
	# Load the image and the WCS
	hdu = fits.open(inim)[0]
	wcs = WCS(hdu.header)
	# Make the cutout, including the WCS
	cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs, mode='partial',fill_value=1.0e-30)
	# Put the cutout image in the FITS HDU
	hdu.data = cutout.data
	# Update the FITS header with the cutout WCS
	hdu.header.update(cutout.wcs.to_header())
	# Write the cutout to a new FITS file
	hdu.writeto(outim, overwrite=True)




# running time check

import time
starttime=time.time()
# job running
endtime=time.time()
duration= endtime-starttime


starttime=time.time();
ls saCalib* | wc -l
endtime=time.time();os.system('ls saCalib* | wc -l')




# single image check from lines (image-set-group.py)
for ii in lines:
    iii=ii[:-1].split(',')
    if len(iii)==1:
        print (iii)

# stdout
'''
f = io.StringIO()
with redirect_stdout(f):
	func(arg)
s = f.getvalue()
'''


#header keyword remove


#h.remove('PV1_0')

kwds=['PV1_0', 'PV1_1', 'PV1_2', 'PV1_3', 'PV1_4', 'PV1_5',
		'PV1_6', 'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10',
		'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3', 'PV2_4', 'PV2_5',
		'PV2_6', 'PV2_7', 'PV2_8', 'PV2_9', 'PV2_10']

kwds=['PV1_5', 'PV1_6', 'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10']
def delhdrkey(im,kwds=kwds):
	h=fits.getheader(im)
	for k in kwds:
		if k in list(h.keys()):
			print('removing',k,h[k])
			fits.delval(im,k)


badlist=[]
for a,b zip(imlist,result):
	if b==None:
		print(a,b)
		badlist.append(a)
print badlist



# map
for result in map(func, arg0)
	print result

result=map(func,args)
print(list(result))
# lambda


from itertools import repeat
