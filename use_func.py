import os,sys
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
	pixscale=round(np.sqrt(cd11**2 + cd21**2) *3600 ,4)
	puthdr(i,'PSCALE',round(pixscale,3))
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
