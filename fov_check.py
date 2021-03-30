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
		print(im,'Too far away',round(sep.arcmin,0))

ra,dec=178.206042,   44.120722 # NGC3938

for im in caliblist : imcenter_offset(im,ra=ra, dec=dec, seplimit=10)
