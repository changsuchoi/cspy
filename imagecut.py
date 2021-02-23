# code for image cut and copy from original image to small images
# in rectangluar region of input size
# and input center position (RA,DEC) ind degree unit
# python imagecut.py test.fits
# if you want to cut image in pixel coordinates, then use 'image' NOT 'physical' in DS9 window.

# test.fits -> test-cut.fits

# Changsu Choi 2017/8/22

import astropy.io.fits as fits
import numpy as np
from astropy.wcs import WCS
from pyraf import iraf
import astropy.wcs.utils as wcsutils
import os,sys
import glob


#im=sys.argv[1]
#infile=sorted(glob.glob(im))
#size=sys.argv[2] # arcmin
# input
#im=  'test.fits'
#ra, dec = 308.71805, 60.15368 #NGC6946
#ra, dec = 178.20602, 44.12025
#ra, dec = 161.64562673, 13.75085444
#ra, dec = 14.95871, -07.57796
#197.44879, -23.38383, #'13:09:47.71', '-23:23:01.8' # NGC4993 center coordinates J2000 deg unit
#ra, dec = 197.47311, -23.368377
size = 10 # arcmin unit, length of square side
ra,dec=161.645641, 13.750859
def imcopy(im,size=10,ra=False,dec=False):
	hdr= fits.getheader(im)
	w = WCS(hdr)
	xpscale,ypscale=wcsutils.proj_plane_pixel_scales(w)*60
	pixscale=(xpscale+ypscale)/2.
	if (ra==False) or (dec==False) :
		print('RA or DEC input, False, position will be center of',im)
		px, py= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
		ax,bx=px-size/2/pixscale,px+size/2/pixscale
		ay,by=py-size/2/pixscale,py+size/2/pixscale
	else:
		px, py = w.wcs_world2pix(ra, dec, 1)
		print ('center pixel coordinates', int(px), int(py) )
		ax,bx=px-size/2/pixscale,px+size/2/pixscale
		ay,by=py-size/2/pixscale,py+size/2/pixscale
	print ('pixel scale =', '%.3f'% (pixscale*60), size,
			'arcmin rectangular cut =',int(bx - ax),'pixels')
	region='['+str(int(ax))+':'+str(int(bx))+','+str(int(ay))+':'+str(int(by))+']'
	outname=os.path.splitext(im)[0]+'_'+str(size)+'min_cut.fits'
	print (outname,'will be created')
	#region='[200:2048,60:2048]'
	chinim=im+region
	iraf.imcopy(chinim,output=outname)

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

def xyimcopy(inim,sizex,sizey):
	hdr= fits.getheader(inim)
	w = WCS(hdr)
	px, py = w.wcs_world2pix(ra, dec, 1)
	print ('center pixel coordinates', int(px), int(py) )
	xpscale,ypscale=wcsutils.proj_plane_pixel_scales(w)*60 # pixel scale armin unit
	pixscale=(xpscale+ypscale)/2.

	ax,bx=px-sizex/2/pixscale,px+sizex/2/pixscale
	ay,by=py-sizey/2/pixscale,py+sizey/2/pixscale
	print ('pixel scale =', '%.3f'% (pixscale*60), sizex, 'arcmin rectangular cut =',int(bx - ax),'pixels')
	print ('pixel scale =', '%.3f'% (pixscale*60), sizey, 'arcmin rectangular cut =',int(by - ay),'pixels')
	region='['+str(int(ax))+':'+str(int(bx))+','+str(int(ay))+':'+str(int(by))+']'
	outname=inim[:-5]+'-'+str(sizex)+'+'+str(sizey)+'-arcmin-cut.fits'
	print (outname,'will be created')
	#region='[200:2048,60:2048]'
	chinim=inim+region
	iraf.imcopy(chinim,output=outname)

def imcopypix(inim,region):
	#hdr= fits.getheader(inim)
	#w = WCS(hdr)
	#px, py = w.wcs_world2pix(ra, dec, 1)
	#print 'center pixel coordinates', int(px), int(py)
	#xpscale,ypscale=wcsutils.proj_plane_pixel_scales(w)*60 # pixel scale armin unit
	#pixscale=(xpscale+ypscale)/2.

	#ax,bx=px-size/2/pixscale,px+size/2/pixscale
	#ay,by=py-size/2/pixscale,py+size/2/pixscale
	#print 'pixel scale =', '%.3f'% (pixscale*60), size, 'arcmin rectangular cut =',int(bx - ax),'pixels'
	#region='['+str(int(ax))+':'+str(int(bx))+','+str(int(ay))+':'+str(int(by))+']'
	outname=inim[:-5]+'-cut.fits'
	print outname,'will be created'
	#region='[200:2048,60:2048]'
	chinim=inim+region
	iraf.imcopy(chinim,output=outname)



#2D cutout WCS example
'''
>>> from astropy.coordinates import SkyCoord
>>> from astropy.wcs import WCS
>>> position = SkyCoord('13h11m29.96s -01d19m18.7s', frame='icrs')
>>> wcs = WCS(naxis=2)
>>> rho = np.pi / 3.
>>> scale = 0.05 / 3600.
>>> wcs.wcs.cd = [[scale*np.cos(rho), -scale*np.sin(rho)],
...               [scale*np.sin(rho), scale*np.cos(rho)]]
>>> wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
>>> wcs.wcs.crval = [position.ra.to_value(u.deg),
...                  position.dec.to_value(u.deg)]
>>> wcs.wcs.crpix = [50, 100]

# Download an example FITS file, create a 2D cutout, and save it to a
# new FITS file, including the updated cutout WCS.
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS


def download_image_save_cutout(url, position, size):
    # Download the image
    filename = download_file(url)

    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    # Make the cutout, including the WCS
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    cutout_filename = 'example_cutout.fits'
    hdu.writeto(cutout_filename, overwrite=True)


if __name__ == '__main__':
    url = 'https://astropy.stsci.edu/data/photometry/spitzer_example_image.fits'

    position = (500, 300)
    size = (400, 400)
    download_image_save_cutout(url, position, size)
'''
