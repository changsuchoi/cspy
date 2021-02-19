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


im=sys.argv[1]
infile=sorted(glob.glob(im))
#size=sys.argv[2] # arcmin
# input
#im=  'test.fits'
ra, dec = 308.71805, 60.15368 #NGC6946
ra, dec = 178.20602, 44.12025
ra, dec = 161.64562673, 13.75085444
#ra, dec = 14.95871, -07.57796
#197.44879, -23.38383, #'13:09:47.71', '-23:23:01.8' # NGC4993 center coordinates J2000 deg unit
#ra, dec = 197.47311, -23.368377
size = 120 # arcmin unit, length of square side

def imcopy(inim,size):
	hdr= fits.getheader(inim)
	w = WCS(hdr)
	px, py = w.wcs_world2pix(ra, dec, 1)
	print ('center pixel coordinates', int(px), int(py) )
	xpscale,ypscale=wcsutils.proj_plane_pixel_scales(w)*60 # pixel scale armin unit
	pixscale=(xpscale+ypscale)/2.

	ax,bx=px-size/2/pixscale,px+size/2/pixscale
	ay,by=py-size/2/pixscale,py+size/2/pixscale
	print ('pixel scale =', '%.3f'% (pixscale*60), size, 'arcmin rectangular cut =',int(bx - ax),'pixels')
	region='['+str(int(ax))+':'+str(int(bx))+','+str(int(ay))+':'+str(int(by))+']'
	outname=inim[:-5]+'-'+str(size)+'-arcmin-cut.fits'
	print (outname,'will be created')
	#region='[200:2048,60:2048]'
	chinim=inim+region
	iraf.imcopy(chinim,output=outname)

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


imcopy(im,size)


for i in infile : imcopy(i,size)

print 'done.\n'
