#python scamp astrometry

from astropy.io import ascii
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, Column
from astropy.time import Time
from pyraf import iraf
import os, sys
import numpy as np


#astrompar_direc='/data0/astrom/'
#infile=sys.argv[1]
#files=np.genfromtxt(infile,usecols=(0),dtype=str)
#objlist=list(files)

inimage='test.fit'
inhead=fits.getheader(inimage)
'''
wcsreset fr_fdobj.M82.20140228.0032.fits physical
wcsreset fr_fdobj.M82.20140228.0032.fits world

hedit fr_fdobj.M82.20140228.0032.fits WAT0_001 'system=image' ver-
hedit fr_fdobj.M82.20140228.0032.fits WAT1_001 'wtype=tan axtype=ra' ver-
hedit fr_fdobj.M82.20140228.0032.fits WAT2_001 'wtype=tan axtype=dec' ver-

hedit fr_fdobj.M82.20140228.0032.fits RADECSYS 'FK5'   add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits EQUINOX 2000. add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CTYPE1 'RA---TAN' add+ ver- 
hedit fr_fdobj.M82.20140228.0032.fits CTYPE2 'DEC--TAN' add+ ver-

hedit fr_fdobj.M82.20140228.0032.fits CRVAL1 148.89732 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CRVAL2 69.64831 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CRPIX1 911 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CRPIX2 1149 add+ ver-

hedit fr_fdobj.M82.20140228.0032.fits CD1_1 2.2E-4 add+ ver- 
hedit fr_fdobj.M82.20140228.0032.fits CD1_2 1.9E-5 add+ ver- 
hedit fr_fdobj.M82.20140228.0032.fits CD2_1 1.9E-5 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CD2_2 -2.2E-4 add+ ver-

CD1_1   =  -1.3363857514000E-4 / Transformation matrix
CD1_2   =  -3.4441449067900E-6 / no comment
CD2_1   =  3.51080909959001E-6 / no comment
CD2_2   =  -1.3355607527000E-4 / no comment
'''


def wcsresetphysical(group):
#	group=(",".join(group))
	iraf.wcsreset.setParam('image',group)
	iraf.wcsreset.setParam('wcs','physical')
	iraf.wcsreset(group)

def wcsresetworld(group):
#	group=(",".join(group))
	iraf.wcsreset.setParam('image',group)
	iraf.wcsreset.setParam('wcs','world')
	iraf.wcsreset(group)

def hedit(image,par,val)




presecom='sex -c prese.sex '+inimage+' -DETECT_MINAREA=5 -DETECT_THRESH=2.0 -ANALYSIS_THRESH=1.5 -PARAMETERS_NAME=prese.param -CATALOG_NAME=prese.cat -CATALOG_TYPE=ASCII'

presecat=ascii.read('prese.cat')


astromse='sex -c astrom.sex '+inimage+' -DETECT_MINAREA=5 -DETECT_THRESH=2.0 -ANALYSIS_THRESH=2.0 -PARAMETERS_NAME=astrom.param -CATALOG_NAME=astrom.cat -CATALOG_TYPE=ASCII'
scampcom='scamp astrom.cat -c astrom.scamp'
swarpcom='swarp '+inimage+' astrom.swarp -IMAGEOUT_NAME='+inimage+'s'




header




