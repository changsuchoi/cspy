#this code will registrate all images with names of "cCalibrate*.fits"
#be sure to have files those are comicray calibrated
#be sure to same fields
#usage : python wregister_itelescope.py ref.fits

import numpy as np
import os,sys
from astropy.io import fits
from astropy.io import ascii
import astropy.units as u
import astropy.coordinates as coord
import astropy.units as u
import subprocess
from astropy.table import Table, Column
from astropy.time import Time
from pyraf import iraf
import numpy as np
import matplotlib.pyplot as plt

#ref="b_ref.fits"
ref=sys.argv[1]
os.system('rm rah*.fits')
subprocess.call('ls ahcCalibrate*.fits > obj.list',shell=True)
subprocess.call('ls ahcCalibrate*.fits > obj.list.done',shell=True)
subprocess.call('gethead ahcCal*.fits DATE-OBS > obj.list.date',shell=True)
resnum=subprocess.check_output('ls cCalibrate*.fits | wc -l',shell=True)
print "files astrometry done : ", resnum
colnames=('files','obstime')

filename='obj.list.date'
info=ascii.read(filename,Reader=ascii.NoHeader,names=colnames)
info=Table(info)
info.sort('obstime')
files=info['files']
obstime=info['obstime']

objfile='obj.list.done'
objlist=np.genfromtxt(objfile,usecols=(0),dtype=str)
pathname=os.getcwd()




#"""
#wregister afdobj.iPTF14yb.20140227.0148.fits afdobj.iPTF14yb.20140227.0140.fits #rafdobj.iPTF14yb.20140227.0148.fits
#"""

def wregister(infile,ref,output) :
	#inlist=(",".join(inlist))
	#output=(",".join(output))
	iraf.wregister.unlearn()
	iraf.wregister.setParam('input',infile)
	iraf.wregister.setParam('output',output)
	iraf.wregister.setParam('reference',ref)
	#iraf.wregister.setParam('mode','')
	iraf.wregister(infile,ref,output)


#print gg
#inlist=[]
#ref='ref.fits'
#ref=inlist[0]

inlist=list(files)
for mm in range(len(inlist)) :
	output='r'+inlist[mm]
	wrcom="wregister "+inlist[mm]+" "+output+" "+ref
	infile=inlist[mm]		
	print wrcom
	wregister(infile,ref,output)

#os.system('rm ahcCalibrate*.fits')
print 'all done'	
os.system('ds9 rah*.fits &')	
