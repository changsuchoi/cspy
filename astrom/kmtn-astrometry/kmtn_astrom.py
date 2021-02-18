# python code for KMTNET data astrometry (specially CTIO data)
# read kmtnet_astrom.txt first to understand the order and process
# 2015.09.17 Changsu Choi 

from astropy.io import ascii
import numpy as np
import os,sys
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
import astropy.units as u
import subprocess
import pp


#os.system('gethead kmtc.20150218.00503*.fits ra dec filter object exptime date-obs > info.txt')

info=ascii.read('info.txt')
addlist=info['col1']
ra=info['col2']
dec=info['col3']
filters=info['col4']
obj=info['col5']
exptime=info['col6']
dateobs=info['col7']

'''
def mefcr :
	num=addlist[n][14:-5]
	rad=coord.Angle(ra[n],unit=u.hour)	
	radd=rad.degree
	decd=coord.Angle(dec[n],unit=u.deg)
	decdd=decd.degree	
	# python 
	makemef='python kmtn_makemef.py '+num
	resetcrval='python kmtn_resetcrval.py '+ num+'.fits -c '+str(radd)+','+str(decdd)
	os.system(makemef)
	os.system(resetcrval)
'''
for n in range(len(addlist)):
	num=addlist[n][14:-5]
	rad=coord.Angle(ra[n],unit=u.hour)	
	radd=rad.degree
	decd=coord.Angle(dec[n],unit=u.deg)
	decdd=decd.degree	
	# python 
	makemef='python kmtn_makemef.py '+num
	resetcrval='python kmtn_resetcrval.py '+ num+'.fits -c '+str(radd)+','+str(decdd)
	os.system(makemef)
	os.system(resetcrval)


# sextractor
#sexcom= 'sex '+num+'.fits -c  kmtnet.sex -CATALOG_NAME '+num+'.cat -HEADER_SUFFIX NONE -DETECT_THRESH 50.0 -ANALYSIS_THRESH 50.0 -SATUR_LEVEL 60000.0 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE weight.fits'
#scampcom='scamp '+num+'.cat -c kmtnet.scamp -ASTREF_CATALOG 2MASS -POSITION_MAXERR 20.0 -CROSSID_RADIUS 5.0 -DISTORT_DEGREES 3 -PROJECTION_TYPE TPV -AHEADER_GLOBAL kmtnet_global_ctio.ahead -CHECKPLOT_TYPE NONE' 

def sexscamp(files) :
	threshold='30'
	num=files[14:-5]
	sexcom= 'sex '+num+'.fits -c kmtnet.sex -CATALOG_NAME '+num+'.cat -HEADER_SUFFIX NONE -DETECT_THRESH 20.0 -ANALYSIS_THRESH 20.0 -SATUR_LEVEL 60000.0 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE weight.fits'
	scampcom='scamp '+num+'.cat -c kmtnet.scamp -ASTREF_CATALOG UCAC-3 -POSITION_MAXERR 20.0 -CROSSID_RADIUS 5.0 -DISTORT_DEGREES 3 -PROJECTION_TYPE TPV -AHEADER_GLOBAL kmtnet_global_ctio.ahead -CHECKPLOT_TYPE ASTR_REFSYSMAP,FGROUPS,DISTORTION,ASTR_REFERROR2D,ASTR_REFERROR1D -CHECKPLOT_NAME astr_refsysmap,fgroups,distort,astr_referror2d,astr_referror1d -STABILITY_TYPE INSTRUMENT' 
	os.system(sexcom)
	os.system(scampcom)



def set4astrom(files) : 
	from astropy.io import fits
	from astropy.io import ascii

	print files
	num=files[14:-5]
	f=open(num+'.head','r')
	lines=f.readlines()
	f.close()

	f=open(num+'.kk.head','w')
	for line in lines[0:51] : f.write(line)
	f.close()
	data=fits.getdata(num+'.kk.fits')
	hdr=fits.getheader(num+'.kk.fits')
	hdr.fromTxtFile(num+'.kk.head')
	newfile='a'+num+'.kk.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	f=open(num+'.mm.head','w')
	for line in lines[51:102] : f.write(line)
	f.close()
	data=fits.getdata(num+'.mm.fits')
	hdr=fits.getheader(num+'.mm.fits')
	hdr.fromTxtFile(num+'.mm.head')
	newfile='a'+num+'.mm.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	f=open(num+'.nn.head','w')
	for line in lines[102:153] : f.write(line)
	f.close()
	data=fits.getdata(num+'.nn.fits')
	hdr=fits.getheader(num+'.nn.fits')
	hdr.fromTxtFile(num+'.nn.head')
	newfile='a'+num+'.nn.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	f=open(num+'.tt.head','w')
	for line in lines[153:] : f.write(line)
	f.close()
	data=fits.getdata(num+'.tt.fits')
	hdr=fits.getheader(num+'.tt.fits')
	hdr.fromTxtFile(num+'.tt.head')
	newfile='a'+num+'.tt.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	os.system('swarp -c kmtnet.swarp a'+num+'*.fits -IMAGEOUT_NAME a'+files)		
	os.system('ds9 a'+files+' &')

##  sextractor and scamp
for n in range(len(addlist)):
	sexscamp(addlist[n])

## final file making
for files in addlist : set4astrom(files)


'''
##  header edition
for n in range(len(addlist)):
	num=addlist[n][14:-5]
	hdr=fits.getheader(addlist[n])
	data=fits.getdata(num+'.fits')
	hdr.fromTxtFile('006022.head')
	newfile='a'+addlist[n]
	fits.writeto(newfile,data,hdr,clobber=True)
'''
