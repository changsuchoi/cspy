# read list and make images combined for oneday 3 fits files
# mean exposure start time will be update to DATE-OBS keyword in combined image header
# python imagecombine.py 'gre*.fits' median
# python 3 porting completed, 2019.05.28, Changsu Choi

from astropy.io import ascii
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, Column
from astropy.time import Time
from pyraf import iraf
import os, sys
import numpy as np
import matplotlib.pyplot as plt
#import ds9

# os.system('rm *_com.fits')
# str1=sys.argv[1]
#str1='Calib*ter.fits'
combinestr='median'  # sys.argv[2] #  average|median|lmedian|sum|quadrature|nmodel
#if combinestr not in ['average','median','lmedian','sum','quadrature','nmodel']:
#	print ("give me one of these options,",'average','median','lmedian','sum','quadrature','nmodel')


com='gethead '+str1+' DATE-OBS > obj.list.date'
os.system(com)
#os.system('gethead trreCal*gregister.fits DATE-OBS > obj.list.date')
#os.system('gethead reCal*gregister.fits DATE-OBS > obj.list.date')

filename='obj.list.date'
colnames=('files','obstime')
info=ascii.read(filename,Reader=ascii.NoHeader,names=colnames)
##	info=np.loadtxt(filename)
info=Table(info)
info.sort('obstime')
files=info['files']
obstime=info['obstime']


'''
def imcombine(group,output):
	group=(",".join(group))
#	combine=average
	iraf.imcombine.unlearn()
	iraf.imcombine.setParam('input',group)
	iraf.imcombine.setParam('output',output)
	iraf.imcombine.setParam('combine','median')
	iraf.imcombine.project='no'
#	iraf.imcombine.setParam('combine','average')
	iraf.imcombine.setParam('reject','none')
	iraf.imcombine.setParam('scale','mode')
	iraf.imcombine.setParam('zero','mode')
	iraf.imcombine(group,output=output)
'''
def imcombine(group,output):
	#group=(",".join(group))
	iraf.imcombine(group,output=output,combine=combinestr,project="no",reject="none",scale="none",zero="mode")

obsdt=[]
for n in range(len(files)):
	header=fits.getheader(files[n])
	hobsdate=header['DATE-OBS']
	obsdt.append(hobsdate)
	print (hobsdate	+'\n')

#time conversion to mjd
t = Time(obstime, format='isot', scale='utc')
tjd=t.mjd
#-------------------------------------------------------------
#center time calculation and put it to header
def centertimeheader(inim,putim) :
	obsdt=[]
	for n in range(len(inim)):
		header=fits.getheader(inim[n])
		hobsdate=header['DATE-OBS']
		obsdt.append(hobsdate)
	tt = Time(obsdt, format='isot', scale='utc')
	ttjd=tt.jd
	ttjdmean=np.mean(ttjd)

	print (ttjd)
	print (ttjdmean)
	ttjdmeanutc=Time(ttjdmean,format='jd',scale='utc')

	putdata,putheader=fits.getdata(putim, header=True)
	os.system('rm '+putim)
	putheader['DATE-OBS']=ttjdmeanutc.isot
	fits.writeto(putim, putdata, putheader, overwrite=True)
#---------------------------------------------------------------

def exposuresum(ims,expkey='EXPTIME'):
	exps=[fits.getheader(i)[expkey] for i in ims]
	expsum=np.sum(exps)
	return str(int(expsum))

com=[]
for i in range(len(files)) :
	#com.append(files[i])
	if i==0 : com.append(files[0])

	else :
		# if time between two exposures less than 10 min, then append new file to group	successive files will be combined together
		if (tjd[i]-tjd[i-1]) < (20/1440.) :
			com.append(files[i])
		else :

			if len(com) == 1 :
				output=com[0][:-5]+'_'+str(len(com))+'_com.fits'
				print ('these ',str(len(com)),' files will be combined ',com,' output = ',output+'\n')
				os.system('cp '+com[0]+' '+output)
				com=[]
				com.append(files[i])

			else :
				f=open('tmpcom.list','w')
				for n in com : f.write(n+'\n')
				f.close()
				output=com[0][:-5]+'_'+str(len(com))+'_com.fits'
				imcombine('@tmpcom.list',output)
				centertimeheader(com,output)
				print ('these ',str(len(com)),' files will be combined ',com,' output = ',output+'\n')
				com=[]
				com.append(files[i])

f=open('tmpcom.list','w')
for n in com : f.write(n+'\n')
f.close()
output=com[0][:-5]+'_'+str(len(com))+'_com.fits'
imcombine('@tmpcom.list',output)
centertimeheader(com,output)
print ('these ',str(len(com)),' files will be combined ',com,' output = ',output+'\n')
os.system('rm tmpcom.list')
#output=com[0][:-5]+'_'+str(len(com))+'_com.fits'
#imcombine(,output)

print ('all done, check it out \n')
