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

# os.system('rm *_com.fits')
# str1=sys.argv[1]
combinestr='median'  # sys.argv[2] #  average|median|lmedian|sum|quadrature|nmodel
#if combinestr not in ['average','median','lmedian','sum','quadrature','nmodel']:
#	print ("give me one of these options,",'average','median','lmedian','sum','quadrature','nmodel')



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
	group=(",".join(group))
	iraf.imcombine(group,output=output,combine=combinestr,
					project="no",reject="none",scale="none",zero="mode")

#-------------------------------------------------------------
#center time calculation and put it to header
def centertimeheader(inim):
	obsdt=[]
	for n in range(len(inim)):
		header=fits.getheader(inim[n])
		hobsdate=header['DATE-OBS']
		obsdt.append(hobsdate)
	tt = Time(obsdt, format='isot', scale='utc')
	ttjd=tt.jd
	ttjdmean=np.mean(ttjd)
	#print (ttjd)
	#print (ttjdmean)
	ttjdmeanutc=Time(ttjdmean,format='jd',scale='utc')
	datetimestr=ttjdmeanutc.isot[:4]+ttjdmeanutc.isot[5:7]+\
				ttjdmeanutc.isot[8:10]+'-'+ttjdmeanutc.isot[11:13]+\
				ttjdmeanutc.isot[14:16]+ttjdmeanutc.isot[17:19]
	nameset=os.path.splitext(inim[0])[0].split('-')
	output=nameset[0]+'-'+nameset[1]+'-'+nameset[2]+'-'+\
			datetimestr+'-'+nameset[5]+'-'+exposuresum(inim)+\
			'_com.fits'
	#putdata,putheader=fits.getdata(putim, header=True)
	#os.system('rm '+putim)
	#putheader['DATE-OBS']=ttjdmeanutc.isot
	#fits.writeto(putim, putdata, putheader, overwrite=True)
	return ttjdmeanutc.isot, output
#---------------------------------------------------------------
def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'

def exposuresum(ims,expkey='EXPTIME'):
	exps=[fits.getheader(i)[expkey] for i in ims]
	expsum=np.sum(exps)
	return str(int(expsum))

glist=glob.glob('Calib*gregister.fits')
glist.sort()
glines=epoch_group(glist)

# EXP_TIME EXPTIME header keyword
# alllist=glob.glob('*Calib*.fits')
# for m in alllist:
# 		puthdr(m,'EXPTIME', fits.getheader(m)['EXP_TIME'])

for i in glines:
	ii=i[:-1].split(',')
	if len(ii)==1 :
		if os.path.isfile(centertimeheader(ii)[1]) : pass
		else:
			print("single image")
			os.system('cp '+ii[0]+' '+centertimeheader(ii)[1])
			puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])
	else:
		if os.path.isfile(centertimeheader(ii)[1]) : pass
		else:
			print(len(ii),'images will be combined', centertimeheader(ii)[1])
			imcombine(i[:-1], centertimeheader(ii)[1])
			puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])
os.system('rm *gregister.fits')
