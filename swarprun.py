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

os.system('swarp -d > default.swarp')

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

def swarpcom(imlist):
	newdt, newname=centertimeheader(imlist)
	param_dict={
	'IMAGEOUT_NAME'     : newname,
	'COMBINE'           : 'Y',
	'COMBINE_TYPE'      : 'MEDIAN',
	'SUBTRACT_BACK'     : 'N',
	'WRITE_FILEINFO'    : 'Y',
	'BACK_TYPE'         : 'AUTO',
	'BACK_DEFAULT'      : '0.0',
	'BACK_SIZE'         : '64',
	'BACK_FILTERSIZE'   : '3',
	'CELESTIAL_TYPE'    : 'NATIVE',
	'PROJECTION_TYPE'   : 'TPV',
	'RESAMPLE'          : 'Y',
	'FSCALE_KEYWORD'    : 'FLXSCALE',
	'COPY_KEYWORDS'     : 'OBJECT,DATE-OBS,FILTER',
	'WRITE_FILEINFO'    : 'Y',
	'WRITE_XML'         : 'N',
	}
	#
	#output=os.path.splitext(inlist[0])[0]+'_'+str(len(inlist))+'_com'+os.path.splitext(inlist[0])[1]
	optstr=''
	for i in param_dict:
		#print(' -{} {}'.format(i,param_dict[i]))
		optstr += ' -{} {}'.format(i,param_dict[i])
	swarpcom='swarp -c default.swarp ' + ','.join(imlist) + optstr
	#print(swarpcom)
	os.system(swarpcom)



salist=glob.glob('saCalib*.fits')
salist.sort()
salines= epoch_group(salist)
#salines[-1][:-1].split(',')

for i in salines:
	ii=i[:-1].split(',')
	if len(ii)==1 :
		if os.path.isfile(centertimeheader(ii)[1]) : pass
		else:
			print("single image")
			os.system('cp '+ii[0]+' '+centertimeheader(ii)[1])
			#puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])
	else:
		if os.path.isfile(centertimeheader(ii)[1]) : pass
		else:
			print(len(ii),'images will be combined', centertimeheader(ii)[1])
			swarpcom(ii)
			puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])
