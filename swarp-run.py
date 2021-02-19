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



def hedit(filename, keywords, value) :
	'''
	hedit(image, keyword, value)
	'''
	from pyraf import iraf
	from astropy.io import fits
	iraf.imutil()
	iraf.imutil.hedit.setParam('images',filename)
	iraf.imutil.hedit.setParam('field',keywords)
	iraf.imutil.hedit.setParam('value',value)
	iraf.imutil.hedit.setParam('add','yes')
	iraf.imutil.hedit.setParam('addonly','yes')
	iraf.imutil.hedit.setParam('verify','no')
	iraf.imutil.hedit(mode='h')


# inputlist, centertime
def newlines(linestr):
    newline,saline,grline = [],[],[]
    for i in lines[-1].split(',') :
        if 'Calib'in i :
            print('sa'+i)
            saline.append('sa'+i)
            grline.append('gr'+i)
            newline.append(i)
	print('Calib-',len(newline),'saCalib-',len(saline),'grCalib-',len(grline))
    return newline,saline,grline

def centertimeheader(inlist) :
	outim=os.path.splitext(inlist[0])[0]+'_'+str(len(inlist))+'_com'+os.path.splitext(inlist[0])[1]
	obsdt,obsmjd=[],[]
	for n in range(len(inlist)):
		header=fits.getheader(inlist[n])
		hobsdate=header['DATE-OBS']
		hmjd=header['MJD']
		obsdt.append(hobsdate)
		obsmjd.append(hmjd)
	tt = Time(obsdt, format='isot', scale='utc')
	ttjd=tt.jd
	ttjdmean=np.mean(ttjd)

	print (ttjd)
	print (ttjdmean)
	ttjdmeanutc=Time(ttjdmean,format='jd',scale='utc')

	putdata,putheader=fits.getdata(outim, header=True)
	# os.system('rm '+putim)
	putheader['DATE-OBS']=ttjdmeanutc.isot
	fits.writeto(putim, putdata, putheader, overwrite=True)


param_dict={
'IMAGEOUT_NAME'     : output,
'COMBINE'           : 'Y',
'COMBINE_TYPE'      : 'MEDIAN',
'SUBTRACT_BACK'     : 'N',
'COPY_KEYWORDS'     : 'OBJECT,FILTER',
'WRITE_FILEINFO'    : 'Y',
'BACK_TYPE'         : 'AUTO',
'BACK_DEFAULT'      : 0.0,
'BACK_SIZE'         : 64,
'BACK_FILTERSIZE'   : 3
'CELESTIAL_TYPE'    : 'NATIVE'
'PROJECTION_TYPE'   : 'TAN'
'RESAMPLE'          : 'Y'
'FSCALE_KEYWORD'    : 'FLXSCALE'
'COPY_KEYWORDS'     : 'OBJECT'
'WRITE_FILEINFO'    : 'N'
}

optstr=''
output=os.path.splitext(inlist[0])[0]+'_'+str(len(inlist))+'_com'+os.path.splitext(inlist[0])[1]
for i in param_dict:
    print(' -{} {}'.format(i,param_dict[i]))
    optstr += ' -{} {}'.format(i,param_dict[i])

swarpcom='swarp -c default.swarp ' + ','.join(saline) + optstr
print(swarpcom)
os.system(swarpcom)



with open('comfiles.txt','r') as file_handle: lines = file_handle.read().splitlines()
for ii in lines:
    if len(ii) == 1 :
        print(ii[0], 'single image')
        os.system('mv '+ii[0]+' gr'+ii[0])
    else :
        ref_image=ii.split(',')[0]
        images_to_align=ii.split(',')[:-1]
        print (ii.split(','))
