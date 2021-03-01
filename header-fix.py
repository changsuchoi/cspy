import os
import astropy.io.fits as fits
from astropy.time import Time


def oswalkfunc():
	f=open('oswalk.list','w')
	#workDIr = os.path.abspath(b'.')
	for root, dirs, files in os.walk('.'): # os.walk(".", topdown = False):
	   # all files with path names
	   for name in files:
	      #print(os.path.join(root, name))
	      f.write(os.path.join(root, name)+'\n')
	f.close()
	with open('oswalk.list','r') as file_handle: lines = file_handle.read().splitlines()
	print(len(lines),'files')
	return lines

lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']
print (len(fitslist),'fits files')
hkey=['DATE-OBS', 'EXPTIME', 'MJD']
flist1,flist2=[],[]
for im in fitslist:
	print(im)
	h=fits.getheader(im)
	if ('DATE-OBS' in list(h.keys())) and (len(h['DATE-OBS'])>19): pass
	else : flist1.append(im)
	if 'EXPTIME' not in list(h.keys()): flist2.append(im)








h=fits.getheader(im)
for k in hkey:
	if k in list(h.keys()): print(k)
	else:print(k, 'missing')

 # lines = [line.strip() for line in file_handle]

#header keyword remove (PVx_x)
#h.remove('PV1_0')

kwds=['PV1_0', 'PV1_1', 'PV1_2', 'PV1_3', 'PV1_4', 'PV1_5',
		'PV1_6', 'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10',
		'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3', 'PV2_4', 'PV2_5',
		'PV2_6', 'PV2_7', 'PV2_8', 'PV2_9', 'PV2_10']
# kwds=['PV1_5', 'PV1_6', 'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10']
def delhdrkey(im,kwds=kwds):
	h=fits.getheader(im)
	for k in kwds:
		if k in list(h.keys()):
			print('removing',k,h[k])
			fits.delval(im,k)

#for im in imlist:
#	print(im)
#	delhdrkey(im)

kwds=['REFCAT',  'LOWMAG' ,'HIGHMAG',
'NUM_AUTO','ZP_AUTO' ,'ZPE_AUTO','UL5_AUTO',
'NUM_AP3' ,'ZP_AP3'  ,'ZPE_AP3' ,'UL5_AP3' ,
'NUM_AP5' ,'ZP_AP5'  ,'ZPE_AP5' ,'UL5_AP5' ,
'NUM_AP7' ,'ZP_AP7'  ,'ZPE_AP7' ,'UL5_AP7' ,
'NUM_F10' ,'ZP_F10'  ,'ZPE_F10' ,'UL5_F10' ,
'NUM_F15' ,'ZP_F15'  ,'ZPE_F15' ,'UL5_F15' ,
'NUM_F20' ,'ZP_F20'  ,'ZPE_F20' ,'UL5_F20' ,
'NUM_F25' ,'ZP_F25'  ,'ZPE_F25' ,'UL5_F25' ,
'NUM_F30' ,'ZP_F30'  ,'ZPE_F30' ,'UL5_F30' ,
'NUM_OPTA','ZP_OPTA' ,'ZPE_OPTA','UL5_OPTA',
'NUM_PSF' ,'ZP_PSF'  ,'ZPE_PSF' ,'UL5_PSF' ]
#for im in imlist:
#	print(im)
#	delhdrkey(im,kwds=sekwd)
'''
import parmap
imlist.sort()
stime=time.time()
cpunum=2
result=parmap.map(delhdrkey, imlist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('header fix done', len(imlist), duration/60)
'''

# header check and update
'''
import astropy.io.fits as fits
import numpy as np
os.system('gethead *Calib*.fits DATE-OBS EXPTIME')
keyword='EXPTIME'
# keyword = 'DATE-OBS'
for im in caliblist:
	hdr=fits.getheader(im)
	#hdr.get(keyword,default=False)
	if keyword not in list(hdr.keys()):
		hdrval=hdr['EXP_TIME']
		print('header update',keyword,hdrval)
		puthdr(im, keyword, hdrval, hdrcomment='')
	else:pass
keyword='DATE-OBS'
fixlist=[]
for im in flist:
	hdr=fits.getheader(im)
	if len(hdr[keyword]) < 19 :
		print(im)
		fixlist.append(im)


# MCD30ICH DATE-OBS fix
dateobs=fits.getheader(im)['DATE-OBS']+'T'+fits.getheader(im)['UT']
puthdr(im,'DATE-OBS',dateobs)
# SOAO DATE-OBS fix
dateobs=fits.getheader(im)['DATE-OBS']+'T'+fits.getheader(im)['TIME-OBS']
puthdr(im,'DATE-OBS',dateobs)
'''
