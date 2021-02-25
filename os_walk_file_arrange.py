# file arrange, remove, rename

import os
import astropy.io.fits as fits

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
 # lines = [line.strip() for line in file_handle]

def fnamechange(ii):
	#for CCA250
	i=ii.split('/')[-1]
	head=fits.getheader(ii)
	objname=head['OBJECT']
	dateobs=head['DATE-OBS']
	datestr=dateobs[:4]+dateobs[5:7]+dateobs[8:10]+'-'+dateobs[11:13]+dateobs[14:16]+dateobs[17:20]
	filterstr=head['FILTER']
	exptimestr=str(int(head['EXPTIME']))
	newname='Calib-CCA250-'+objname+'-'+datestr+'-'+filterstr+'-'+exptimestr+'.fits'
	print('cp '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+newname)
	os.system('cp '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+newname)

def LSGTfilechange(ii):
	# From Calib-LSGT-NGC3367-20180519-220208-g-BIN1-W-180-003.fits
	# To   Calib-LSGT-NGC3367-20180519-220208-g-180.fits
	i=ii.split('/')[-1]
	frag=i.split('-')
	frag[0]=='Calib'
#	if frag[1]=='T52' : obs='LSGT'
#	else : obs=frag[1]
	finalname='Calib-LSGT'+'-'+frag[2]+'-'+frag[3]+'-'+frag[4]+'-'+frag[5]+'-'+frag[8]+'.fits'
	os.system('mv '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+finalname)

def iTelfilechange(ii):
	# From Calib-T21-ceouobs.changsu-NGC3367-20161130-042831-R-BIN1-E-180-003.fits
	# To   Calib-T21-NGC3367-20161130-042831-R-180.fits
	i=ii.split('/')[-1]
	frag=i.split('-')
	frag[0]=='Calib'
#	if frag[1]=='T52' : obs='LSGT'
#	else : obs=frag[1]
	#finalname='Calib-'+ frag[1] +'-'+frag[2]+'-'+frag[3]+'-'+frag[4]+'-'+frag[5]+'-'+frag[8]+'.fits'
	finalname='Calib-'+ frag[1] +'-'+frag[3]+'-'+frag[4]+'-'+frag[5]+'-'+frag[6]+'-'+frag[9]+'.fits'
	os.system('mv '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+finalname)

def simplerename(ii,a,b):
	'''
	simplerename(filename, from, to)
	'''
	import os
		#i=ii.split('/')[-1]
	os.system('rename '+a+' '+b+' '+ii)

def oswalknamesep(i):
	filename=i.split('/')[-1]
	head='/'.join(i.split('/')[:-1])+'/'
	return filename, head

###########################################################################
lines= oswalkfunc()
lines.sort()

fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']
files=os.listdir('.')
dirs=[i for i in files if os.path.isdir(i)]

lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']
for i in lines :
	if ('Cal' in i and 'psf' in i) or ('merge.cat' in i) or ('Cal' in i and '.xml' in i) or ('Cal' in i and '.png' in i) or ('Cal' in i and '.cat' in i) or ('Cal' in i and 'seg' in i)  or ('hdre' in i ) or ('reCal' in i ) or ('recCal' in i) or ('wr' in i and '.fit' in i) or ('gregister' in i) :
#	if 'com.cat' in i :
		print('remove', i)
		os.remove(i)
## LSGT
lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']
for i in fitslist :
	if 'cCalib' in i :
		print('rename', i)
		os.system('rename cCalib Calib '+i)

lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']
for i in fitslist :
	if 'Calibrated' in i :
		print('rename', i)
		os.system('rename Calibrated Calib '+i)


lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']

for i in fitslist :
	if 'T52-ceouobs.changsu' in i :
			print('rename', i)
			os.system('rename  T52-ceouobs.changsu LSGT '+i)
	if 'T52-ceouobs.joonho' in i  :
			print('rename', i)
			os.system('rename  T52-ceouobs.joonho LSGT '+i)

lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']

for i in fitslist :
	if ('LSGT' in i) and ('BIN' in i) :
		print('rename', i)
		LSGTfilechange(i)

## CCA250

lines= oswalkfunc()
lines.sort()
for i in lines:
	if 'CCA250' in i and '.new' in i :
		print('rename & remove', i)
		fnamechange(i)
		os.remove(i)

lines= oswalkfunc()
lines.sort()
for i in lines :
	if 'CCA250' in i:
		os.system('rename NGC3367-18 NGC3367-2018 '+i)
		os.system('rename NGC3367-17 NGC3367-2017 '+i)
		os.system('rename Calibrated Calib '+i)
		os.system('rename 0.0.fits 0.fits '+i)
		os.system('rename 00.fits .fits '+i)
		os.system('rename ..fits .fits '+i)

## CCA250 directory and files

os.chdir('CCA250')
os.system('rename 100-c 100c Calib*.fits')
os.system('mv *-m575-* m575/')
os.system('mv *-m625-* m625/')
os.system('mv *-m675-* m675/')
os.system('mv *-m725-* m725/')
os.system('mv *-m775-* m775/')
os.system('mv *-V-* V/')
os.system('mv *-R-* R/')

os.chdir('c')
os.system('rename 100-c 100c Calib*.fits')
os.system('mv *-100c-* ../100c')
os.chdir('..')
os.rmdir('c')

os.system('rename NGC3367-18 NGC3367-2018 Calib*.fits')
os.system('rename NGC3367-17 NGC3367-2017 Calib*.fits')
os.system('rename 0.0.fits 0.fits Calib*.fits')
os.system('rename 00.fits .fits Calib*.fits')
os.system('rename ..fits .fits Calib*.fits')



## itelescope T21
lines= oswalkfunc()
lines.sort()

for i in lines :
	if 'Calib-T21-ceou' in i:
		print('file name :',i)
		iTelfilechange(i)

## MAO SNUCAM
lines= oswalkfunc()
lines.sort()
for i in lines :
	if 'SNUCAM' in i :
		if ('reaCal' in i) or ('reCal' in i) or ('aCalib' in i) or('Calib-MAIDANAK' in i):
			print('remove',i)
			os.remove(i)

## MCD30INCH
lines= oswalkfunc()
lines.sort()
for i in lines:
	if 'MCD30INCH' in i :
		print(i)
		if not 'Calib-MCD30INCH' in i:
			print( 'rename ',i)



		simplerename(i,'Cal-30inch','Calib-MCD30INCH')

'''
!rename Cal-30inch Calib-MCD30INCH Cal*.fits

!rename Calib-30inch Calib-MCD30INCH Cal*.fits

!rename Calib-MCD30inch Calib-MCD30INCH Cal*.fits
'''

## SOAO
lines= oswalkfunc()
lines.sort()
for i in lines:
	if 'SOAO' in i and 'SOAO_FLI' in i:
		print ('rename',i)
		simplerename(i,'SOAO_FLI','SOAO')
	if 'SOAO' in i and 'SOAO_FLI4k' in i:
		print ('rename',i)
		simplerename(i,'SOAO_FLI4k','SOAO')
	if 'SOAO' in i and 'SOAO4k' in i:
		print ('rename',i)
		simplerename(i,'SOAO4k','SOAO')


## DOAO
lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']
print(len(fitslist))
for i in fitslist:
	if 'gregister' in i: os.remove(i)

lines= oswalkfunc()
lines.sort()
