# photometry process script based on CSPY files
# test on NGC3367/LOAO/R
# remove unnecessary files first

''' 0. file rename and remove '''
# os_walk_file_arrange.py
# remove and rename, leave only Calib-OBS-OBJ-DATE-TIME-FILTER-EXPTIME.fits
# return Calib*.fits only


''' 1. calibfitslist.py '''
caliblist=glob.glob("Calib*.fits")
oklist, nolist=calibfitslist(caliblist)
oklist.sort()
nolist.sort()
print ('oklist', len(oklist),'nolist', len(nolist))
# return oklist

'''  1.a FOV check '''
# fov-check.py field separation > 10 arcmin
# move them to bad/ directory
ra,dec=71.4270833, -59.2471806 # NGC1672

ra,dec=178.206042,   44.120722 # NGC3938
for im in caliblist : imcenter_offset(im,ra=ra,dec=dec,seplimit=10)
'''  1.b header check '''
# DATE-OBS, EXPTIME, FILTER
# header-fix.py
!gethead Calib*.fits DATE-OBS EXPTIME

''' 2. imagesetgroup.py '''
# return clines, salines
caliblist.sort()
#oklist.sort()
#caliblist==oklist
lines=epoch_group(caliblist)
len(lines)
#salines=epoch_group(salist)
#len(salines)

''' 3. alipy-python3-epoch.py '''
# gregister for each epoch
# single image check
for ii in lines:
    iii=ii[:-1].split(',')
    if len(iii)==1:
        print (iii)

import time
starttime=time.time()
alipy_epoch(lines)
endtime=time.time()
duration= endtime-starttime
glist=glob.glob('Calib*gregister.fits')
glist.sort()
print('alipy work done', duration, len(glist))

''' 4. imagecombine.py '''
# image combine for each epoch with Calib*gregister.fits files
# return Calib*-180_com.fits
glist=glob.glob('Calib*gregister.fits')
glist.sort()
glines=epoch_group(glist)
imagecombine_epoch(glines)
os.system('ls Calib*com.fits | wc -l')
calcomlist=glob.glob('Calib*com.fits')
print('imagecombine done')
print(len(glines), len(calcomlist))

''' 5. scamp-astrometry.net-result.py '''
# multiprocessing cpu=4 projection='TPV'
# return saCalib*.fits, salist
# with cpu=4, 3 sec for 1 image process
# after scamp, remove PVX_X keywords in header,
# cut image less than the size of ref.fits
import time
import parmap
imlist=glob.glob('Calib*.fits')
imlist.sort()
print(len(imlist))

for j,im in enumerate(imlist) :
	contnum=scamp_net(im, astref=False, projection='TAN',DETECT_THRESH='3')
	headmerge(im)
	fits.setval('sa'+im, 'FLXSCALE', value=1)
	fits.setval('sa'+im, 'SCAMPCON', value=contnum, hdrcomment='SCAMP cont. num')
	print(j+1, 'of', len(imlist) )
	print('='*60,'\n')

# multi core with parmap
stime=time.time()
cpunum=5
result=parmap.map(scamp_astrometry_net, imlist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('scamp_astrometry.net done', len(imlist), duration/60)
#cpunum=4
#if __name__ == '__main__' :
#	p=Pool(cpunum)
#	p.map(scamp_astrometry_net,oklist)

salist=glob.glob('saCalib*.fits')
salist.sort()
print(len(oklist),len(salist))

imlist=glob.glob('saCalib*.fits')
for im in imlist:
	print(im)
	delhdrkey(im)

''' 5.a. remove PVx_x keyword for saCalib*.fits '''
# header-fix.py delhdrkey()
kwds=['PV1_0', 'PV1_1', 'PV1_2', 'PV1_3', 'PV1_4', 'PV1_5',
		'PV1_6', 'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10','PV1_11'
		'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3', 'PV2_4', 'PV2_5',
		'PV2_6', 'PV2_7', 'PV2_8', 'PV2_9', 'PV2_10','PV2_11']
for im in salist:
	delhdrkey(im,kwds=kwds)

''' 6. swarprun.py '''
# image combine for each opoch which saCalib*fits files
# EXPTIME keyword
# return saCalib*180_com.fits
salist=glob.glob('saCalib*.fits')
salist.sort()
salines=epoch_group(salist)
print(len(salines))
for ii in salines:
    iii=ii[:-1].split(',')
    if len(iii)==1:
        print (iii)
starttime=time.time()
swarp_epoch(salines)
endtime=time.time()
duration= endtime-starttime
print('swarp_epoch done', duration)
os.system('ls saCalib*com.fits | wc -l')

''' 7. psfexrun.py '''
# psfex run, get .psf files
# note. MCD30INCH not good for psfex run
# return *.psf

'''
cpunum=2
if __name__ == '__main__' :
	p=Pool(cpunum)
	p.map(psfexrun,fitslist)
	p.close()
	p.join()
endtime=time.time()
duration= endtime-starttime
'''
imlist=glob.glob('saCalib*.fits')
imlist.sort()
print(len(imlist))
stime=time.time()
cpunum=3
result=parmap.map(psfexrun, imlist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('psfexrun done', len(imlist), duration/60)

''' 8. se_1st '''
# optimal_aper.py
# se-1st-phot-function.py
# se-1st-command.py
# sextractor run, zp, zperr, plots, se1, FWHM, pscale, ul, skyval, skysig
# return *.se1, *.png, header update
DETECT_MINAREA = str(5)
DETECT_THRESH  = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)
lowmag=13
highmag=19
filname,filerr='R','Rerr'
filname,filerr='r','rerr'
filname,filerr='r_psf','e_r_psf'
filname,filerr='B','Berr'
filname,filerr='V','Verr'
filname,filerr='I','Ierr'
filname,filerr='g','gerr'
filname,filerr='u_psf','e_u_psf'
filname,filerr='g_psf','e_g_psf'
filname,filerr='i','ierr'
filname,filerr='i_psf','e_i_psf'
filname,filerr='z_psf','e_z_psf'

magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
magtype=magtypes[0]
refcat='../../ps1-Tonry-NGC3367.cat'
refcat='../../ps1-Tonry-NGC3938.cat'
psf=True
imlist=glob.glob('saCalib*.fits')

# single core
stime=time.time()
imlist.sort()
badlist=[]
for nn,im in enumerate(imlist) :
	print('=' *60,'\n')
	print(nn+1, 'of ',len(imlist),im)
	if fits.getheader(im).get('UL5_OPTA',default=None)==None:
		try : se1st(im)
		except:
			print(im,'has a problem, add it to badlist')
			badlist.append(im)
	else: pass
print(badlist)
etime=time.time()
duration=etime-stime
print('se_1st done', len(imlist), duration/60)

# multiprocess
import parmap
imlist.sort()
stime=time.time()
cpunum=5
result=parmap.map(se1st, imlist, psf=True, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('se_1st done', len(imlist), duration/60)

badlist=[]
for a,b in zip(imlist,result):
	if b==None:
		print(a,b)
		badlist.append(a)
print (badlist)

for i in badlist:
	print(i)
	os.system('mv *'+os.path.splitext(i)[0][2:]+'* bad/')


# 16sec for 1 image

# return badlist, and zp related values in header, plots(zp, error, ul5, fov)
# in badlist, saCalib*.fits are major files because 'PV1_5' ~ 'PV1_10' keyword
# remove keywords in header and run se1st again
# number of images must be same (Calib = saCalib) only due to image quality.

# after removing PV Keys, run se1st again
# check badlist, then choose to move badlist to bad/

''' 9.se_2nd '''
# se-2nd-phot-function.py
# return *.sef (threshold 1.5), *.dat (final catalog)
os.system('rm *Calib*_seg.fits *Calib*_ap.fits')
imlist=glob.glob('saCalib*.fits')
imlist.sort()
print(len(imlist))
#single
badlist=[]
for j,im in enumerate(imlist):
	print('='*60,'\n')
	print(j+1,'of',len(imlist))
	if os.path.isfile(os.path.splitext(im)[0]+'.dat'):pass
	else:
		try: se2nd(im)
		except:badlist.append(im)
	#se2nd(im)
print (badlist)
print('se2nd finished')

# multi core
import parmap
imlist.sort()
stime=time.time()
cpunum=5
result=parmap.map(se2nd, imlist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('se_2nd done', len(imlist), duration/60)

badlist=[]
for a,b in zip(imlist,result):
	if b=='Pass':
		print(a,b)
		badlist.append(a)
print (badlist)
print('se_2nd done', len(imlist), duration/60)
os.system('rm *ap.fits *seg.fits')

''' 10. imagecut.py '''
# image cut size 10 arcmin
# return Calib*_10min_cut.fits
# imcopy
#for im in imlist:
#	imcopy(im,position=(ra,dec))
os.system('rm *Calib*cut.fits')
size = 10 # arcmin unit, length of square side
ra,dec=161.645641, 13.750859 #NGC3367
ra,dec=71.4270833, -59.2471806 # NGC1672
ra,dec=71.45625, -59.245139 #NGC1672
ra,dec=178.206042,   44.120722 #NGC3938

positions=(ra,dec)
stime=time.time()
cpunum=3
result=parmap.map(imcopy, imlist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
cutcallist=glob.glob('*Calib*cut.fits')
print('imcopy image cut done', len(imlist), len(cutcallist), duration/60)
'''
IrafError: Error running IRAF task imcopy
Error in write: Select error for <Subprocess '/data7/cschoi/util/iraf-2.16.1-2018.11.01/bin.linux64/x_images.e -c', at 7f2948107070>: file descriptors [119]
[Errno 32] Broken pipe
'''
# trim
sz=(20,20)
sz=(10,10)
pos=(ra,dec) # NGC3367
#salist=glob.glob('saCalib*.fits')
#salist.sort()
stime=time.time()
cpunum=2
result=parmap.map(trim, imlist,sizes=sz,positions=pos, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
cutlist=glob.glob('saCalib*cut.fits')
print('trim image cut done', len(imlist),len(cutlist), duration/60)
#trim(im, positions=(False,False), sizes=(10,10))
imlist=cutlist +imlist

#single core
badlist=[]
sz=(10,10)
pos=(ra,dec)
for j,im in enumerate(imlist):
	print('='*60,'\n')
	print(j+1,'of',len(imlist))
	#if os.path.isfile(os.path.splitext(im)[0]+'.dat'):pass
	#else:
	try: trim(im, positions=pos, sizes=sz)
	except:badlist.append(im)
	#se2nd(im)
print (badlist)

''' 11. registration, align '''
scamp astrometry first ('ref.fits')
os.system('mv ref.fits ref.fits.old')
os.system('mv saref.fits ref.fits')


''' 11.1 alipy-ref2in '''
# registration ref.fits to iuput image, single 2.5sec
# return reg_Calib*.fits
#salist=glob.glob('saCalib*.fits')
#salist.sort()
stime=time.time()
cpunum=1
result=parmap.map(alipy_ref2in, imlist,head='reg_', pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
regimlist=glob.glob('reg_*Calib*cut.fits')
print('alipy_ref2in done', len(imlist),len(regimlist), duration/60)

#single core
#salist=glob.glob('saCalib*.fits')
#salist.sort()
result,nolist=[],[]
for n,im in enumerate(imlist):
	print(n+1, 'of', len(imlist))
	r=alipy_ref2in(im)
	if r==None:nolist.append(im)
	result.append(r)
print(len(imlist),len(nolist))

''' 11.2 wcsremap.py '''
# registration ref.fits to input cut file
# return rew_*Calib*cut.fits
# remap2min : remap ref.fits to have the the pixelscale value (least pixel scale -0.001)
#imlist=glob.glob('saCalib*.fits')# + glob.glob('saCalib*cut.fits')
imlist.sort()
print(len(imlist))
remap2min(imlist, refim='ref.fits')
stime=time.time()
cpunum=4
result=parmap.map(wcsremap, imlist, refim='ref.fits',pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
rewlist=glob.glob('rew_*Calib*.fits')
print('wcsremap done', len(imlist),len(rewlist) ,duration/60)

for im in imlist:
	print(im)
	wcsremap(im, refim='saref.fits')
rewlist=glob.glob('rew_*Calib*.fits')
print('wcsremap done', len(imlist),len(rewlist) ,duration/60)

''' 11.3 spalipy '''
# registration ref.fits to input cut file
# return res_*Calib*cut.fits
salist=glob.glob('saCalib*.fits')
salist.sort()
#remap2min(salist, refim='ref.fits')
stime=time.time()
cpunum=5
result=parmap.map(spalipy, salist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
rewlist=glob.glob('res_saCalib*.fits')
print('spalipy done', len(salist),len(rewlist) ,duration/60)

''' 11.4 wcstools remap '''
salist=glob.glob('saCalib*.fits')
salist.sort()
stime=time.time()
cpunum=1
result=parmap.map(remap_ref2in, imlist, refim='ref.fits',pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
remlist=glob.glob('rem_*Calib*.fits')
print('remap_ref2in done', len(imlist),len(remlist), duration/60)

''' 11.5 remap input images '''
imlist=

# list comprehension
#fitslist=glob.glob('*Calib*fits')
#salist=[s for s in fitslist if 'saCalib' in s and 'cut' not in s]
#sacutlist=[s for s in fitslist if 'saCalib' in s and 'cut' in s]

#for im in cutlist:
#	wcsremap(im)

''' 12. hotpants.py '''
# subtraction
# return hd*,hc*.fits
# imlist, cutlist
#salist=glob.glob('saCalib*.fits')
#salist.sort()
head='reg_'
head='rew_'
head='rem_'
# reflist=reflist=['rew_'+s for s in salist]
# reference image = multi registered reference images
stime=time.time()
cpunum=4
result=parmap.map(hprun, imlist,head='rew_',tl=0 ,tu=65000, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
hdlist=glob.glob('hd'+head+'*.fits')
print('hotpants run done', len(imlist),len(hdlist), duration/60)

# reference image = 'ref.fits' only
stime=time.time()
cpunum=2
result=parmap.map(hprun_ref, imlist,refim='rem_ref_05min_cut.fits',tl=-10000 ,tu=200000, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
hdlist=glob.glob('hd*.fits')
print('hotpants run done', len(imlist),len(hdlist), duration/60)

for im in imlist:
	print(im)
	hprun(im,head=head,tl=0, tu=65000,ng=False)
hdlist=glob.glob('hd'+head+'*.fits')
print('hotpants run done', len(imlist),len(hdlist), duration/60)

for im in imlist:
	print(im)
	hprun_ref(im,tl=0, tu=65000,ng=False)
hdlist=glob.glob('hd'+head+'*.fits')
print('hotpants run done', len(imlist),len(hdlist), duration/60)

'''
cpunum=2
if __name__ == '__main__' :
    p=Pool(cpunum)
    p.map(hprun,imlist)
	p.close()
	p.join()
'''
# single





''' 13. se_sub.py '''
# se-target.py
# sub image photometry,
# return mag, magerr or ul5, detection plot

hdlist=glob.glob('hdrew_*Calib*.fits')
hdlist.sort()
print(len(hdlist))

targetname='SN2018kp'
tra,tdec=161.637792, 13.741869
targetname='SN2017gax'
tra,tdec=71.45625, -59.245139
targetname='SN2017ein'
tra,tdec=178.221875, 44.123919
pos=(tra,tdec)
sizes=(5,5)
sz=sizes
thead=''
for j in tblcols: thead+=j+'\t'
for h in hdrkey: thead+=h+'\t'
resulthead= '#FILE'+'\t'+'MJD'+'\t'+'OBSERVATORY'+'\t'+'FILTER'+'\t'+thead+'\t'+'sep'+'\t'+'detect'+'\n'
print(resulthead)
f=open(targetname+'-'+os.getcwd().split('/')[-2]+'-'+os.getcwd().split('/')[-1]+'.dat','w')
f.write(resulthead)
badsesublist=[]
for n,im in enumerate(hdlist):
	print('='*60,'\n')
	try:
		body,detect=sub_target_phot(im, tra, tdec, sep=5.0)
		print( im, detect,'\n\n',body)
		f.write(body)
		trimstamp(im, positions=pos, sizes=sz)
		fitplot_target_stamp(im,tra,tdec,sizes)
	except:
		print('error, adding to badsesublist')
		badsesublist.append(im)
f.close()
print(badsesublist)
print('and finally comes to the end ...')





''' 14. final LC photometry data file '''


def ds9list(imlist):
	ll=''
	for l in imlist: ll+=l +' '
	os.system('ds9 '+ ll+' &')


































import time
starttime=time.time()
psfexrun(im)
endtime=time.time()
duration= endtime-starttime


starttime=time.time();
ls saCalib* | wc -l
endtime=time.time();os.system('ls saCalib* | wc -l')


from multiprocessing import Process,Pool
cpunum=4
if __name__ == '__main__' :
p=Pool(cpunum)
p.map(hotpantsrun,infile)
