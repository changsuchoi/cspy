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

'''  1.b header check '''

''' 2. scamp-astrometry.net-result.py '''
# multiprocessing cpu=4 projection='TPV'
# return saCalib*.fits, salist
# with cpu=4, 3 sec for 1 image process
# after scamp, remove PVX_X keywords in header,
# cut image less than the size of ref.fits
for j,im in enumerate(oklist) :
	contnum=scamp_net(im, astref=False, projection='TAN',DETECT_THRESH='3')
	headmerge(im)
	fits.setval('sa'+im, 'FLXSCALE', value=1)
	fits.setval('sa'+im, 'SCAMPCON', value=contnum, hdrcomment='SCAMP cont. num')
	print(j+1, 'of', len(oklist) )
	print('='*60,'\n')

# multi core with parmap
import parmap
stime=time.time()
cpunum=20
result=parmap.map(scamp_astrometry_net, oklist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('scamp_astrometry.net done', len(imlist), duration)
#cpunum=4
#if __name__ == '__main__' :
#	p=Pool(cpunum)
#	p.map(scamp_astrometry_net,oklist)

salist=glob.glob('saCalib*.fits')
salist.sort()
print(len(oklist),len(salist))


''' 2.a. remove PVx_x keyword for saCalib*.fits '''
# header-fix.py delhdrkey()
kwds=['PV1_0', 'PV1_1', 'PV1_2', 'PV1_3', 'PV1_4', 'PV1_5',
		'PV1_6', 'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10',
		'PV2_0', 'PV2_1', 'PV2_2', 'PV2_3', 'PV2_4', 'PV2_5',
		'PV2_6', 'PV2_7', 'PV2_8', 'PV2_9', 'PV2_10']
for im in salist:
	delhdrkey(im,kwds=kwds)

''' 3. imagesetgroup.py '''
# return clines, salines
lines=epoch_group(oklist)
salines=epoch_group(salist)
len(lines),len(salines)

''' 4. alipy-python3-epoch.py '''
# gregister for each epoch
# single image check
for ii in lines:
    iii=ii[:-1].split(',')
    if len(iii)==1:
        print (iii)
for ii in salines:
    iii=ii[:-1].split(',')
    if len(iii)==1:
        print (iii)

starttime=time.time()
alipy_epoch(lines)
endtime=time.time()
duration= endtime-starttime
glist=glob.glob('Calib*gregister.fits')
glist.sort()
print('alipy work done', duration, len(glist))

''' 5. imagecombine.py '''
# image combine for each epoch with Calib*gregister.fits files
# return Calib*-180_com.fits
glist=glob.glob('Calib*gregister.fits')
glist.sort()
glines=epoch_group(glist)
imagecombine_epoch(glines)
os.system('ls Calib*com.fits | wc -l')
print('imagecombine done')

''' 6. swarprun.py '''
# image combine for each opoch which saCalib*fits files
# EXPTIME keyword
# return saCalib*180_com.fits
salist=glob.glob('saCalib*.fits')
salist.sort()
salines=epoch_group(salist)
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
imlist=glob.glob('*Calib*.fits')
imlist.sort()
stime=time.time()
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
stime=time.time()
cpunum=20
result=parmap.map(scamp_astrometry_net, oklist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('scamp_astrometry.net done', len(imlist), duration)

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
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
magtype=magtypes[0]
refcat='../../ps1-Tonry-NGC3367.cat'
psf=True
imlist=glob.glob('*Calib*.fits')

# single core
'''
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
print('se_1st done', len(imlist), duration)
'''
# multiprocess
import parmap
imlist.sort()
stime=time.time()
cpunum=20
result=parmap.map(se1st, imlist, pm_pbar=True, pm_processes=cpunum)
etime=time.time()
duration=etime-stime
print('se_1st done', len(imlist), duration/60)

badlist=[]
for a,b in zip(imlist,result):
	if b==None:
		print(a,b)
		badlist.append(a)
print (badlist)

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
imlist=glob.glob('*Calib*.fits')

#single
for j,im in enumerate(imlist):
	print(j,'of',len(imlist))
	if os.path.isfile(os.path.splitext(im)[0]+'.dat'):pass
	else: se2nd(im)

stime=time.time()
cpunum=2
if __name__ == '__main__' :
	p=Pool(cpunum)
	p.map(se2nd,imlist)
	p.close()
	p.join()
etime=time.time()
duration=etime-stime
print('se_2nd done', len(imlist), duration)


''' 10. imagecut.py '''
# image cut size 10 arcmin
#return Calib*_10min_cut.fits
size = 10 # arcmin unit, length of square side
ra,dec=161.645641, 13.750859
for im in imlist:
	imcopy(im,position=(ra,dec))

''' 11. alipy-ref2in.py '''
# registration ref.fits to iuput image, single 2.5sec
# return reg_Calib*.fits
for im in imlist:
	alipy_ref2in(im)

''' 11.1 alipy-ref2in.py '''
# registration ref.fits to input cut file
# return reg_Calib*cut.fits
cutlist=glob.glob('Calib*cut.fits')
cutlist.sort()
for im in cutlist:
	alipy_ref2in(im)

''' 12. hotpants.py '''
# subtraction
# return hd*,hc*.fits
# imlist, cutlist

cpunum=2
if __name__ == '__main__' :
    p=Pool(cpunum)
    p.map(hprun,imlist)
	p.close()
	p.join()


cpunum=2
if __name__ == '__main__' :
    p=Pool(cpunum)
    p.map(hprun,cutlist)
	p.close()
	p.join()

''' 13. se_sub.py '''



''' 14. final LC photometry data file '''





































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
