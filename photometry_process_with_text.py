# photometry process script based on CSPY files
# test on NGC3367/LOAO/R
# remove unnecessary files first

# 0. file rename and remove
# os_walk_file_arrange.py
# remove and rename, leave only Calib-OBS-OBJ-DATE-TIME-FILTER-EXPTIME.fits
# return Calib*.fits only


# 1. calibfitslist.py
caliblist=glob.glob("Calib*.fits")
oklist, nolist=calibfitslist(caliblist)
oklist.sort()
nolist.sort()
print ('oklist', len(oklist),'nolist', len(nolist))
# return oklist

# 1.a FOV check
# fov-check.py field separation > 10 arcmin
# move them to bad/ directory

# 2. scamp-astrometry.net-result.py
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

#cpunum=4
#if __name__ == '__main__' :
#	p=Pool(cpunum)
#	p.map(scamp_astrometry_net,oklist)

salist=glob.glob('saCalib*.fits')
salist.sort()
print(len(oklist),len(salist))

# 3. imagesetgroup.py
# return clines, salines
lines=epoch_group(oklist)
salines=epoch_group(salist)
len(lines),len(salines)
# 4. alipy-python3-epoch.py
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

# 5. imagecombine.py
# image combine for each epoch with Calib*gregister.fits files
# return Calib*-180_com.fits
glist=glob.glob('Calib*gregister.fits')
glist.sort()
glines=epoch_group(glist)
imagecombine_epoch(glines)
os.system('ls Calib*com.fits | wc -l')
print('imagecombine done')

# 6. swarprun.py
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

# 7. psfexrun.py
# psfex run, get .psf files
# note. MCD30INCH not good for psfex run
# return *.psf
fitslist=glob.glob('*Calib*.fits')
fitslist.sort()
starttime=time.time()
cpunum=4
if __name__ == '__main__' :
	p=Pool(cpunum)
	p.map(psfexrun,fitslist)
endtime=time.time()
duration= endtime-starttime

# 8. se_1st
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
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
magtype=magtypes[0]
refcat='../../ps1-Tonry-NGC3367.cat'
psf=True
imlist=glob.glob('*Calib*.fits')
# multiprocess
'''
imlist.sort()
cpunum=2
if __name__ == '__main__' :
	p=Pool(cpunum)
	p.map(se1st,imlist)
'''
#single
stime=time.time()
badlist=[]
for nn,im in enumerate(imlist) :
	print('=' *60,'\n')
	print(nn+1, 'of ',len(imlist),im)
	if fits.getheader(im).get('UL5_OPTA',default=True):
		try : se1st(im)
		except:
			print(im,'has a problem, add it to badlist')
			badlist.append(im)
	else: pass
print(badlist)
etime=time.time()
duration=etime-stime
print('se_1st done', len(imlist), duration)
#starttime=time.time()
#se1st(im)
#endtime=time.time()
#duration= endtime-starttime
# 16sec for 1 image

# return badlist, and zp related values in header, plots(zp, error, ul5, fov)
# in badlist, saCalib*.fits are major files because 'PV1_5' ~ 'PV1_10' keyword
# remove keywords in header and run se1st again
# number of images must be same (Calib = saCalib) only due to image quality.
delhdrkey


# 9.se_2nd
# se-2nd-phot-function.py
# return *.sef (threshold 1.5), *.dat (final catalog)

#for im in imlist:
#	psf=True
#	se2com(im)
#	secat_zp(im)

cpunum=4
if __name__ == '__main__' :
	p=Pool(cpunum)
	p.map(se2nd,imlist)

# 10. imagecut.py
# image cut size 10 arcmin
#return Calib*_10min_cut.fits
size = 10 # arcmin unit, length of square side
ra,dec=161.645641, 13.750859
for im in imlist:
	imcopy(im,position=(ra,dec))

# 11. alipy-ref2in.py
# registration ref.fits to iuput image, single 2.5sec
# return reg_Calib*.fits
for im in imlist:
	alipy_ref2in(im)

# 11.1 alipy-ref2in.py
# registration ref.fits to input cut file
# return reg_Calib*cut.fits
cutlist=glob.glob('Calib*cut.fits')
cutlist.sort()
for im in cutlist:
	alipy_ref2in(im)

# 12. hotpants.py
# subtraction
# return hd*,hc*.fits
# imlist, cutlist

cpunum=5
if __name__ == '__main__' :
    p=Pool(cpunum)
    p.map(hprun,imlist)
cpunum=5
if __name__ == '__main__' :
    p=Pool(cpunum)
    p.map(hprun,cutlist)

# 13. se_sub.py
# se_sub on subtracted image





































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
