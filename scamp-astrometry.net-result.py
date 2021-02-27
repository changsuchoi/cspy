# scamp astrometry for astrometry.net result file
# first version 2020.3.18 Changsu Choi
# function scamp_net(i) 2020.3.26
# header remove astrometry.net head and add scamp .head 2020.3.26
# refering to
# https://github.com/mommermi/2017Spring_Astroinformatics
# To do :




import os
import astropy.io.fits as fits
import glob
import astropy.io.ascii as ascii
import subprocess
from multiprocessing import Process,Pool


astcodedirec   = '/data7/cschoi/code/cspy/astrom/astromtest/'
astseconfig    = astcodedirec+'astrom.sex'
astseparam     = astcodedirec+'astrom.par'
astscampconfig = astcodedirec+'astrom.scamp'
astswarpconfig = astcodedirec+'astrom.swarp'

# input file name
#i = 'Calibrated-LOAO-NGC3367-20180707-034519-R-60.fits'
#i = 'Calibrated-MCD30inch-NGC3367-20181220-115218-R-60.fits'

# os.system('ls *.fits')
#imlist		= glob.glob(input('image to process\t: '))
#imlist		= glob.glob('Calib*.fits')
#imlist.sort()
#for img in imlist: print(img)

def scamp_net(i,astref=False,projection='TAN',DETECT_THRESH='3'):
	newname='sa'+i
	print('='*60, '\n')
	#os.system('cp '+i+' '+newname)
	iname = i.split('.')[0]
	# source extractor
	secom = 'sex '+i+' -c '+astseconfig+' -PARAMETERS_NAME '+astseparam+\
			' -CATALOG_NAME '+iname+'.ldac '+' -DETECT_THRESH '+DETECT_THRESH
	seout = subprocess.getoutput(secom)
	line1 = [s for s in seout.split('\n') if 'RMS' in s]
	line2 = [s for s in seout.split('\n') if 'Objects: detected' in s]
	skymed, skysig = float(line1[0].split('Background:')[1].split('RMS:')[0]), float(line1[0].split('RMS:')[1].split('/')[0])
	nobj, nse      = float(line2[0].split('Objects: detected ')[1].split('/')[0]), float(line2[0].split('Objects: detected ')[1].split('sextracted')[1])
	print('sextractor working ...')
	print('detected',nobj,'sextracted',nse)
	# scamp
	print('scamp working ...')
	opt1= ' -ASTREF_CATALOG GAIA-EDR3 -SAVE_REFCATALOG Y'
	opt1a=' -ASTREFCAT_NAME astrefcat.cat'
	opt2= ' -CROSSID_RADIUS 2.0'
	opt3= ' -PIXSCALE_MAXERR 1.2'             # Max scale-factor uncertainty
	opt4= ' -POSANGLE_MAXERR 5.0'            # Max position-angle uncertainty (deg)
	opt5= ' -POSITION_MAXERR 1.0'           # Max positional uncertainty (arcmin)
	opt6= ' -SN_THRESHOLDS 10.0,100.0'      # S/N thresholds (in sigmas) for all and high-SN sample
	opt7= ' -FWHM_THRESHOLDS 0.0,100.0'       # FWHM thresholds (in pixels) for sources
	opt8= ' -ELLIPTICITY_MAX 0.5'             # Max. source ellipticity
	opt9= ' -PROJECTION_TYPE '+projection+' '			 # SAME, TPV or TAN
	opt10=' -ASTREF_CATALOG GAIA-EDR3 '
	scampcom='scamp -c '+astscampconfig+' '+iname+'.ldac'+' -ASTREF_CATALOG 2MASS'
	scampcom='scamp -c '+astscampconfig+' '+iname+'.ldac'+' -ASTREF_CATALOG GAIA-DR2 -SAVE_REFCATALOG Y'
	scampcoma='scamp -c '+astscampconfig+' '+iname+'.ldac'+' -ASTREF_CATALOG GAIA-EDR3'+ opt9
	scampcomb='scamp -c '+astscampconfig+' '+iname+'.ldac'+opt1+opt9
	scampcom1='scamp -c '+astscampconfig+' '+iname+'.ldac'+' -ASTREF_CATALOG GAIA-EDR3 -SAVE_REFCATALOG Y '+opt9
	if astref==False: scampcom=scampcoma
	else :
		opt1= ' -ASTREF_CATALOG FILE -ASTREFCAT_NAME '+astref
		scampcomb='scamp -c '+astscampconfig+' '+iname+'.ldac'+opt1+opt9
		scampcom=scampcomb
	print(scampcom)
	scampout=subprocess.getoutput(scampcom)
	line1=[s for s in scampout.split('\n') if 'cont.' in s]
	contnum = scampout.split(line1[0])[1].split('\n')[1].split(' ')[11]
	contnum = scampout.split(line1[0])[1].split('\n')[1].split('"')[1].split(' ')[3]
	print('cont.',contnum)
	print('='*60, '\n')
	os.system('rm '+iname+'.ldac')
	return(contnum)

def scamp_net0(i,projection='TAN',DETECT_THRESH='3'):
	#save refcat used for scamp
	newname='sa'+i
	print('='*60, '\n')
	#os.system('cp '+i+' '+newname)
	iname = i.split('.')[0]
	# source extractor
	secom = 'sex '+i+' -c '+astseconfig+' -PARAMETERS_NAME '+astseparam+\
			' -CATALOG_NAME '+iname+'.ldac '+' -DETECT_THRESH '+DETECT_THRESH
	seout = subprocess.getoutput(secom)
	line1 = [s for s in seout.split('\n') if 'RMS' in s]
	line2 = [s for s in seout.split('\n') if 'Objects: detected' in s]
	skymed, skysig = float(line1[0].split('Background:')[1].split('RMS:')[0]), float(line1[0].split('RMS:')[1].split('/')[0])
	nobj, nse      = float(line2[0].split('Objects: detected ')[1].split('/')[0]), float(line2[0].split('Objects: detected ')[1].split('sextracted')[1])
	print('sextractor working ...')
	print('detected',nobj,'sextracted',nse)
	# scamp
	print('scamp working ...')
	opt1= ' -ASTREF_CATALOG GAIA-EDR3 -SAVE_REFCATALOG Y'
	opt1a=' -ASTREFCAT_NAME astrefcat.cat'
	opt2= ' -CROSSID_RADIUS 2.0'
	opt3= ' -PIXSCALE_MAXERR 1.2'             # Max scale-factor uncertainty
	opt4= ' -POSANGLE_MAXERR 5.0'            # Max position-angle uncertainty (deg)
	opt5= ' -POSITION_MAXERR 1.0'           # Max positional uncertainty (arcmin)
	opt6= ' -SN_THRESHOLDS 10.0,100.0'      # S/N thresholds (in sigmas) for all and high-SN sample
	opt7= ' -FWHM_THRESHOLDS 0.0,100.0'       # FWHM thresholds (in pixels) for sources
	opt8= ' -ELLIPTICITY_MAX 0.5'             # Max. source ellipticity
	opt9= ' -PROJECTION_TYPE '+projection+' '			 # SAME, TPV or TAN
	opt10=' -ASTREF_CATLOG GAIA-EDR3 '
	scampcom1='scamp -c '+astscampconfig+' '+iname+'.ldac'+\
				' -CHECKPLOT_TYPE NONE -ASTREF_CATALOG GAIA-EDR3 -SAVE_REFCATALOG Y '+opt9
	print(scampcom1)
	scampout=subprocess.getoutput(scampcom1)
	line1=[s for s in scampout.split('\n') if 'cont.' in s]
	contnum = scampout.split(line1[0])[1].split('\n')[1].split(' ')[11]
	contnum = scampout.split(line1[0])[1].split('\n')[1].split('"')[1].split(' ')[3]
	print('cont.',contnum,'if greater than 2.5, it was successful')
	print('Save the astrometry reference cat file, find GAIA*.cat')
	print('='*60, '\n')
	os.system('rm '+iname+'.ldac')
	return(contnum)

def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'

def headmerge(i):
	newname='sa'+i
	os.system('cp '+i+' '+newname)
	iname = i.split('.')[0]
	hdr=fits.getheader(i)
	# merging header
	print(i, newname)
	newhead = open(iname+'.head', 'r').readlines()
	#line1 = [s for s in newhead.split('\n') if 'Sorbonne' in s]
	for m in range(len(newhead)):
		if 'Sorbonne' in newhead[m]:
			print(m)
	newhead.remove(newhead[1])
	ff=open(iname+'.head','w')
	ff.writelines(newhead)
	ff.close()
	for m in range(len(hdr)):
		if hdr[m] == '--Start of Astrometry.net WCS solution--':
			n=m
			print(n)
	hdr1=hdr
	del hdr1[n:]
	hdr1.extend(fits.Header.fromtextfile(iname+'.head'), update=True, update_first=True)
	# hdr1.fromtextfile(iname+'.head')
	fits.writeto('sa'+i,fits.getdata(i),hdr1,overwrite=True)
	print(i,' is done.')

'''
	for m in range(len(hdr)):
		if hdr[m] == '--Start of Astrometry.net WCS solution--':
			n=m
			print(n)
	hdr1=hdr
	del hdr1[n:]

#hdr1.fromTxtFile('astromtest.head')
hdr1.extend(fits.Header.fromtextfile(iname+'.head'), update=True, update_first=True)
#hdr1.fromtextfile('astromtest.head',update=True,update_first=True)
fits.writeto('a'+inim,fits.getdata(inim),hdr1)
'''


#f=open('scamp_net_result.txt','w')
# scamp_net(iii)
# headmerge(iii)
def scamp_astrometry_net(im, astref=False, projection='TAN',DETECT_THRESH='3'):
	contnum=scamp_net(im, astref=astref,projection=projection,DETECT_THRESH='3')
	headmerge(im)
	fits.setval('sa'+im, 'FLXSCALE', value=1)
	fits.setval('sa'+im, 'SCAMPCON', value=contnum, hdrcomment='SCAMP cont. num')
	return 'done'
	
def scamp_astrometry_net0(im,projection='TAN',DETECT_THRESH='3'):
	contnum=scamp_net0(im,projection=projection,DETECT_THRESH='3')
	headmerge(im)
	fits.setval('sa'+im, 'FLXSCALE', value=1)
	fits.setval('sa'+im, 'SCAMPCON', value=contnum, hdrcomment='SCAMP cont. num')

#using astrometry reference file
# save ASTREF_CATLOG GAIA-EDR3*.cat
'''
scamp_astrometry_net0('ref.fits',projection='TAN')
gaia=glob.glob('GAIA*.cat')[0]
for i,im in enumerate(oklist[1:]):
	print('='*60)
	print (i+1, 'of' ,len (oklist[1:]))
	scamp_astrometry_net(im,astref=gaia,projection='TAN')
'''

# for non astrometry.net result, like scamp result already done before
def scamp_non_astrometry_net(im, projection='TAN'):
	contnum=scamp_net(im, astref=False, projection=projection)
	swarpcom='swarp -c '+astswarpconfig+' '+im+' -IMAGEOUT_NAME sa'+im
	os.system(swarpcom)
	fits.setval('sa'+im, 'FLXSCALE', value=1)
	fits.setval('sa'+im, 'SCAMPCON', value=contnum, hdrcomment='SCAMP cont. num')


# !swarp t.fits -c default.swarp -PROJECTION_TYPE TPV -IMAGEOUT_NAME=sat.fits

'''
for j,im in enumerate(oklist) :
	contnum=scamp_net(im,projection='TAN')
	headmerge(im)
	fits.setval('sa'+im, 'FLXSCALE', value=1)
	fits.setval('sa'+im, 'SCAMPCON', value=contnum, hdrcomment='SCAMP cont. num')
	print(j+1, 'of', len(oklist) )
	print('='*60,'\n')

salist=glob.glob('saCalib*.fits')
salist.sort()
print ('from oklist', len(oklist), 'salist',len(salist))
'''

'''
if __name__ == '__main__':
	cpunum=4
	if __name__ == '__main__' :
		p=Pool(cpunum)
		p.map(scamp_astrometry_net,oklist)
'''

#f.close()
#hdr1.fromTxtFile('astromtest.head')
#hdr1.extend(fits.Header.fromtextfile(iname+'.head'), update=True, update_first=True)
#hdr1.fromtextfile('astromtest.head',update=True,update_first=True)
#fits.writeto('a'+inim,fits.getdata(inim),hdr1)




#header merge
'''

hdu = fits.open(newname, mode='update', verify='silentfix',
				ignore_missing_end=True)

for m in range(len(hdu[0].header)):
	if hdu[0].header[m] == '--Start of Astrometry.net WCS solution--':
		n=m
		print(n)

del hdu[0].header[n:]
for line in newhead:
	key = line[:8].strip()
	try:
		value = float(line[10:30].replace('\'', ' ').strip())
	except ValueError:
		value = line[10:30].replace('\'', ' ').strip()
	comment = line[30:].strip()
	# print(key,'|',value,'|',comment)
	if key.find('END') > -1:
		break
	hdu[0].header[key] = (str(value), comment)
		#hdu[0].header[key] = (value, comment)

hdu.flush()
'''

#for i in fitslist :
#	fname,head=oswalknamesep(i)






'''
# Default configuration file for SCAMP 2.7.8
# EB 2020-03-17
#

#---------------------------- Reference catalogs ------------------------------

REF_SERVER         vizier.unistra.fr   # Internet addresses of catalog servers
                                       # Possible mirrors include:
                                       # vizier.nao.ac.jp,
                                       # vizier.hia.nrc.ca,
                                       # vizier.ast.cam.ac.uk,
                                       # vizier.iucaa.in,
                                       # vizier.china-vo.org,
                                       # vizier.cfa.harvard.edu and
                                       # viziersaao.chpc.ac.za
ASTREF_CATALOG         2MASS           # NONE,FILE,USNO-A2,USNO-B1,GSC-2.3,
                                       # TYCHO-2,UCAC-4,URAT-1,NOMAD-1,PPMX,
                                       # CMC-15,2MASS,DENIS-3,SDSS-R9,SDSS-R12,
                                       # IGSL,GAIA-DR1,GAIA-DR2,PANSTARRS-1,
                                       # or ALLWISE
ASTREF_BAND            DEFAULT         # Photom. band for astr.ref.magnitudes
                                       # or DEFAULT, BLUEST, or REDDEST
ASTREFMAG_LIMITS       -99.0,99.0      # Select magnitude range in ASTREF_BAND
SAVE_REFCATALOG        N               # Save ref catalogs in FITS-LDAC format?
REFOUT_CATPATH         .               # Save path for reference catalogs

#--------------------------- Merged output catalogs ---------------------------

MERGEDOUTCAT_TYPE      NONE            # NONE, ASCII_HEAD, ASCII, FITS_LDAC
MERGEDOUTCAT_NAME      merged.cat      # Merged output catalog filename

#--------------------------- Full output catalogs ---------------------------

FULLOUTCAT_TYPE        NONE            # NONE, ASCII_HEAD, ASCII, FITS_LDAC
FULLOUTCAT_NAME        full.cat        # Full output catalog filename

#----------------------------- Pattern matching -------------------------------

MATCH                  Y               # Do pattern-matching (Y/N) ?
MATCH_NMAX             0               # Max.number of detections for MATCHing
                                       # (0=auto)
PIXSCALE_MAXERR        1.2             # Max scale-factor uncertainty
POSANGLE_MAXERR        5.0             # Max position-angle uncertainty (deg)
POSITION_MAXERR        1.0             # Max positional uncertainty (arcmin)
MATCH_RESOL            0               # Matching resolution (arcsec); 0=auto
MATCH_FLIPPED          N               # Allow matching with flipped axes?
MOSAIC_TYPE            UNCHANGED       # UNCHANGED, SAME_CRVAL, SHARE_PROJAXIS,
                                       # FIX_FOCALPLANE or LOOSE

#---------------------------- Cross-identification ----------------------------

CROSSID_RADIUS         2.0             # Cross-id initial radius (arcsec)

#---------------------------- Astrometric solution ----------------------------

SOLVE_ASTROM           Y               # Compute astrometric solution (Y/N) ?
PROJECTION_TYPE        SAME            # SAME, TPV or TAN
ASTRINSTRU_KEY         FILTER,QRUNID   # FITS keyword(s) defining the astrom
STABILITY_TYPE         INSTRUMENT      # EXPOSURE, PRE-DISTORTED or INSTRUMENT
CENTROID_KEYS          XWIN_IMAGE,YWIN_IMAGE # Cat. parameters for centroiding
CENTROIDERR_KEYS       ERRAWIN_IMAGE,ERRBWIN_IMAGE,ERRTHETAWIN_IMAGE
                                       # Cat. params for centroid err ellipse
DISTORT_KEYS           XWIN_IMAGE,YWIN_IMAGE # Cat. parameters or FITS keywords
DISTORT_GROUPS         1,1             # Polynom group for each context key
DISTORT_DEGREES        3               # Polynom degree for each group

#---------------------------- Photometric solution ----------------------------

SOLVE_PHOTOM           Y               # Compute photometric solution (Y/N) ?
MAGZERO_OUT            0.0             # Magnitude zero-point(s) in output
MAGZERO_INTERR         0.01            # Internal mag.zero-point accuracy
MAGZERO_REFERR         0.03            # Photom.field mag.zero-point accuracy
PHOTINSTRU_KEY         FILTER          # FITS keyword(s) defining the photom.
MAGZERO_KEY            PHOT_C          # FITS keyword for the mag zero-point
EXPOTIME_KEY           EXPTIME         # FITS keyword for the exposure time (s)
AIRMASS_KEY            AIRMASS         # FITS keyword for the airmass
EXTINCT_KEY            PHOT_K          # FITS keyword for the extinction coeff
PHOTOMFLAG_KEY         PHOTFLAG        # FITS keyword for the photometry flag
PHOTFLUX_KEY           FLUX_AUTO       # Catalog param. for the flux measurement
PHOTFLUXERR_KEY        FLUXERR_AUTO    # Catalog parameter for the flux error

#----------------------------- Source selection -------------------------------

SN_THRESHOLDS          10.0,100.0      # S/N thresholds (in sigmas) for all and
                                       # high-SN sample
FWHM_THRESHOLDS        0.0,100.0       # FWHM thresholds (in pixels) for sources
ELLIPTICITY_MAX        0.5             # Max. source ellipticity
FLAGS_MASK             0x00f0          # Global rejection mask on SEx FLAGS

#------------------------------- WCS headers ----------------------------------

AHEADER_SUFFIX         .ahead          # Filename extension for additional
                                       # input headers
HEADER_SUFFIX          .head           # Filename extension for output headers

#------------------------------- Check-plots ----------------------------------

CHECKPLOT_DEV          PNG             # NULL, XWIN, TK, PS, PSC, XFIG, PNG,
                                       # JPEG, AQT, PDF or SVG
CHECKPLOT_TYPE         FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR
CHECKPLOT_NAME         fgroups,distort,astr_interror2d,astr_interror1d,astr_referror2d,astr_referror1d,astr_chi2,psphot_error # Check-plot filename(s)

#------------------------------ Miscellaneous ---------------------------------

VERBOSE_TYPE           NORMAL          # QUIET, NORMAL, LOG or FULL
WRITE_XML              Y               # Write XML file (Y/N)?
XML_NAME               scamp.xml       # Filename for XML output
NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SCAMP
                                       # 0 = automatic
'''
