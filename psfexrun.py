# psfexrun.py
# usage : python psfexrun.py imagelist
# Changsu Choi

import os,sys
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.io.votable import parse
from astropy.time import Time
import matplotlib.pyplot as plt
from pyraf import iraf
from astropy.nddata import Cutout2D
from multiprocessing import Process,Pool


#os.system("rm psf*.fits *.psf *.xml")

# infile= sys.argv[1]
# inlist= np.genfromtxt(infile,usecols=(0),dtype=str)
# inlist=list(inlist)
# os.system("cp /data0/code/psfex.config/* .")

psfexconfigdir="/data7/cschoi/code/cspy/psfex.config/"
#configdir='/home/changsu/code/psfex.config/'
def psfexxml(xmlfile):
	votable=parse(xmlfile)
	table=votable.get_first_table()
	data= table.array
	data['FWHM_Mean']
	fmp=data['FWHM_Mean'][0]
	fmw=data['FWHM_WCS_Mean'][0]
	ps=data['PixelScale_WCS_Mean'][0]
	em=data['Ellipticity_Mean'][0]
#	print "FWHM in pixel, ",fmp # , " FWHM in Arcsec, ",fmw, " pixel scale, ",ps
	return fmp,fmw,ps,em

def pixelscale(i):
	cd11 = fits.getheader(i)['CD1_1']
	cd12 = fits.getheader(i)['CD1_2']
	cd21 = fits.getheader(i)['CD2_1']
	cd22 = fits.getheader(i)['CD2_2']
	pixscale=round(np.sqrt(cd11**2 + cd21**2) *3600 ,4)
	print('Pixel scale =', pixscale,'\"')
	return pixscale

def imcopy(inim,outname):
	region='[100:125,100:125]'
	#outname='tr'+inim
	chinim=inim+region
	iraf.imcopy(chinim,output=outname)

def psfex(i):
	arcsec5 = str(round(5 / pixelscale(i),2))
	fn=os.path.splitext(i)[0]
	presecom = 'sex -c '+psfexconfigdir+'prepsfex.sex '+i+\
			' -CATALOG_NAME '+fn+'.cat -PARAMETERS_NAME '+psfexconfigdir+'prepsfex.param '+\
			'-FILTER_NAME '+psfexconfigdir+'default.conv' + ' -PHOT_APERTURES '+arcsec5
	opt1 = ' -SAMPLE_FWHMRANGE 1.0,30.0 -SAMPLE_VARIABILITY 0.5 -SAMPLE_MINSN 5 -SAMPLE_MAXELLIP 1.0 '
	#opt2 = ' -CHECKPLOT_DEV PNG -CHECKPLOT_TYPE FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS '
	#opt3 = ' -CHECKPLOT_NAME fwhm, ellipticity, counts, countfrac, chi2, resi '
	#opt4 = ' -CHECKIMAGE_TYPE CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS,MOFFAT,-MOFFAT,-SYMMETRICAL '
	opt4 = ' -CHECKIMAGE_TYPE SNAPSHOTS '
	opt4 = ' -CHECKIMAGE_TYPE NONE '
	#opt5 = ' -CHECKIMAGE_NAME chi,proto,samp,resi,snap,moffat,submoffat,subsym '
	opt5 = ' -CHECKIMAGE_NAME snap '
	#+fn+'.psfex_chi.fits '+fn+'.psfex_proto.fits '+fn+'.psfex_sample.fits '+fn+'.psfex_resi.fits '+fn+'.psfex_snap.fits '
	opt6 = ' -XML_NAME '+fn+'.psfex.xml '
	psfexcom = 'psfex -c '+psfexconfigdir+'default.psfex '+fn+'.cat'+opt1+opt6+opt4
	os.system(presecom)
	os.system(psfexcom)
	os.system('rm '+fn+'.cat)
	print (i, ' fwhm value is ', psfexxml(fn+'.psfex.xml'))
	#imcopy('snap_'+i,'psf-'+i)
	#imcopycom='imcopy snap_'+i+'[100:125,100:125] psf-'+i
	#os.system(imcopycom)

def psfexrun(im):
	print('='*60,'\n')
	if os.path.isfile(os.path.splitext(im)[0]+'.psf') : pass
	else: psfex(im)
	os.system('rm *psfex.xml')
	return 'Done'

inlist=glob.glob("*Calib*.fits")
inlist.sort()

for j,inimage in enumerate(inlist) :
	print('='*60,'\n')
	print(j+1,'of ',len(inlist))
	if os.path.isfile(os.path.splitext(inimage)[0]+'.psf') : pass
	else: psfex(inimage)



os.system('rm *psfex.xml')


print ('all done')

#pre source extractor

'''
# Simple configuration file for SExtractor prior to PSFEx use
# only non-default parameters are present.
# EB 2007-08-01
#

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME     prepsfex.cat   # Catalog filename
CATALOG_TYPE     FITS_LDAC      # FITS_LDAC format
PARAMETERS_NAME  prepsfex.param # name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_MINAREA   5              # minimum number of pixels above threshold
DETECT_THRESH    3              # a fairly conservative threshold
ANALYSIS_THRESH  3              # idem

FILTER           Y              # apply filter for detection ("Y" or "N")?
FILTER_NAME      default.conv   # name of the file containing the filter

#-------------------------------- WEIGHTing ----------------------------------
#-------------------------------- FLAGging -----------------------------------
#------------------------------ Photometry -----------------------------------

PHOT_APERTURES   10             # <- put the referrence aperture diameter here
SATUR_LEVEL      60000.0        # <- put the right saturation threshold here
GAIN             1.0              # <- put the detector gain in e-/ADU here

#------------------------- Star/Galaxy Separation ----------------------------
#------------------------------ Background -----------------------------------
#------------------------------ Check Image ----------------------------------
#--------------------- Memory (change with caution!) -------------------------
#------------------------------- ASSOCiation ---------------------------------
#----------------------------- Miscellaneous ---------------------------------
'''

# psfex config
'''
# Default configuration file for PSFEx 3.9.0
# EB 2010-10-10
#

#-------------------------------- PSF model ----------------------------------

BASIS_TYPE      PIXEL_AUTO      # NONE, PIXEL, GAUSS-LAGUERRE or FILE
BASIS_NUMBER    20              # Basis number or parameter
PSF_SAMPLING    0.0             # Sampling step in pixel units (0.0 = auto)
PSF_ACCURACY    0.01            # Accuracy to expect from PSF "pixel" values
PSF_SIZE        25,25           # Image size of the PSF model
CENTER_KEYS     X_IMAGE,Y_IMAGE # Catalogue parameters for source pre-centering
PHOTFLUX_KEY    FLUX_APER(1)    # Catalogue parameter for photometric norm.
PHOTFLUXERR_KEY FLUXERR_APER(1) # Catalogue parameter for photometric error

#----------------------------- PSF variability -----------------------------

PSFVAR_KEYS     X_IMAGE,Y_IMAGE # Catalogue or FITS (preceded by :) params
PSFVAR_GROUPS   1,1             # Group tag for each context key
PSFVAR_DEGREES  2               # Polynom degree for each group

#----------------------------- Sample selection ------------------------------

SAMPLE_AUTOSELECT  Y            # Automatically select the FWHM (Y/N) ?
SAMPLEVAR_TYPE     SEEING       # File-to-file PSF variability: NONE or SEEING
SAMPLE_FWHMRANGE   1.0,30.0     # Allowed FWHM range
SAMPLE_VARIABILITY 0.5          # Allowed FWHM variability (1.0 = 100%)
SAMPLE_MINSN       3           # Minimum S/N for a source to be used
SAMPLE_MAXELLIP    0.3          # Maximum (A-B)/(A+B) for a source to be used

#------------------------------- Check-plots ----------------------------------

CHECKPLOT_DEV       PNG         # NULL, XWIN, TK, PS, PSC, XFIG, PNG,
                                # JPEG, AQT, PDF or SVG
CHECKPLOT_TYPE      FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS        # FWHM,ELLIPTICITY,COUNTS, COUNT_FRACTION, CHI2, RESIDUALS # or NONE
CHECKPLOT_NAME      NONE		#
#CHECKPLOT_NAME      fwhm.png, ellipticity.png, counts.png, countfrac.png, chi2.png, resi.png

#------------------------------ Check-Images ---------------------------------

CHECKIMAGE_TYPE CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS,MOFFAT,-MOFFAT,-SYMMETRICAL       # CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS,MOFFAT,-MOFFAT,-SYMMETRICAL
                                # Check-image types
CHECKIMAGE_NAME snap.fits       # chi.fits,proto.fits,samp.fits,resi.fits,snap.fits,moffat.fits,submoffat.fits,subsym.fits
                                # Check-image filenames

#----------------------------- Miscellaneous ---------------------------------

PSF_DIR                         # Where to write PSFs (empty=same as input)
PSF_SUFFIX      .psf            # Filename extension for output PSF filename
VERBOSE_TYPE    NORMAL          # can be QUIET,NORMAL,LOG or FULL
WRITE_XML       Y               # Write XML file (Y/N)?
XML_NAME        psfex.xml       # Filename for XML output
NTHREADS        0               # Number of simultaneous threads for
                                # the SMP version of PSFEx
                                # 0 = automatic
'''
