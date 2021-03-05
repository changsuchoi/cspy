# 1st photometry
# estimate FWHM value and auto mag zp
import os
import glob
import astropy.io.fits as fits
import numpy as np
import subprocess
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clipping
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Process,Pool


# input files, config and params
seconfigdir ='/data7/cschoi/code/cspy/sex.config/'
seconfig    ='se1.sex'
separam     ='se1.param'
separam_noPSF = 'se1_noPSF.param'
seconv      ='default.conv'
sennw       ='default.nnw'
DETECT_MINAREA = str(3)
DETECT_THRESH  = str(1.5)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)

psf=True
# source extractor command
def se2com(im,psf=psf):
	#os.system('rm *Calib*_seg.fits *Calib*_ap.fits')
	PSCALE=fits.getheader(im)['PSCALE']
	fwhm_pix=fits.getheader(im)['FWHM_PIX']
	opt_ap=fits.getheader(im)['OPT_AP']
	aper_list=[3,5,7]
	aper_list_fwhm=[1.0,1.5,2.0,2.5,3.0]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	#aper_input = aper_input[:-1]
	for i in aper_list_fwhm: aper_input += '{},'.format(round(i*fwhm_pix,1))
	aper_input = aper_input+str(opt_ap)
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.sef'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	opt5a=' -CHECKIMAGE_TYPE NONE '
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	opt9=' -SEEING_FWHM '+str(round(fwhm_pix*PSCALE))+ ' '
	if psf==True:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 + im
	else:
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 + im
	os.system(secommand)

def secat_zp(im):
	from astropy.table import Table
	from astropy import units as u
	from astropy.table import Column
	fn = os.path.splitext(im)[0]
	hdr=fits.getheader(im)
	sefcat=ascii.read(fn+'.sef')
	print('Making final catalog',fn+'.dat')
	sefcat['MAG_AUTO']=sefcat['MAG_AUTO']+hdr['ZP_AUTO']
	sefcat['MAGERR_AUTO']=np.sqrt(sefcat['MAGERR_AUTO']**2 + hdr['ZPE_AUTO']**2)
	sefcat['MAG_APER']=sefcat['MAG_APER']+hdr['ZP_AP3']
	sefcat['MAGERR_APER']=np.sqrt(sefcat['MAGERR_APER']**2 + hdr['ZPE_AP3']**2)
	sefcat['MAG_APER_1']=sefcat['MAG_APER_1']+hdr['ZP_AP5']
	sefcat['MAGERR_APER_1']=np.sqrt(sefcat['MAGERR_APER_1']**2 + hdr['ZPE_AP5']**2)
	sefcat['MAG_APER_2']=sefcat['MAG_APER_2']+hdr['ZP_AP7']
	sefcat['MAGERR_APER_2']=np.sqrt(sefcat['MAGERR_APER_2']**2 + hdr['ZPE_AP7']**2)
	sefcat['MAG_APER_3']=sefcat['MAG_APER_3']+hdr['ZP_F10']
	sefcat['MAGERR_APER_3']=np.sqrt(sefcat['MAGERR_APER_3']**2 + hdr['ZPE_F10']**2)
	sefcat['MAG_APER_4']=sefcat['MAG_APER_4']+hdr['ZP_F15']
	sefcat['MAGERR_APER_4']=np.sqrt(sefcat['MAGERR_APER_4']**2 + hdr['ZPE_F15']**2)
	sefcat['MAG_APER_5']=sefcat['MAG_APER_5']+hdr['ZP_F20']
	sefcat['MAGERR_APER_5']=np.sqrt(sefcat['MAGERR_APER_5']**2 + hdr['ZPE_F20']**2)
	sefcat['MAG_APER_6']=sefcat['MAG_APER_6']+hdr['ZP_F25']
	sefcat['MAGERR_APER_6']=np.sqrt(sefcat['MAGERR_APER_6']**2 + hdr['ZPE_F25']**2)
	sefcat['MAG_APER_7']=sefcat['MAG_APER_7']+hdr['ZP_F30']
	sefcat['MAGERR_APER_7']=np.sqrt(sefcat['MAGERR_APER_7']**2 + hdr['ZPE_F30']**2)
	sefcat['MAG_APER_8']=sefcat['MAG_APER_8']+hdr['ZP_OPTA']
	sefcat['MAGERR_APER_8']=np.sqrt(sefcat['MAGERR_APER_8']**2 + hdr['ZPE_OPTA']**2)
	if 'MAG_PSF' in sefcat.colnames:
		sefcat['MAG_PSF']=sefcat['MAG_PSF']+hdr['ZP_PSF']
		sefcat['MAGERR_PSF']=np.sqrt(sefcat['MAGERR_PSF']**2 + hdr['ZPE_PSF']**2)
	sefcat.write(fn+'.dat',format='ascii.commented_header',overwrite=True)

def se2nd(im):
	psf=True
	se2com(im)
	secat_zp(im)
	return 'Done'


'''
	param_dict={
	'CATALOG_NAME'     : fn+'.se1',
	'PARAMETERS_NAME'  : seconfigdir+'se1.param',
	'DETECT_MINAREA'   : DETECT_MINAREA,
	'DETECT_THRESH'    : DETECT_THRESH,
	'FILTER_NAME'      : seconfigdir+seconv,
	'DEBLEND_NTHRESH'  : DEBLEND_NTHRESH,
	'DEBLEND_MINCONT'  : DEBLEND_MINCONT,
	'BACK_TYPE'        : 'AUTO',
	'BACK_DEFAULT'     : '0.0',
	'BACK_SIZE'        : '64',
	'BACK_FILTERSIZE'  : '3',
	'BACKPHOTO_TYPE'   : 'LOCAL',
	'PHOT_APERTURES'   : aper_input,
	'SATUR_LEVEL'      : '60000',
	'GAIN'             : '1.0',
	'PIXEL_SCALE'      : str(PSCALE),
	'SEEING_FWHM'      : '1.2',
	'STARNNW_NAME'     :seconfigdir+sennw,
	'CHECKIMAGE_TYPE'  : 'SEGMENTATION,OBJECTS,BACKGROUND',
	'CHECKIMAGE_NAME'  : fn+'_seg.fits,'+fn+'_obj.fits,'+fn+'bg.fits',
	'PSF_NAME'         : fn+'.psf'
	}
	optstr=''
	for i in param_dict:
		#print(' -{} {}'.format(i,param_dict[i]))
		optstr += ' -{} {}'.format(i,param_dict[i])
	secom='sex -c '+seconfigdir+seconfig +' '+ im + optstr
'''
'''
# Default configuration file for SExtractor 2.19.5
# EB 2014-03-19
#

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME     test.cat       # name of the output catalog
CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
PARAMETERS_NAME  default.param  # name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
DETECT_MINAREA   5              # min. # of pixels above threshold
DETECT_MAXAREA   0              # max. # of pixels above threshold (0=unlimited)
THRESH_TYPE      RELATIVE       # threshold type: RELATIVE (in sigmas)
                                # or ABSOLUTE (in ADUs)
DETECT_THRESH    1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

FILTER           Y              # apply filter for detection (Y or N)?
FILTER_NAME      default.conv   # name of the file containing the filter
FILTER_THRESH                   # Threshold[s] for retina filtering

DEBLEND_NTHRESH  32             # Number of deblending sub-thresholds
DEBLEND_MINCONT  0.005          # Minimum contrast parameter for deblending

CLEAN            Y              # Clean spurious detections? (Y or N)?
CLEAN_PARAM      1.0            # Cleaning efficiency

MASK_TYPE        CORRECT        # type of detection MASKing: can be one of
                                # NONE, BLANK or CORRECT

#-------------------------------- WEIGHTing ----------------------------------

WEIGHT_TYPE      NONE           # type of WEIGHTing: NONE, BACKGROUND,
                                # MAP_RMS, MAP_VAR or MAP_WEIGHT
RESCALE_WEIGHTS  Y              # Rescale input weights/variances (Y/N)?
WEIGHT_IMAGE     weight.fits    # weight-map filename
WEIGHT_GAIN      Y              # modulate gain (E/ADU) with weights? (Y/N)
WEIGHT_THRESH                   # weight threshold[s] for bad pixels

#-------------------------------- FLAGging -----------------------------------

FLAG_IMAGE       flag.fits      # filename for an input FLAG-image
FLAG_TYPE        OR             # flag pixel combination: OR, AND, MIN, MAX
                                # or MOST

#------------------------------ Photometry -----------------------------------

PHOT_APERTURES   5              # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,
                                # <min_radius>
PHOT_AUTOAPERS   0.0,0.0        # <estimation>,<measurement> minimum apertures
                                # for MAG_AUTO and MAG_PETRO
PHOT_FLUXFRAC    0.5            # flux fraction[s] used for FLUX_RADIUS

SATUR_LEVEL      50000.0        # level (in ADUs) at which arises saturation
SATUR_KEY        SATURATE       # keyword for saturation level (in ADUs)

MAG_ZEROPOINT    0.0            # magnitude zero-point
MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
GAIN             0.0            # detector gain in e-/ADU
GAIN_KEY         GAIN           # keyword for detector gain in e-/ADU
PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)

#------------------------- Star/Galaxy Separation ----------------------------

SEEING_FWHM      1.2            # stellar FWHM in arcsec
STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------

BACK_TYPE        AUTO           # AUTO or MANUAL
BACK_VALUE       0.0            # Default background value in MANUAL mode
BACK_SIZE        64             # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>

BACKPHOTO_TYPE   GLOBAL         # can be GLOBAL or LOCAL
BACKPHOTO_THICK  24             # thickness of the background LOCAL annulus
BACK_FILTTHRESH  0.0            # Threshold above which the background-
                                # map filter operates

#------------------------------ Check Image ----------------------------------

CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                # or APERTURES
CHECKIMAGE_NAME  check.fits     # Filename for the check-image

#--------------------- Memory (change with caution!) -------------------------

MEMORY_OBJSTACK  3000           # number of objects in stack
MEMORY_PIXSTACK  300000         # number of pixels in stack
MEMORY_BUFSIZE   1024           # number of lines in buffer

#------------------------------- ASSOCiation ---------------------------------

ASSOC_NAME       sky.list       # name of the ASCII file to ASSOCiate
ASSOC_DATA       2,3,4          # columns of the data to replicate (0=all)
ASSOC_PARAMS     2,3,4          # columns of xpos,ypos[,mag]
ASSOCCOORD_TYPE  PIXEL          # ASSOC coordinates: PIXEL or WORLD
ASSOC_RADIUS     2.0            # cross-matching radius (pixels)
ASSOC_TYPE       NEAREST        # ASSOCiation method: FIRST, NEAREST, MEAN,
                                # MAG_MEAN, SUM, MAG_SUM, MIN or MAX
ASSOCSELEC_TYPE  MATCHED        # ASSOC selection type: ALL, MATCHED or -MATCHED

#----------------------------- Miscellaneous ---------------------------------

VERBOSE_TYPE     NORMAL         # can be QUIET, NORMAL or FULL
HEADER_SUFFIX    .head          # Filename extension for additional headers
WRITE_XML        N              # Write XML file (Y/N)?
XML_NAME         sex.xml        # Filename for XML output
XSL_URL          file:///usr/share/sextractor/sextractor.xsl
                                # Filename for XSL style-sheet
NTHREADS         1              # 1 single thread

FITS_UNSIGNED    N              # Treat FITS integer values as unsigned (Y/N)?
INTERP_MAXXLAG   16             # Max. lag along X for 0-weight interpolation
INTERP_MAXYLAG   16             # Max. lag along Y for 0-weight interpolation
INTERP_TYPE      ALL            # Interpolation type: NONE, VAR_ONLY or ALL

#--------------------------- Experimental Stuff -----------------------------

PSF_NAME         default.psf    # File containing the PSF model
PSF_NMAX         1              # Max.number of PSFs fitted simultaneously
PATTERN_TYPE     RINGS-HARMONIC # can RINGS-QUADPOLE, RINGS-OCTOPOLE,
                                # RINGS-HARMONICS or GAUSS-LAGUERRE
SOM_NAME         default.som    # File containing Self-Organizing Map weights
'''
