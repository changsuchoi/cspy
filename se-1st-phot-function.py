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
from scipy.interpolate import UnivariateSpline

def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'

def pixelscale(i):
	cd11 = fits.getheader(i)['CD1_1']
	cd12 = fits.getheader(i)['CD1_2']
	cd21 = fits.getheader(i)['CD2_1']
	cd22 = fits.getheader(i)['CD2_2']
	pixscale=round(np.sqrt(cd11**2 + cd21**2) *3600 ,4)
	puthdr(i,'PSCALE',round(pixscale,3))
	#print('Pixel scale =', pixscale,'\"')
	return pixscale

def fwhm_img(im,mtbl1):
	fwhm_img=sigma_clipped_stats(mtbl1['FWHM_IMAGE'], sigma=3, maxiters=10)
	filtered_data=sigma_clip(mtbl1['FWHM_IMAGE'],sigma=3,maxiters=10)
	selected, nonselected= ~filtered_data.mask, filtered_data.mask
	print('FWHM_IMAGE','{}'.format(round(fwhm_img[0],3)),
		len(mtbl1[selected]),'stars from',len(mtbl1))
	puthdr(im, 'FWHM_PIX', round(fwhm_img[0],3), hdrcomment='FWHM PIXEL VALUE')
	return round(fwhm_img[0],3)

# input files, config and params
seconfigdir ='/data7/cschoi/code/cspy/sex.config/'
seconfig    ='se1.sex'
separam     ='se1.param'
separam_noPSF = 'se1_noPSF.param'
growthparam = 'growth.param'
seconv      ='default.conv'
sennw       ='default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH  = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)
lowmag=14
highmag=19
filname,filerr='R','Rerr'
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
magtype=magtypes[0]
refcat='../../ps1-Tonry-NGC3367.cat'
# source extractor command

def secom(im,psf=False):
    #PSCALE=fits.getheader(i)['PSCALE']
	PSCALE=pixelscale(im)
	skyval, skysig,fwhm,opt_ap=opt_ap_fwhm(im)
	aper_list,aper_list2=[3,5,7],[1.0*fwhm, 1.5*fwhm, 2.0*fwhm, 2.5*fwhm, 3.0*fwhm, opt_ap]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	for i in aper_list2: aper_input += '{},'.format(round(i,1))
	aper_input = aper_input[:-1]
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.se1'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	opt5=' -CHECKIMAGE_TYPE NONE '
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	opt9=' -SEEING_FWHM '+str(round(PSCALE*fwhm,3))+' '
	if psf==True:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 +im
	else:
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 +im
	print(secommand)
	#sexout = subprocess.getoutput(secommand)
	#line = [s for s in sexout.split('\n') if 'RMS' in s]
	#skymed, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	os.system(secommand)
	return skyval, skysig

#macthing
def matching(intbl, reftbl, inra, indec, refra, refdec, sep=2.0):
    """
    MATCHING TWO CATALOG WITH RA, Dec COORD. WITH python
    INPUT   :   SE catalog, SDSS catalog file name, sepertation [arcsec]
    OUTPUT  :   MATCED CATALOG FILE & TABLE
    """
    import numpy as np
    import astropy.units as u
    from astropy.table import Table, Column
    from astropy.coordinates import SkyCoord
    from astropy.io import ascii
    incoord     = SkyCoord(inra, indec, unit=(u.deg, u.deg))
    refcoord    = SkyCoord(refra, refdec, unit=(u.deg, u.deg))
    #   INDEX FOR REF.TABLE
    indx, d2d, d3d  = incoord.match_to_catalog_sky(refcoord)
    mreftbl         = reftbl[indx]
    mreftbl['sep']  = d2d
    mergetbl        = intbl
    for col in mreftbl.colnames:
        mergetbl[col]   = mreftbl[col]
    indx_sep        = np.where(mergetbl['sep']*3600.<sep)
    mtbl            = mergetbl[indx_sep]
    #mtbl.write(mergename, format='ascii', overwrite=True)
    return mtbl

def limitmag(N, zp, aper, skysigma):			# 3? 5?, zp, diameter [pixel], skysigma
	import numpy as np
	R           = float(aper)/2.				# to radius
	braket      = N*skysigma*np.sqrt(np.pi*(R**2))
	upperlimit  = float(zp)-2.5*np.log10(braket)
	return round(upperlimit, 3)

def starcut(mtbl,lowmag=lowmag,highmag=highmag,filname=filname,magtype=magtype):
	idx=np.where( (mtbl['SNR_WIN'] >20) &
				(mtbl['FLAGS'] == 0) &
				(mtbl[filname] < highmag) &
				(mtbl[filname] > lowmag) &
				(mtbl[magtype[:3]+'ERR'+magtype[3:]]<0.1)		)
	return mtbl[idx]

#zp calculation
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']

def zpcal(mtbl1,filname, magtype):
    zp=mtbl1[filname]-mtbl1[magtype]
    #zp3=sigma_clipped_stats(zp, sigma=3, maxiters=10)
    zp2=sigma_clipped_stats(zp, sigma=2, maxiters=10)
    print ('zp ',zp2[0], 'zp err',zp2[2])
    filtered_data=sigma_clip(zp,sigma=2,maxiters=10)
    selected, nonselected= ~filtered_data.mask, filtered_data.mask
    zperrp=np.sqrt( np.sum(mtbl1[filerr][selected]**2 + \
						mtbl1[magtype[:3]+'ERR'+magtype[3:]][selected]**2)\
						/ len(mtbl1) )
    print(magtype, 'zp', '{},'.format(round(zp2[0],3)),
	 	'zperr', '{},'.format(round(zperrp,3)),
		len(mtbl1[selected]),'stars from',len(mtbl1))
    return zp2, selected,zperrp

def fwhm_img(im,mtbl1):
	fwhm_img=sigma_clipped_stats(mtbl1['FWHM_IMAGE'], sigma=3, maxiters=10)
	filtered_data=sigma_clip(mtbl1['FWHM_IMAGE'],sigma=3,maxiters=10)
	selected, nonselected= ~filtered_data.mask, filtered_data.mask
	print('FWHM_IMAGE','{}'.format(round(fwhm_img[0],3)),
		len(mtbl1[selected]),'stars from',len(mtbl1))
	puthdr(im, 'FWHM_PIX', round(fwhm_img[0],3), hdrcomment='FWHM PIXEL VALUE')
	return fwhm_img[0]

'''
def fwhm_wcs(mtbl1):
	fwhm_wcs=sigma_clipped_stats(mtbl1['FWHM_WORLD'], sigma=3, maxiters=10)
	filtered_data=sigma_clip(mtbl1['FWHM_WORLD'],sigma=3,maxiters=10)
	selected, nonselected= ~filtered_data.mask, filtered_data.mask
	print('FWHM_WORLD','{},'.format(round(fwhm_wcs[0]*3600,3)),
		len(mtbl1[selected]),'stars from',len(mtbl1))
	return fwhm_wcs[0]*3600
'''
# 5sigma detection limit estimate for MAG_AUTO, MAG_PSF
# error fitting polinomial
def UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2):
	fn=os.path.splitext(im)[0]
	from astropy.modeling import models, fitting
	import numpy as np
	magerrtype = magtype[:3]+'ERR'+magtype[3:]
	x,y = mtbl1[magtype],mtbl1[magerrtype]
	#x,y=setbl['MAG_AUTO'],setbl['MAGERR_AUTO']
	x,y=x[np.where(y<1)],y[np.where(y<1)]
	#fit_init=models.Polynomial1D(7)
	fit_init=models.Exponential1D()
	fit_t=fitting.LevMarLSQFitter()
	#fit_t=fitting.LinearLSQFitter()
	t=fit_t(fit_init,x,y)
	plt.plot(setbl[magtype]+zp2[0],setbl[magerrtype],'ro')
	plt.plot(mtbl[magtype]+zp2[0],mtbl[magerrtype],'ko')
	plt.plot(x+zp2[0],y,'bo')
	plt.xlabel('Mag')
	plt.ylabel('Error')
	plt.xlim(10,25)
	plt.ylim(-0.1,0.5)
	xp=np.linspace(-20,0,20001)
	plt.plot(xp+zp2[0],t(xp),'--')
	plt.hlines(0.2,10,25)
	idx_min=np.where(np.abs(t(xp)-0.198) == np.min(np.abs(t(xp)-0.198)))
	xp[idx_min]+zp2[0] #
	plt.vlines(xp[idx_min]+zp2[0],-0.1,0.5)
	# result print, file, header
	plt.text(12,0.4,'5 sigma Detection Limit error=0.198')
	plt.text(12,0.3,'5sig_UL = '+'{}'.format(round((xp[idx_min]+zp2[0])[0],3)))
	plt.title(fn+ ' '+ magtype+' '+'5 sig Detection Limit')
	plt.savefig(fn+'_'+magtype+'_'+'5sigUL.png')
	plt.close()
	return round((xp[idx_min]+zp2[0])[0],3)


def zp_plot(mtbl1, zp2, selected, magtype, im, filname=filname, filerr=filerr):
	fn=os.path.splitext(im)[0]
	zp=mtbl1[filname]-mtbl1[magtype]
	xr=np.linspace(np.min(mtbl1[filname]), np.max(mtbl1[filname]), len(zp))
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	zperrp=np.sqrt( np.sum(mtbl1[filerr][selected]**2 + \
						mtbl1[magtype[:3]+'ERR'+magtype[3:]][selected]**2)\
						/ len(mtbl1) )
	plt.plot(mtbl1[filname],zp,'o')
	plt.errorbar(mtbl1[filname],zp,yerr=zperrp,fmt='o',capsize=1)
	plt.ylim(zp2[0]-2,zp2[0]+2)
	plt.xlim(np.min(mtbl1[filname]),np.max(mtbl1[filname]))
	#plt.hlines(zp3[0],xmin=12,xmax=20,color='b')
	plt.hlines(zp2[0],xmin=12,xmax=20,color='r')
	sig2=np.ones(len(mtbl1))*zp2[2]
	zp2a=np.ones(len(mtbl1))*zp2[0]
	#plt.fill_between(xr,zp3a+sig3,zp3a-sig3,color='b',alpha=0.5)
	plt.fill_between(xr,zp2a+sig2,zp2a-sig2,color='r',alpha=0.5)
	plt.plot(mtbl1[filname][~selected],zp[~selected],'go')
	plt.plot(mtbl1[filname][selected],zp[selected],'ko')
	plt.title(fn+', '+filname+' '+magtype)
	plt.ylabel('Zeropoint (AB Mag)')
	plt.xlabel(filname +' Reference Mag (AB)')
	plt.savefig(fn+'_'+filname+'_'+magtype+'_zp.png')
	plt.close()

# fits to png with regions
def fitplot(im, mtbl1, magtype, selected):
	import matplotlib.pyplot as plt
	from astropy.wcs import WCS
	from astropy.io import fits
	from astropy.visualization import MinMaxInterval,ZScaleInterval,PercentileInterval
	from astropy.visualization import SqrtStretch,LinearStretch
	from astropy.visualization import ImageNormalize
	imdata,imhdr=fits.getdata(im,header=True)
	norm = ImageNormalize(imdata,
			interval=PercentileInterval(99.5),
			stretch=LinearStretch()
							)
	wcs=WCS(imhdr)
	fig,ax=plt.subplots()
	ax=plt.subplot(projection=wcs)
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')
	#ax.invert_xaxis()
	ax.set_title(im+' '+magtype)
	ax.scatter(mtbl1[selected]['ALPHA_J2000'],mtbl1[selected]['DELTA_J2000'],
		transform=ax.get_transform('fk5'),s=20, edgecolor='green',facecolor='none')
	ax.scatter(mtbl1[~selected]['ALPHA_J2000'],mtbl1[~selected]['DELTA_J2000'],
		transform=ax.get_transform('fk5'),s=20, edgecolor='red',facecolor='none')
	img=ax.imshow(imdata,cmap='gray',norm=norm,origin='lower')
	ax.invert_yaxis()
	fig.colorbar(img)
	plt.savefig(os.path.splitext(im)[0]+'_'+magtype+'_FOV.png')
	plt.close()














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
