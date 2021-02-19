# 1st photometry 
# estimate FWHM value and auto mag zp
import os
import glob
import astropy.io.fits as fits
import numpy as np

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
	print('Pixel scale =', pixscale,'\"')
	return pixscale

# input files, config and params
seconfigdir ='/data7/cschoi/code/sex.config/'
seconfig    ='phot.sex'
separam     ='phot.param'
seconv      ='default.conv'
sennw       ='default.nnw' 
DETECT_MINAREA = str(5)
DETECT_THRESH  = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)
aper_low, aper_high = 0.1,10 
aper_list=np.linspace(aper_low,aper_high,31)
aper_low, aper_high = 0.1,12 
aper_list=np.linspace(aper_low,aper_high,61)

aper_input = ''
for i in aper_list: aper_input = PSCALE*(aper_input)+'{},'.format(round(i,1))
aper_input = aper_input[:-1]

# source extractor command

def secom(i):
    PSCALE=fits.getheader(i)['PSCALE']
    aper_input = ''
    for m in aper_list: aper_input = (aper_input)+'{},'.format(round(m/PSCALE,1))
    aper_input1 = aper_input[:-1]

    fn = os.path.splitext(i)[0]
    opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.se'
    opt2=' -PARAMETERS_NAME '+seconfigdir+separam +' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
    opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
    opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
    opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' + ' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits' 
    opt6=' -PHOT_APERTURES '+aper_input1
    opt7=' -PSF_NAME '+fn+'.psf '
    secommand= 'sex -c '+opt1+opt2+opt3+opt4+opt5+opt6+opt7 + i
    print(secommand)
    os.system(secommand)

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


setbl=ascii.read('test.se')
pstbl=ascii.read('../../ps1-Tonry-NGC3367.cat')
tbl=matching(setbl, pstbl, setbl['ALPHA_J2000'],setbl['DELTA_J2000'],pstbl['ra'],pstbl['dec'])
idx=np.where((tbl['FLAGS']==0) & (tbl['SNR_WIN']>10))
tbl1=tbl[idx]


#zp calculation
def zpcal(tbl,filname, magtype):
    zp=tbl1[filname]-tbl1['MAG_AUTO']
    #zp3=sigma_clipped_stats(zp, sigma=3, maxiters=10)
    zp2=sigma_clipped_stats(zp, sigma=2, maxiters=10)
    print ('zp ',zp2[0], 'zp err',zp2[2])
    sig2,sig3=[],[]
    zp3a,zp2a=[],[]
    for i in range(len(zp)): zp2a.append(zp2[0])
    for i in range(len(zp)): sig2.append(zp2[2])
    sig2=np.asarray(sig2)
    zp2a=np.asarray(zp2a)

    filtered_data=sigma_clip(zp,sigma=2,maxiters=10)
    selected, nonselected= filtered_data.mask, ~filtered_data.mask

#def zp_plot():
    xr=np.linspace(np.min(tbl1['R']), np.max(tbl1['R']), len(zp))
    plt.plot(tbl1['R'],zp,'o')
    plt.ylim(zp2[0]-1,zp2[0]+1)
    plt.xlim(np.min(tbl1['R']),np.max(tbl1['R']))
    #plt.hlines(zp3[0],xmin=12,xmax=20,color='b')
    plt.hlines(zp2[0],xmin=12,xmax=20,color='r')
#plt.fill_between(xr,zp3a+sig3,zp3a-sig3,color='b',alpha=0.5)
plt.fill_between(xr,zp2a+sig2,zp2a-sig2,color='r',alpha=0.5)
plt.plot(tbl1['R'][nonselected],zp[nonselected],'go')
plt.plot(tbl1['R'][selected],zp[selected],'ko')
plt.title(fn)
plt.ylabel('Zeropoint (AB Mag)')
plt.xlabel(filname +' Reference Mag (AB)')
plt.savefig(fn+'_zp.png')


# fits to png with regions
def fitplot(im):
    import matplotlib.pyplot as plt
    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.visualization import MinMaxInterval,ZScaleInterval,PercentileInterval
    from astropy.visualization import SqrtStretch,LinearStretch
    from astropy.visualization import ImageNormalize

    imdata,imhdr=fits.getdata(im,header=True)
    norm = ImageNormalize(imdata, interval=PercentileInterval(99.5), stretch=LinearStretch())
    wcs=WCS(imhdr)
    
    fig,ax=plt.subplot(projection=wcs)
    ax.imshow(data,cmap='gray',norm=norm,origin='lower')
    ax.invert_yaxis() 

    fig.colorbar(ax)


# result print, file, header




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




