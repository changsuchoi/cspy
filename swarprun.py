#swarp running warpper

from astropy.io import ascii
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, Column
from astropy.time import Time
from pyraf import iraf
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process,Pool

os.system('swarp -d > default.swarp')

#-------------------------------------------------------------
#center time calculation and put it to header
def centertimeheader(inim):
	obsdt=[]
	for n in range(len(inim)):
		header=fits.getheader(inim[n])
		hobsdate=header['DATE-OBS']
		obsdt.append(hobsdate)
	tt = Time(obsdt, format='isot', scale='utc')
	ttjd=tt.jd
	ttjdmean=np.mean(ttjd)
	#print (ttjd)
	#print (ttjdmean)
	ttjdmeanutc=Time(ttjdmean,format='jd',scale='utc')
	datetimestr=ttjdmeanutc.isot[:4]+ttjdmeanutc.isot[5:7]+\
				ttjdmeanutc.isot[8:10]+'-'+ttjdmeanutc.isot[11:13]+\
				ttjdmeanutc.isot[14:16]+ttjdmeanutc.isot[17:19]
	nameset=os.path.splitext(inim[0])[0].split('-')
	output=nameset[0]+'-'+nameset[1]+'-'+nameset[2]+'-'+\
			datetimestr+'-'+nameset[5]+'-'+exposuresum(inim)+\
			'_com.fits'
	#putdata,putheader=fits.getdata(putim, header=True)
	#os.system('rm '+putim)
	#putheader['DATE-OBS']=ttjdmeanutc.isot
	#fits.writeto(putim, putdata, putheader, overwrite=True)
	return ttjdmeanutc.isot, output
#---------------------------------------------------------------
def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'

def exposuresum(ims,expkey='EXPTIME'):
	exps=[fits.getheader(i)[expkey] for i in ims]
	expsum=np.sum(exps)
	return str(int(expsum))

def swarpcom(imlist):
	newdt, newname=centertimeheader(imlist)
	param_dict={
	'IMAGEOUT_NAME'     : newname,
	'COMBINE'           : 'Y',
	'COMBINE_TYPE'      : 'MEDIAN',
	'SUBTRACT_BACK'     : 'N',
	'WRITE_FILEINFO'    : 'Y',
	'BACK_TYPE'         : 'AUTO',
	'BACK_DEFAULT'      : '0.0',
	'BACK_SIZE'         : '64',
	'BACK_FILTERSIZE'   : '3',
	'CELESTIAL_TYPE'    : 'NATIVE',
	'PROJECTION_TYPE'   : 'TAN',
	'RESAMPLE'          : 'Y',
	'FSCALE_KEYWORD'    : 'FLXSCALE',
	'COPY_KEYWORDS'     : 'OBJECT,DATE-OBS,FILTER',
	'WRITE_FILEINFO'    : 'Y',
	'WRITE_XML'         : 'N',
	}
	#
	#output=os.path.splitext(inlist[0])[0]+'_'+str(len(inlist))+'_com'+os.path.splitext(inlist[0])[1]
	optstr=''
	for i in param_dict:
		#print(' -{} {}'.format(i,param_dict[i]))
		optstr += ' -{} {}'.format(i,param_dict[i])
	swarpcom='swarp -c default.swarp ' + ','.join(imlist) + optstr
	#print(swarpcom)
	os.system(swarpcom)

def radec_center(im):
	hdr = fits.getheader(im)
	#	RA, Dec center for reference catalog query
	xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
	w = WCS(im)
	racent, deccent = w.all_pix2world(xcent, ycent, 1)
	c=SkyCoord(racent,deccent,unit="deg")
	rastr=c.ra.to_string(unit=u.hourangle,sep=':')
	decstr=c.dec.to_string(unit=u.deg,sep=':')
	racent, deccent = racent.item(), deccent.item()
	return rastr,decstr,racent,deccent

salist=glob.glob('saCalib*.fits')
salist.sort()
salines=epoch_group(salist)
#salines[-1][:-1].split(',')
def swarp_epoch(salines):
	for n,i in enumerate(salines):
		print('='*60,'\n')
		print(n,'of',len(salines))
		ii=i[:-1].split(',')
		if len(ii)==1 :
			print("single image, PASS")
			pass
			#if os.path.isfile(centertimeheader(ii)[1]) : pass
			#else:
				#	os.system('cp '+ii[0]+' '+centertimeheader(ii)[1])
				#puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])
		else:
			if os.path.isfile(centertimeheader(ii)[1]) : pass
			else:
				print(len(ii),'images will be combined', centertimeheader(ii)[1])
				swarpcom(ii)
				puthdr(centertimeheader(ii)[1],'DATE-OBS', centertimeheader(ii)[0])


# swarpregister in development
'''
def swarpregister(im,refim='ref.fits'):
	rahms,decdms,rac,dec=radec_center(im)
	PSCALE=fits.getheader(im)['PSCALE']
	inputs=' '+refim+' '
	outname='regs_'+im
	swarpcom0='swarp -c default.swarp '
	opt1=' -IMAGEOUT_NAME '+ outname+' '
	opt2=' -COMBINE N '
	opt3=' -CELESTIAL_TYPE NATIVE '    # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
	opt4=' -PROJECTION_TYPE TAN '      # Any WCS projection code or NONE
	opt5=' -PROJECTION_ERR 0.001 '     # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
	opt6=' -CENTER_TYPE MANUAL '          # MANUAL, ALL or MOST
	opt7=' -CENTER '+rahms+','+decdms+' ' # Coordinates of the image center
	opt8=' -PIXELSCALE_TYPE MANUAL '      # MANUAL,FIT,MIN,MAX or MEDIAN
	opt9=' -PIXEL_SCALE '+str(PSCALE)+' ' # Pixel scale
	opt10=' -IMAGE_SIZE 0 '               # Image size (0 = AUTOMATIC)
	opt11=' -RESAMPLE Y '
	opt12=' -SUBTRACT_BACK N '
	opt13=' -DELETE_TMPFILES N '
	opt14=' -COPY_KEYWORDS OBJECT '
	opt15=' -WRITE_FILEINFO Y '
	swarpcom=swarpcom0+inputs+opt1+opt2+opt3+opt4+opt5+opt6+opt7+opt8+opt9+opt10+opt11+opt12+opt13+opt14+opt15
	print(swarpcom)

	os.system(swarpcom)
'''

# swarp -d

'''
# Default configuration file for SWarp 2.41.4
# EB 2021-02-04
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME          coadd.fits      # Output filename
WEIGHTOUT_NAME       coadd.weight.fits # Output weight-map filename

HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers

#------------------------------- Input Weights --------------------------------

WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
WEIGHT_SUFFIX          .weight.fits    # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)

#------------------------------- Co-addition ----------------------------------

COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED
                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,
                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,
                                       # AND,NAND,OR or NOR

#-------------------------------- Astrometry ----------------------------------

CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            ALL             # MANUAL, ALL or MOST
CENTER         00:00:00.0, +00:00:00.0 # Coordinates of the image center
PIXELSCALE_TYPE        MEDIAN          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            0.0             # Pixel scale
IMAGE_SIZE             0               # Image size (0 = AUTOMATIC)

#-------------------------------- Resampling ----------------------------------

RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images

RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)

FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header

GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found

#--------------------------- Background subtraction ---------------------------

SUBTRACT_BACK          Y               # Subtraction sky background (Y/N)?
                                       # (all or for each image)

BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE              128             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)

#------------------------------ Memory management -----------------------------

VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               2047            # Maximum amount of virtual memory (MB)
MEM_MAX                256             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        256             # RAM dedicated to co-addition(MB)

#------------------------------ Miscellaneous ---------------------------------

DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         N               # Write information about each input
                                       # file in the output image header?
WRITE_XML              Y               # Write XML file (Y/N)?
XML_NAME               swarp.xml       # Filename for XML output
VERBOSE_TYPE           NORMAL          # QUIET,LOG,NORMAL, or FULL

NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic

'''


# swarp -dd
'''
	# Default configuration file for SWarp 2.41.4
# EB 2021-02-04
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME          coadd.fits      # Output filename
WEIGHTOUT_NAME       coadd.weight.fits # Output weight-map filename
HEADEROUT_NAME                         # Out. header filename (overrides suffix)

HEADER_NAME                            # Header filename if suffix not used
HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers
TILE_COMPRESS          N               # Write tile compressed output image (Y/N)?

#------------------------------- Input Weights --------------------------------

WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
RESCALE_WEIGHTS        Y               # Rescale input weights/variances (Y/N)?
WEIGHT_SUFFIX          .weight.fits    # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)
WEIGHT_THRESH                         # Bad pixel weight-threshold

#------------------------------- Co-addition ----------------------------------

COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED
                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,
                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,
                                       # AND,NAND,OR or NOR
CLIP_AMPFRAC           0.3             # Fraction of flux variation allowed
                                       # with clipping
CLIP_SIGMA             4.0             # RMS error multiple variation allowed
                                       # with clipping
CLIP_WRITELOG          N               # Write output file with coordinates of
                                       # clipped pixels (Y/N)
CLIP_LOGNAME           clipped.log     # Name of output file with coordinates
                                       # of clipped pixels
BLANK_BADPIXELS        N               # Set to 0 pixels having a weight of 0

#-------------------------------- Astrometry ----------------------------------

CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            ALL             # MANUAL, ALL or MOST
CENTER         00:00:00.0, +00:00:00.0 # Coordinates of the image center
PIXELSCALE_TYPE        MEDIAN          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            0.0             # Pixel scale
IMAGE_SIZE             0               # Image size (0 = AUTOMATIC)

#-------------------------------- Resampling ----------------------------------

RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images

RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)

FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header

GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found
                                       # 0 = infinity (all or for each image)
SATLEV_KEYWORD         SATURATE        # FITS keyword for saturation level (ADU)
SATLEV_DEFAULT         50000.0         # Default saturation if no FITS keyword

#--------------------------- Background subtraction ---------------------------

SUBTRACT_BACK          Y               # Subtraction sky background (Y/N)?
                                       # (all or for each image)

BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE              128             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)
BACK_FILTTHRESH        0.0             # Threshold above which the background-
                                       # map filter operates

#------------------------------ Memory management -----------------------------

VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               2047            # Maximum amount of virtual memory (MB)
MEM_MAX                256             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        256             # RAM dedicated to co-addition(MB)

#------------------------------ Miscellaneous ---------------------------------

DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         N               # Write information about each input
                                       # file in the output image header?
WRITE_XML              Y               # Write XML file (Y/N)?
XML_NAME               swarp.xml       # Filename for XML output
XSL_URL                file:///data7/cschoi/local/share/swarp/swarp.xsl
                                       # Filename for XSL style-sheet
VERBOSE_TYPE           NORMAL          # QUIET,LOG,NORMAL, or FULL
NNODES                 1               # Number of nodes (for clusters)
NODE_INDEX             0               # Node index (for clusters)

NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic
NOPENFILES_MAX         512             # Maximum number of files opened by SWarp

'''
