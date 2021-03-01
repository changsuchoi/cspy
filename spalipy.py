# spalipy


from spalipy import Spalipy as sp
import os,sys
import sep
from astropy.io import fits
from astropy.table import Table

seconfigdir = '/data7/cschoi/code/cspy/sex.config/'
seconfig = 'se1.sex'
separam = 'spalipy.param'
#growthparam = 'growth.param'
#separam_noPSF = 'se1_noPSF.param'
seconv = 'default.conv'
sennw = 'default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)

def secom_spalipy(im):
    #PSCALE=fits.getheader(i)['PSCALE']
	PSCALE=pixelscale(im)
	#skyval, skysig,fwhm,opt_ap=opt_ap_fwhm(im)
	#aper_list,aper_list2=[3,5,7],[1.0*fwhm, 1.5*fwhm, 2.0*fwhm, 2.5*fwhm, 3.0*fwhm, opt_ap]
	#aper_input = ''
	#for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	#for i in aper_list2: aper_input += '{},'.format(round(i,1))
	#aper_input = aper_input[:-1]
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.spl'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	#opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	#opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	# 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	opt5=' -CHECKIMAGE_TYPE NONE '
	#opt6=' -PHOT_APERTURES '+aper_input+' '
	#opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	#opt9=' -SEEING_FWHM '+str(round(PSCALE*fwhm,3))+' '
	#if psf==True:
	#	secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 +im
	#else:
	#	secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 +im
	secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt8+im
	print(secommand)
	#sexout = subprocess.getoutput(secommand)
	#line = [s for s in sexout.split('\n') if 'RMS' in s]
	#skymed, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	os.system(secommand)
	return ascii.read(fn+'.spl')


# source_data = fits.getdata("source.fits")
# template_data = fits.getdata("template.fits")
source_data = fits.getdata(refim)
template_data = fits.getdata(im)

# Run sep on each set of data
# ...
# source_extracted = sep.extract(...)
# template_extracted = sep.extract(...)
#source_det = Table(source_extracted)
#template_det = Table(template_extracted)
source_det = secom_spalipy(refim)
template_det = secom_spalipy(im)
source_det.rename_column('X_IMAGE','x')
source_det.rename_column('Y_IMAGE','y')
source_det.rename_column('FLUX_AUTO','flux')
source_det.rename_column('FWHM_IMAGE','fwhm')
source_det.rename_column('FLAGS','flag')

template_det.rename_column('X_IMAGE','x')
template_det.rename_column('Y_IMAGE','y')
template_det.rename_column('FLUX_AUTO','flux')
template_det.rename_column('FWHM_IMAGE','fwhm')
template_det.rename_column('FLAGS','flag')

sp0 = sp(source_data, source_det=source_det, template_det=template_det)
sp0.align()
fits.writeto("source_aligned.fits", data=sp.aligned_data)




'''
Quick run

If you have (geometrically) well-behaved images with a significant overlap, then good results can usually be obtained with a call such as:

align-fits-simple source.fits source_aligned.fits template.fits

To take advantage of all the dials and sliders to tweak the alignment, take a look at the entire parameter descriptions via:

align-fits -h

Alternatively, one can pass lower level objects to perform an alignment interactively or within an external script, see running spalipy interactively.
Description

A source image is transformed to the pixel-coordinate system of a template image using their respective detections as tie-points.

Matching quads of detections between the two catalogues are used to match corresponding detections in the two images. An initial affine transformation is calculated from a quad match, and is applied to source image detections. Following this, cross-matching is performed within some tolerance to find corresponding detections across the image. The remaining residuals between the matched detection coordinates are used to construct 2D spline surfaces that represent the spatially-varying residuals in x and y axes. These surfaces are used to calculate the correction needed to properly register the images even in the face of non-homogeneous coordinate transformation between the images. Flux conservation is relatively robust so long as the pixel scale between source and template is the same. Proper investigation with different pixel scales has not been performed.

Note: the affine transformation uses scipy.interpolation.affine_transform which doesn't handle nans properly, therefore replace all nan values in the source image prior to running spalipy.
Examples

spalipy can be run in two modes - via the command-line scripts or interactively. The second big choice is to either provide your own detection catalogues or let spalipy perform its own detection. Each of these scenarios is shown below.
via the command-line
Using the internal detection routine

When using the internal detection routines, there are two command-line scripts: align-fits and align-fits-simple. For narrow field-of-view images without significant distortions, align-fits-simple is probabably entirely sufficient to get a good alignment. (align-fits-simple has a significantly reduced parameter list and sets some automatically, for example it will always switch off spline fitting and does not allow the user to pass existing detection catalogues.)

align-fits-simple source.fits source_aligned.fits template.fits

or

align-fits source.fits source_aligned.fits -tf template.fits

Take notice of the -tf argument in the second example, this is because align-fits offers multiple ways to provide detections, as shown in the next section.
Passing existing SExtractor detection catalogues

If one already has detection catalogues from a SExtractor run, then these can be used to save repetition.

e.g. create two SExtractor catalogues for the image:

sex -c /path/to/my/sex.config source.fits -CATALOG_NAME source.cat
sex -c /path/to/my/sex.config template.fits -CATALOG_NAME template.cat

Note: At a minimum, the SExtractor catalogues must contain the columns X_IMAGE, Y_IMAGE, FLUX, FWHM_IMAGE, FLAGS.

We must use align-fits here since align-fits-simple does not allow us to pass catalogues:

align-fits source.fits source_aligned.fits -sc source.cat -tc template.cat

Interactively
Using the internal detection routine

from astropy.io import fits
from spalipy import Spalipy

source_data = fits.getdata("source.fits")
template_data = fits.getdata("template.fits")

sp = Spalipy(source_data, template_data=template_data)
sp.align()
fits.writeto("source_aligned.fits", data=sp.aligned_data)

Passing existing detection tables

Analagously to passing SExtractor catalogues, one can pass existing astropy.Table objects when calling spalipy interactively, for examples as the output of a prior sep.extract() call.

Note: At a minimum, the detection tables must contain the columns x, y, flux, fwhm, flag.

import sep
from astropy.io import fits
from astropy.table import Table
from spalipy import Spalipy

source_data = fits.getdata("source.fits")
template_data = fits.getdata("template.fits")
# Run sep on each set of data
# ...
# source_extracted = sep.extract(...)
# template_extracted = sep.extract(...)
source_det = Table(source_extracted)
template_det = Table(template_extracted)

sp = Spalipy(source_data, source_det=source_det, template_det=template_det)
sp.align()
fits.writeto("source_aligned.fits", data=sp.aligned_data)

Logging

When running interactively, all information is output in logging. To see these one can do

import logging
logging.getLogger().setLevel(logging.INFO)  # or logging.DEBUG for more messages

prior to any of the interactive example calls.

Statistics for the transformation goodness can also be accessed via:

sp.log_transform_stats()

Parameter tuning

Several parameters should have the main focus of attention if an acceptable alignment is not being found.

    If you have a small number of detections in your image overlap then min_n_match will need to be lowered from its default of 100, but it is also worth raising n_det so that the alignment uses all of your sources. See n_det docstring on its float vs int format, but it is safe/easy to just set to some overly large value such that it won't limit your detection tables, e.g. n_det=10000.
    sub_tile at a default of 2 will effectively fit an affine transformation in each quart of the image. On extremely distorted images even this may not be enough and so it can be raised to 3 (or even 4). It is a balancing act that there must still be sufficient detections in each image in each sub tile region from which to make a fit. If you have a low number of detections, or detections are spread strongly unevenly, this should be set to 1.
    spline_order should generally only be lowered from its default of 3. Setting it to zero might actually be preferable for simple alignment tasks. Also, with a low number of detections, and particularly with regions of low number of detections, the splines may misbehave.
    max_match_dist is the tolerance in pixels when considering a source and template detection as matched after the affine transform. One may increase this in the case of poorly centred detections. Note that this has an indirect impact on min_sep (set to 2 * max_match_dist by default) - when raising max_match_dist then min_sep correspondingly increases, offering some guard against ambguous cross-matching in crowded regions. However, raising it too high may mean that too few detections pass the min_sep criterion. In crowded fields and with well-behaved detection centres, reducing max_match_dist may be advisable.
'''
