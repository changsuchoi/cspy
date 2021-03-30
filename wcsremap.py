# wcsremap by Andrew Becker script

# wcsremap -template template.fits -source input.fits -outIm input_remapped.fits
# usage : python wcsremap.py 'Calibrated*.fits' ref.fits
# Changsu Choi 2017/6/28
# multi processing perfomance is not good as expected, so use single cpu run 'for' script part



import glob
import os,sys
from multiprocessing import Process, Pool
import astropy.io.fits as fits
import numpy as np
import alipy
import glob

refim='ref.fits'

def identify_transform(ref_image, images_to_align,rad= 5, nb=500,verbose=False, visual=False):
    identifications = alipy.ident.run(ref_image, images_to_align,r=rad, n=nb, visu=visual)
    if verbose:
        for id in identifications:  # list of the same length as images_to_align.
            if id.ok is True:  # i.e., if it worked
                print("%20s : %20s, flux ratio %.2f" % (id.ukn.name,
                                                        id.trans,
                                                        id.medfluxratio))

            else:
                print( "%20s : no transformation found !" % (id.ukn.name))
    return identifications

# identifications= identify_transform(ref_image, images_to_align, rad= 5, nb=500, verbose=False, visual=False)

def align_images(ref_image, identifications, iraf=True, outdir='alipy_out'):
    outputshape = alipy.align.shape(ref_image)
    for id in identifications:
            if id.ok is True:
                if iraf:
                    alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars,
                                          id.refmatchstars, shape=outputshape,
                                          makepng=False, outdir=outdir)
                else:
                    alipy.align.affineremap(id.ukn.filepath, id.trans,
                                            shape=outputshape, makepng=True,
                                            outdir=outdir)
# register ref to im
def alipyrun(im,ref_image):
	identifications= identify_transform(im, ref_image,
					rad= 5, nb=500, verbose=False, visual=False)
	align_images(im, identifications, iraf=True, outdir='alipy_out')

def alipy_in2ref(ref_image,images_to_align):
    '''
    images_to_align should be a list!
    '''
    identifications= identify_transform(ref_image, images_to_align, rad= 5, nb=500, verbose=False, visual=False)
    align_images(ref_image, identifications, iraf=True, outdir='alipy_out')
    for i in images_to_align :
        newname = 'reg_'+os.path.splitext(i)[0]+os.path.splitext(i)[1]
        os.system('mv alipy_out/'+os.path.splitext(i)[0]+'_gregister'+os.path.splitext(i)[1]+' '+newname)

# match ref image to input image
# ref image will be aligned to input image, no change in input image
# im='Calib-LSGT-NGC3367-20190630-084141-r-180.fits'
# ref_image=['ref.fits']
# output name = reg_Calib-LSGT-NGC3367-20190630-084141-r-180.fits

def alipy_ref2in(im,head='reg_',ref_image=[refim]):
	'''
	ref_image should be a list!
	like ref_image=['ref.fits']
	'''
	newname = head+os.path.splitext(im)[0]+os.path.splitext(im)[1]
	print(im, newname)
	if os.path.isfile(newname):
		os.system('rm '+newname)
		#return 'pass'
	alipyrun(im,ref_image)
	s=os.system('mv alipy_out/'+os.path.splitext(ref_image[0])[0]+
				'_gregister'+os.path.splitext(ref_image[0])[1]+' '+newname)
	if s==0 :
		print('Well done')
	else:
		print('Alipy did not worked')
		print(im, 'stopped')
		return None
	return 'Done'

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
	#pixscale=round(np.sqrt(cd11**2 + cd21**2) *3600 ,4)
	pixscale=np.sqrt(cd11**2 + cd21**2) *3600
	puthdr(i,'PSCALE',pixscale)
	#print('Pixel scale =', pixscale,'\"')
	return pixscale

def remap2min(imlist, refim='ref.fits'):
	os.system('rm remap_'+refim)
	import astropy.io.fits as fits
	import numpy as np
	pslist=[fits.getheader(s)['PSCALE'] for s in imlist]
	print('the least pixelscale =',np.min(pslist))
	os.system('remap -v -p '+str(np.min(pslist)-0.001)+' -o remap_'+refim+' '+refim)
	print('remap_'+refim, pixelscale('remap_'+refim))

def remap2max(imlist, refim='ref.fits'):
	os.system('rm remap_'+refim)
	import astropy.io.fits as fits
	import numpy as np
	pslist=[fits.getheader(s)['PSCALE'] for s in imlist]
	print('the largest pixelscale =',np.max(pslist))
	os.system('remap -v -p '+str(np.max(pslist)+0.001)+' -o remap_'+refim+' '+refim)
	print('remap_'+refim, pixelscale('remap_'+refim))

#remap reference image to input image

def remap_ref2in(im, refim='ref.fits'):
	print ('='*60,'\n','refim', 'to', im)
	os.system('rm rem_'+im)
	remapcom='remap -v '+' -f '+im+' -p '+str(pixelscale(im))+' -o '+'rem_'+im+' '+refim
	os.system(remapcom)
'''
	REMAP WCSTools 3.9.6, 14 May 2020, Jessica Mink (jmink@cfa.harvard.edu)
ERROR: h Argument not recognized
Remap FITS or IRAF images into single FITS image using WCS
Usage: remap [-v][-f WCSfile][-a rot][[-b][-j] ra dec][-i bits][-l num] file1.fit file2.fit ... filen.fit
  or : remap [-v][-f WCSfile][-a rot][[-b][-j] ra dec][-i bits][-l num] @filelist
  -a: Output rotation angle in degrees (default 0)
  -b ra dec: Output center in B1950 (FK4) RA and Dec
  -e long lat: Output center in ecliptic longitude and latitude
  -f file: Use WCS from this file as output WCS
  -g long lat: Output center in galactic longitude and latitude
  -i num: Number of bits per output pixel (default is input)
  -j ra dec: center in J2000 (FK5) RA and Dec
  -l num: Log every num rows of output image
  -m mode: c closest pixel (more to come)
  -n num: integer pixel value for blank pixel
  -o name: Name for output image
  -p secpix: Output plate scale in arcsec/pixel (default =input)
  -s: Set BZERO and BSCALE in output file from input file
  -t: Number of samples per linear output pixel
  -u: Delete distortion keywords from output file
  -v: Verbose
  -w type: Output WCS type (input is default)
  -x x y: Output image reference X and Y coordinates (default is center)
  -y nx ny: Output image dimensions (default is first input image)
'''
# remap input images to reference image
def remap_in2ref(files,ref='ref.fits'):
	remapcom='remap -v -f '+ref+' -p '+str(pixelscale(ref))+' -o rem_'+files+' '+files
	os.system(remapcom)
	inheader=fits.getheader(files)
	output='rem'+files
	outdata,outheader=fits.getdata(output, header=True)
	outheader['DATE-OBS'] = inheader['DATE-OBS']
	fits.writeto(output, outdata, outheader, clobber=True)

def wcsremap(im,refim='ref.fits'):
	#print inim
	outim='rew_'+im
	if os.path.isfile(outim):
		os.system('rm '+outim)
		#print(outim,'exist, Pass')
		#return 'Pass'
		#wcsremapstr='wcsremap -template '+refim+' -source '+im+' -outIm '+outim
		# match to ref image
	if pixelscale(im) >= pixelscale(refim):
		wcsremapstr='wcsremap -template '+im+' -source '+refim+' -outIm '+outim
		print('normal wcsremap run')
		try:
			os.system(wcsremapstr)
			print (outim,'is created')
		except: return None
		return 'Normal'
	else :
		print('WCSTools remap and run wcsremap')
		# os.system('rm '+'remap_'+refim)
		if os.path.isfile('remap_'+refim):
			if pixelscale(im) >= pixelscale('remap_'+refim):
				print('image pixel scale < ref pixelscale', pixelscale(im), '>=', pixelscale(refim))
				print('WCSTools remap and run wcsremap', pixelscale(im),pixelscale('remap_'+refim))
				wcsremapstr='wcsremap -template '+im+' -source '+'remap_'+refim+' -outIm '+outim
				try:
					os.system(wcsremapstr)
					print (outim,'is created')
					puthdr(outim,'COMMENT','Remapped reference used')
				except: return None

		else:
			print('make least pixelscale remap_ref.fits first!')
			return None

def wcsremap_im2ref(im,refim='ref.fits'):
	#print inim
	outim='rew_'+im
	if os.path.isfile(outim):
		os.system('rm '+outim)
		#print(outim,'exist, Pass')
		#return 'Pass'
		#wcsremapstr='wcsremap -template '+refim+' -source '+im+' -outIm '+outim
		# match to ref image
	if pixelscale(refim) >= pixelscale(im):
		wcsremapstr='wcsremap -template '+refim+' -source '+im+' -outIm '+outim
		print('normal wcsremap run')
		try:
			os.system(wcsremapstr)
			print (outim,'is created')
		except: return None
		return 'Normal'
	else :
		print('WCSTools remap and run wcsremap')
		# os.system('rm '+'remap_'+refim)
		if os.path.isfile('remap_'+refim):
			if pixelscale('remap_'+refim) >= pixelscale(im):
				print('image pixel scale < ref pixelscale', pixelscale(refim), '>=', pixelscale(im))
				print('WCSTools remap and run wcsremap', 'input',pixelscale(im),'remap ref',pixelscale('remap_'+refim))
				wcsremapstr='wcsremap -template '+'remap_'+refim+' -source '+im+' -outIm '+outim
				try:
					os.system(wcsremapstr)
					print (outim,'is created')
					puthdr(outim,'COMMENT','Remapped reference used')
				except: return None

		else:
			print('make least pixelscale remap_ref.fits first!')
			return None


def wcsremap_alipyref2in(im, refim='ref.fits', gregister=True):
	#print inim
	outim='rew_'+im
	print('='*60,'\n')
	if os.path.isfile(outim):
		os.system('rm '+outim)
		#print(outim,'exist, Pass')
		#return 'Pass'
		#wcsremapstr='wcsremap -template '+refim+' -source '+im+' -outIm '+outim
		# match to ref image
	if pixelscale(im) >= pixelscale(refim):
		wcsremapstr='wcsremap -template '+im+' -source '+refim+' -outIm '+outim
		print('normal wcsremap run')
		try:
			os.system(wcsremapstr)
			print (outim,'is created')
			puthdr(outim,'COMMENT','wcsremap, normal')
		except: return None
		return 'Normal'
	else :
		if gregister!=True:
			print('WCSTools remap and run wcsremap')
			# os.system('rm '+'remap_'+refim)
			if os.path.isfile('remap_'+refim):
				if pixelscale(im) >= pixelscale('remap_'+refim):
					print('image pixel scale < ref pixelscale', pixelscale(im), '>=', pixelscale(refim))
					print('WCSTools remap and run wcsremap', pixelscale(im),pixelscale('remap_'+refim))
					wcsremapstr='wcsremap -template '+im+' -source '+'remap_'+refim+' -outIm '+outim
					try:
						os.system(wcsremapstr)
						print (outim,'is created')
						puthdr(outim,'COMMENT','wcsremap, Remapped reference used')
					except: return None
			else:
				print('make least pixelscale remap_ref.fits first!')
				return None
		else:
			print('Alipy_ref2in will be used')
			alipy_ref2in(im,head='rew_',ref_image=[refim])
			puthdr(outim,'COMMENT','alipy_ref2in, with normal reference')

#wcstools remap
# remap -v -p 0.8 -o remap_r.fits r.fits
# template plate scale > image plate scale : safe result
# template plate scale < image plate scale : resampling or alipy

# for inim in inlist :
# 	wcsremap(im,refim='ref.fits')
def remap_pix(im,pix):
	print('='*60,'\n',im,'remapping','\n')
	remap='remap -v -p '+str(pix)+' -o rem_'+im+' '+im
	os.system(remap)



'''
if __name__ == '__main__' :
	p=Pool(2)
	p.map(wcsremap,inlist)

'''
print('Done.\a')
