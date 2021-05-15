from astropy.io import fits
#============================================================
#	Function
#-------------------------------------------------------------------------#
def reproject(inim, tempim, outim):
	'''
	inim : image to transform
	tempim : template image
	'''
	from reproject import reproject_interp
	#	reference
	#	https://reproject.readthedocs.io/en/stable/
	#	inim2 --> inim1 format
	hdu1 = fits.open(inim1)[0]
	hdu2 = fits.open(inim2)[0]
	array, footprint = reproject_interp(hdu2, hdu1.header)
	fits.writeto(outim, array, hdu1.header, overwrite=True)
#-------------------------------------------------------------------------#
def hotpants(inim, refim, outim='hd.fits', convim='hc.fits'):
	'''
	inim : Science image
	refim : Reference image
	'''
	# path = os.path.dirname(inim)
	# image = os.path.basename(inim)
	# outim = f'hd{image}'
	# convim = f'hc{image}'
	com = f'hotpants -c t -n i -iu 60000 -tu 6000000000 -tl -100000 -v 0 -inim {inim} -tmplim {refim} -outim {outim} -oci {convim}'
	print(com)
	os.system(com)
#============================================================
os.system('ls *.fits')
print('='*60)
print('Input paramters')
print('-'*60)
inim1 = input('Input image\t: ')
if os.path.dirname(inim1) == '':
	inim1 = f'./{inim1}'
inim2 = input('Reference image\t: ')
if inim2 == '':
	inim2 = 'ref.fits'
#-------------------------------------------------------------------------#
#	Path & Name
#-------------------------------------------------------------------------#
path = os.path.dirname(inim1)
image1 = os.path.basename(inim1)
image2 = os.path.basename(inim2)
rpim = f'{path}/rp{image2}'
convim = f'{path}/hc{image1}' 
outim = f'{path}/hd{image1}'
#============================================================
#	Main
#-------------------------------------------------------------------------#
reproject(inim1, inim2, rpim)
hotpants(inim1, rpim, outim, convim)
