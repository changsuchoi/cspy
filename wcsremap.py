# wcsremap by Andrew Becker script

# wcsremap -template template.fits -source input.fits -outIm input_remapped.fits
# usage : python wcsremap.py 'Calibrated*.fits' ref.fits
# Changsu Choi 2017/6/28
# multi processing perfomance is not good as expected, so use single cpu run 'for' script part



import glob
import os,sys
from multiprocessing import Process, Pool
import astropy.io.fits as fits

#os.system('rm wr*.fits')

#instr  = sys.argv[1]
#inlist = glob.glob(instr).sort()
#refim  = sys.argv[2] # ref.fits

#inlist=glob.glob('Calib*.fits')
#refim='ref.fits'

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
		return 'Done'
	else :
		print('WCSTools remap and run wcsremap')
		os.system('rm '+'remap_'+refim)
		os.system('remap -v -p '+str(pixelscale(im)-0.01)+' -o remap_'+refim+' '+refim)
		print('image pixel scale < ref pixelscale', pixelscale(im), '>=', pixelscale(refim))
		print('WCSTools remap and run wcsremap')
		wcsremapstr='wcsremap -template '+im+' -source '+'remap_'+refim+' -outIm '+outim
		try:
			os.system(wcsremapstr)
			print (outim,'is created')
		except: return None
		return 'Done'


#wcstools remap
# remap -v -p 0.8 -o remap_r.fits r.fits
# template plate scale > image plate scale : safe result
# template plate scale < image plate scale : resampling or alipy

# for inim in inlist :
# 	wcsremap(im,refim='ref.fits')

'''
if __name__ == '__main__' :
	p=Pool(2)
	p.map(wcsremap,inlist)

'''
print('Done.\a')
