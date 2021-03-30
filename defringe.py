#	FRINGE PROCESS BASED ON PROF. IM's CODE
#	WORK FOR ONLY LOAO DATA (OR YOU NEED TO MAKE 'fringe_i_ori.fits' and fringe_i.dat)
#	19.03.12	GREGORY S.H. PAEK
#	19.05.21	GREGORY S.H. PAEK
#   19.05.28    changsu fixed some minor points, enabled multi processing
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.io.fits as fits
from multiprocessing import Process,Pool
#============================================================
def fringe_cal(inim):
	data, hdr	= fits.getdata(inim, header=True)
	#subframe	= data[512:1537, 512:1537]
	#skymean, skymed, skysig		= bkgest_mask(inim)
	#subdata		= data-skymed

	dfr_list	= []
	for n in range(len(fringes)):
		'''
		fringe_b= np.median(subdata[yb1[n]:yb2[n], xb1[n]:xb2[n]])
		fringe_f= np.median(subdata[yf1[n]:yf2[n], xf1[n]:xf2[n]])
		'''
		fringe_b= np.median(data[yb1[n]:yb2[n], xb1[n]:xb2[n]])
		fringe_f= np.median(data[yf1[n]:yf2[n], xf1[n]:xf2[n]])
		dfringe	= fringe_b-fringe_f
		dfr_list.append(dfringe)
	return np.array(dfr_list)


path_base	= '/home/sonic/Research/yourpy/fringe'

masterf	= path_base+'/fringe_i.fits'
dataf, hdrf	= fits.getdata(masterf, header=True)

imlist	= glob.glob('Calib*.fits')
imlist.sort()
fringes	= ascii.read(path_base+'/pair_i.dat')

xb1, xb2= fringes['xb']-5, fringes['xb']+5
yb1, yb2= fringes['yb']-5, fringes['yb']+5

xf1, xf2= fringes['xf']-5, fringes['xf']+5
yf1, yf2= fringes['yf']-5, fringes['yf']+5


master_fri		= fringe_cal(masterf)

#for inim in imlist:
def defringe(inim):
	fscale		= np.median(fringe_cal(inim)/master_fri)
	#print(fscale)
	data, hdr	= fits.getdata(inim, header=True)
	fri_scaled	= dataf*fscale
	newim		= data-fri_scaled
	fits.writeto('Df'+inim, newim, hdr, overwrite=True)
	print (inim,fscale,'fringe corrected')
	#print(np.median(data), np.median(fri_scaled), np.median(newim))

for inim in imlist :
    defringe(inim)

'''
cpunum=4
if __name__ == '__main__' :
	p=Pool(cpunum)
	p.map(defringe,imlist)


print ('all done, check it out!')
'''



#==============================================================================
#fringe correction with source extractor
import astropy.io.fits as fits
import numpy as np
import os, sys


def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr = fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)
	comment = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'

def pixelscale(i):
	cd11 = fits.getheader(i)['CD1_1']
	cd12 = fits.getheader(i)['CD1_2']
	cd21 = fits.getheader(i)['CD2_1']
	cd22 = fits.getheader(i)['CD2_2']
	pixscale = round(np.sqrt(cd11**2 + cd21**2) * 3600, 4)
	puthdr(i, 'PSCALE', round(pixscale, 3))
	# print('Pixel scale =', pixscale,'\"')
	return pixscale

# source extractor background image
seconfigdir = '/data7/cschoi/code/cspy/sex.config/'
seconfig = 'defringe.sex'
separam = 'defringe.param'
#growthparam = 'growth.param'
#separam_noPSF = 'se1_noPSF.param'
seconv = 'default.conv'
sennw = 'default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)

def secom_defringe(im, backsize=64, backfiltersize=3, backphototype='LOCAL'):
	#PSCALE = pixelscale(im)
	fn = os.path.splitext(im)[0]
	#aper_list=[s for s in range(1,51)]
	#aper_input = ''
	#for i in aper_list: aper_input += '{},'.format(i,1)
	#aper_input = aper_input[:-1]
	opt1 = seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME ' + fn+'.sec'
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt2a = ' -PARAMETERS_NAME '+seconfigdir+separam
	#opt2b = ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	#opt2b = ' -PARAMETERS_NAME '+seconfigdir+growthparam
	opt2 = ' -FILTER_NAME '+seconfigdir+seconv + ' -STARNNW_NAME '+seconfigdir+sennw
	opt3 = ' -DETECT_MINAREA ' + DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4 = ' -DEBLEND_NTHRESH ' + DEBLEND_NTHRESH + \
	    ' -DEBLEND_MINCONT ' + DEBLEND_MINCONT
	opt5 = ' -CHECKIMAGE_TYPE BACKGROUND,BACKGROUND_RMS,-BACKGROUND '
	opt5a= ' -CHECKIMAGE_NAME '+fn+'_bg.fits,'+fn+'_bg_rms.fits,'+fn+'_sub_bg.fits '
	opt6 = ' -BACK_SIZE '+str(backsize)+' -BACK_FILTERSIZE '+str(backfiltersize)+' -BACKPHOTO_TYPE '+backphototype+' '
	#opt6 = ' -PHOT_APERTURES '+aper_input+' '
	#opt7 = ' -PSF_NAME '+fn+'.psf '
	#opt8 = ' -PIXEL_SCALE '+str(PSCALE)+' '
	secommand = 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt5a+opt6 + im
	print(secommand)
	sexout = subprocess.getoutput(secommand)
	line = [s for s in sexout.split('\n') if 'RMS' in s]
	skymed = float(line[0].split('Background:')[1].split('RMS:')[0])
	skysig= float(line[0].split('RMS:')[1].split('/')[0])
	#os.system(secommand)
	return skymed, skysig

# subtraction
