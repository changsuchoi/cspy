from multiprocessing import Process,Pool
from astropy.io import ascii
import numpy as np
import os
import sys
import glob
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
import astropy.units as u

def hotpantsrun(im, regrefim, il=0, iu=65000, tl=0, tu=65000, sigmatch=False):
	#use sigmatch=True when
	outfile='hd'+im
	convfile='hc'+im
	kernelfile='hk'+im
	# for pan starrs image subtraction set tu, tl more than 100000, -100000
	opt0=' -n t -c i'
	opt0a=' -n i -c t'
	opt1=' -il ' + str(il) +' -iu '+ str(iu)+ ' '
	opt2=' -tl ' + str(tl) +' -tu '+ str(tu)+ ' '
	opt3=' -ng 3 6 0.70 4 1.50 2 3.00'
	opt4=' -oki '+kernelfile+' '
	opt5=' -hki '
	# opt6=' -ig ' + str(ig) +' -ir '+ str(ir)+ ' ' # gain,rdnoise option
	# opt7=' -tg ' + str(tg) +' -tr '+ str(tr)+ ' ' # gain,rdnoise option
	# in the case of Sigma_image > Sigma_template, for better subtraction, you may try this option
	# FWHM = 2.355 sigma
	fwhm_im=fits.getheader(im)['FWHM_PIX']
	fwhm_reg=fits.getheader(regrefim)['FWHM_PIX']
	#com= 'hotpants -v 0 -inim '+regrefim+' -tmplim '+im+\
	#	' -outim '+outfile+opt0 +' -oci '+convfile +opt1+opt2+opt3+opt4+opt5
	com= 'hotpants -v 0 -inim '+im+' -tmplim '+regrefim+\
		' -outim '+outfile + opt0a +' -oci '+ convfile +opt1+opt2+opt3+opt4
	print(com)
	os.system(com)
	#if fwhm_im > fwhm_reg :
	#	com= 'hotpants -v 0 -inim '+regrefim+' -tmplim '+im+\
	#		' -outim '+outfile+' -n t -c i' +' -oci '+convfile +opt1+opt2+opt3+opt4+opt5
	#else:
	#	print('FWHM_INPUT',fwhm_im,'FWHM_REF',fwhm_ref)
	#	com= 'hotpants -v 0 -inim '+regrefim+' -tmplim '+im+\
	#		' -outim '+outfile+' -n t -c i' +' -oci '+convfile +opt1+opt2+opt3+opt4+opt5

def hprun(im, refim='reg_'+im, il=0, iu=65000, tl=0, tu=65000, sigmatch=False):
	print(im, refim)
	#use sigmatch=True when fwhm values are significantly different.
	#refim = 'ref.fits' : default
	outfile='hd'+im
	convfile='hc'+im
	os.system('rm '+outfile+' '+convfile)
	#kernelfile='hk'+im
	# for pan starrs image subtraction set tu, tl more than 100000, -100000
	opt0=' -n t -c i'
	opt0a=' -n i -c t'
	opt1=' -il ' + str(il) +' -iu '+ str(iu)+ ' '
	opt2=' -tl ' + str(tl) +' -tu '+ str(tu)+ ' '
	opt3=' -ng 3 6 0.70 4 1.50 2 3.00'
	#opt4=' -oki '+kernelfile+' '
	opt5=' -hki '
	# opt6=' -ig ' + str(ig) +' -ir '+ str(ir)+ ' ' # gain,rdnoise option
	# opt7=' -tg ' + str(tg) +' -tr '+ str(tr)+ ' ' # gain,rdnoise option
	# in the case of Sigma_image > Sigma_template, for better subtraction, you may try this option
	# FWHM = 2.355 sigma
	#fwhm_im=fits.getheader(im)['FWHM_PIX']
	#fwhm_reg=fits.getheader(refim)['FWHM_PIX']
	#print('input image fwhm',fwhm_im, 'ref image fwhm',fwhm_reg)
	#com= 'hotpants -v 0 -inim '+regrefim+' -tmplim '+im+\
	#	' -outim '+outfile+opt0 +' -oci '+convfile +opt1+opt2+opt3+opt4+opt5
	com= 'hotpants -v 0 -inim '+im+' -tmplim '+refim+\
		' -outim '+outfile + opt0a +' -oci '+ convfile +opt1+opt2+opt3#+opt4
	print(com)
	try:
		os.system(com)
		return 'Done'
	except : return None


	#if fwhm_im > fwhm_reg :
	#	com= 'hotpants -v 0 -inim '+regrefim+' -tmplim '+im+\
	#		' -outim '+outfile+' -n t -c i' +' -oci '+convfile +opt1+opt2+opt3+opt4+opt5
	#else:
	#	print('FWHM_INPUT',fwhm_im,'FWHM_REF',fwhm_ref)
	#	com= 'hotpants -v 0 -inim '+regrefim+' -tmplim '+im+\
	#		' -outim '+outfile+' -n t -c i' +' -oci '+convfile +opt1+opt2+opt3+opt4+opt5

	#print(com)
	#os.system(com)

#hprun(im)

'''
if sigmatch = True :
sigma_im    = fwhm_im / 2.355
sigma_ref   = fwhm_reg/ 2.355
sigma_match = np.sqrt(np.abs(sigma_im**2 - sigma_ref**2))
# E.g. -ng 3 6 0.5*Sigma_match 4 Sigma_match 2 2.0*Sigma_match
ngflag= ' -ng 3 6 '+ '%.3f'%(0.5*sigma_match) + ' 4 '+ '%.3f'%(sigma_match) +' 2 ' +'%.3f'%(2.0*sigma_match)
com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -n i -c t' + ngflag

else : com='hotpants -v 0 -inim '+regrefim+' -tmplim '+im+' -outim '+outfile+' -n t -c i'+' -oci '+convfile
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -c t'
# com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -n i -c t'
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -n i'
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -iu 60000 -tu 60000'
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -iu 60000 -tu 60000 -nrx 2 -nry 2 -nsx 10 -nsy 10 -r 5 '
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile+' -iu 60000 -tu 60000 -nrx 2 -nry 2 -nsx 10 -nsy 10 -r 5 '
#com='hotpants -v 0 -c i -n i -inim '+infile+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
#com='hotpants -v 0 -c i -inim '+infile+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
#com='hotpants -v 0 -inim '+infile+' -tmplim '+refim+' -outim '+outfile+'  -n i -c i -tg 3.6 -tr 4.0 -ig 1.5 -ir 10.0'# KPNO +Maidanak FLI
print (infile)
os.system(com)
'''

'''
for n in range(len(infile)):
print (str(n) + ' of '+ str(len(infile)))
hotpantsrun(infile[n])


cpunum=4
if __name__ == '__main__' :
p=Pool(cpunum)
p.map(hotpantsrun,infile)
'''
