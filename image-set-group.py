# make comfiles.txt for image combine
# 2020.3.26 Changsu.choi

import astropy.io.fits as fits
import glob
import os
from astropy.time import Time


#files = glob.glob('Calib*.fits')
#files.sort()

#files=oklist

def epoch_group(files):
	obsdt=[]
	for n in range(len(files)):
		header=fits.getheader(files[n])
		hobsdate=header['DATE-OBS']
		obsdt.append(hobsdate)
		print (hobsdate	+'\n')

	t = Time(obsdt, format='isot', scale='utc')
	tjd=t.mjd

	f=open('comfiles.txt','w')

	com=[]
	for i in range(len(files)) :
		#com.append(files[i])
		if i==0 : com.append(files[0])  # the 1st image

		else :                          # next images
			# if time between two exposures less than 10 min, then append new file to group	successive files will be combined together
			if (tjd[i]-tjd[i-1]) < (10/1440.) :
				com.append(files[i])               # succeeding images
			else :

				if len(com) == 1 :                 # isolated single image
					output=com[0][:-5]+'_'+str(len(com))+'_com.fits'
					print ('these ',str(len(com)),' files will be combined ',com,' output = ',output+'\n')
					# os.system('cp '+com[0]+' '+output)
					for n in com : f.write(n+'\n')
					com=[]
					com.append(files[i])

				else : 							# made set of serial images
					#f=open('tmpcom.list','w')
					for n in com : f.write(n+',')
					#f.close()
					output=com[0][:-5]+'_'+str(len(com))+'_com.fits'
					#imcombine('@tmpcom.list',output)
					#centertimeheader(com,output)
					print ('these ',str(len(com)),' files will be combined ',com,' output = ',output+'\n')
					f.write('\n')
					com=[]
					com.append(files[i])
	for n in com : f.write(n+',')
	f.write('\n')
	f.close()
	with open('comfiles.txt','r') as file_handle: lines = file_handle.read().splitlines()
	print(lines)
	return (lines)

lines=epoch_group(oklist)




## gregister

'''
for i in lines :
        if len(i.split(',')) <= 2 :
            print(i, 'single image')
        else :
            ref_image=i.split(',')[0]
            images_to_align=i.split(',')[:-1]
            print (i.split(','))
            identifications= identify_transform(ref_image, images_to_align, rad= 5, nb=500, verbose=False, visual=False)
            align_images(ref_image, identifications, iraf=True, outdir='alipy_out')
os.system('mv alipy_out/* .')

'''

#os.system('mv alipy_out/* .')

# scamp
'''
for i in range(len(imlist)) :
	scamp_net(imlist[i])


	print(i+1, 'of', str(len(imlist)) )
'''
