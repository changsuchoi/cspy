# wcsremap by Andrew Becker script

# wcsremap -template template.fits -source input.fits -outIm input_remapped.fits
# usage : python wcsremap.py 'Calibrated*.fits' ref.fits
# Changsu Choi 2017/6/28
# multi processing perfomance is not good as expected, so use single cpu run 'for' script part



import glob
import os,sys
from multiprocessing import Process, Pool

os.system('rm wr*.fits')

instr  = sys.argv[1]
inlist = glob.glob(instr).sort()
refim  = sys.argv[2] # ref.fits

inlist=glob.glob('Calib*.fits')
refim='ref.fits'


def wcsremap(inim,refim):
	#print inim
	outim='ref'+inim
	wcsremapstr='wcsremap -template '+refim+' -source '+inim+' -outIm '+outim
    # match to ref image
	wcsremapstr='wcsremap -template '+inim+' -source '+refim+' -outIm '+outim	
    # match to input image
	os.system(wcsremapstr)
	print (outim)



for inim in inlist : 
	wcsremap(inim,refim)

'''
if __name__ == '__main__' : 
	p=Pool(2)
	p.map(wcsremap,inlist)

'''
print('Done.\a')
