# stamp image cut making 

import mumpy as np
from astropy.io import ascii
from astropy.io import fits

#input

inim='ahacfdobj.IC1222.20150630.0095.fits'

incat=inim[:-5]+'.cat'

# sextractor

# catalog reading
catdata=ascii.read(incat)
num=catdata['NUMBER']
xpix=catdata['X_IMAGE']
ypix=catdata['Y_IMAGE']

# source check

# stamp cut
def stampcut(inim,x,y) :
	os.system('imcopy '+inim+'['+str(x-25)+':'+str(x+25)+','+ str(y-25)+':'+str(y+25)+'] '+ stamp_'+str(num)+inim)

for i in range(len(num)) :
	stampcut(imim,xpix[i],ypix[i])

print 'all done'
os.system('ds9 stamp*.fits &')
