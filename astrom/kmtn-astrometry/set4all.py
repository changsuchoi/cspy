from astropy.io import fits
from astropy.io import ascii
import os,sys
import pp 

ppservers=()
#ncpus
ncpus=4
#make job_server
job_server = pp.Server(ncpus, ppservers=ppservers)
#job_server = pp.Server(ppservers=ppservers)


#os.system('gethead kmtc.201502*.fits ra dec filter object exptime date-obs > info.txt')
info=ascii.read('info.txt')
addlist=info['col1']
ra=info['col2']
dec=info['col3']
filters=info['col4']
obj=info['col5']
exptime=info['col6']
dateobs=info['col7']

set4=['kk','mm','nn','tt']
##  im 9-16, im 1-8, im 17-24, im 25-32  ##

def set4astrom(files) : 
	from astropy.io import fits
	from astropy.io import ascii

	print files
	num=files[14:-5]
	f=open(num+'.head','r')
	lines=f.readlines()
	f.close()

	f=open(num+'.kk.head','w')
	for line in lines[0:51] : f.write(line)
	f.close()
	data=fits.getdata(num+'.kk.fits')
	hdr=fits.getheader(num+'.kk.fits')
	hdr.fromTxtFile(num+'.kk.head')
	newfile='a'+num+'.kk.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	f=open(num+'.mm.head','w')
	for line in lines[51:102] : f.write(line)
	f.close()
	data=fits.getdata(num+'.mm.fits')
	hdr=fits.getheader(num+'.mm.fits')
	hdr.fromTxtFile(num+'.mm.head')
	newfile='a'+num+'.mm.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	f=open(num+'.nn.head','w')
	for line in lines[102:153] : f.write(line)
	f.close()
	data=fits.getdata(num+'.nn.fits')
	hdr=fits.getheader(num+'.nn.fits')
	hdr.fromTxtFile(num+'.nn.head')
	newfile='a'+num+'.nn.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	f=open(num+'.tt.head','w')
	for line in lines[153:] : f.write(line)
	f.close()
	data=fits.getdata(num+'.tt.fits')
	hdr=fits.getheader(num+'.tt.fits')
	hdr.fromTxtFile(num+'.tt.head')
	newfile='a'+num+'.tt.fits'
	fits.writeto(newfile,data,hdr,clobber=True)

	os.system('swarp -c kmtnet.swarp a'+num+'*.fits -IMAGEOUT_NAME a'+files)		
	os.system('ds9 a'+files+' &')

#for files in addlist : set4astrom(files)

#jobs = [(input, job_server.submit(set4astrom, args=(input, ),modules=("fits", ))) for input in addlist]

# for loop of function in job using multi-core
#for input, job in jobs:
#	job()


for files in addlist : set4astrom(files)
