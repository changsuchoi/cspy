import os
import astropy.io.fits as fits

def oswalkfunc():
	f=open('oswalk.list','w')
	#workDIr = os.path.abspath(b'.')
	for root, dirs, files in os.walk('.'): # os.walk(".", topdown = False):
	   # all files with path names
	   for name in files:
	      #print(os.path.join(root, name))
	      f.write(os.path.join(root, name)+'\n')
	f.close()
	with open('oswalk.list','r') as file_handle: lines = file_handle.read().splitlines()
	print(len(lines),'files')
	return lines
 # lines = [line.strip() for line in file_handle]

def fnamechange(ii):
	#for CCA250
	i=ii.split('/')[-1]
	head=fits.getheader(ii)
	objname=head['OBJECT']
	dateobs=head['DATE-OBS']
	datestr=dateobs[:4]+dateobs[5:7]+dateobs[8:10]+'-'+dateobs[11:13]+dateobs[14:16]+dateobs[17:20]
	filterstr=head['FILTER']
	exptimestr=str(int(head['EXPTIME']))
	newname='Calib-CCA250-'+objname+'-'+datestr+'-'+filterstr+'-'+exptimestr+'.fits'
	print('cp '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+newname)
	os.system('cp '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+newname)

def LSGTfilechange(ii):
	# From Calib-LSGT-NGC3367-20180519-220208-g-BIN1-W-180-003.fits
	# To   Calib-LSGT-NGC3367-20180519-220208-g-180.fits
	i=ii.split('/')[-1]
	frag=i.split('-')
	frag[0]=='Calib'
#	if frag[1]=='T52' : obs='LSGT'
#	else : obs=frag[1]
	finalname='Calib-LSGT'+'-'+frag[2]+'-'+frag[3]+'-'+frag[4]+'-'+frag[5]+'-'+frag[8]+'.fits'
	os.system('mv '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+finalname)

def iTelfilechange(ii):
	# From Calib-T21-ceouobs.changsu-NGC3367-20161130-042831-R-BIN1-E-180-003.fits
	# To   Calib-T21-NGC3367-20161130-042831-R-180.fits
	i=ii.split('/')[-1]
	frag=i.split('-')
	frag[0]=='Calib'
#	if frag[1]=='T52' : obs='LSGT'
#	else : obs=frag[1]
	#finalname='Calib-'+ frag[1] +'-'+frag[2]+'-'+frag[3]+'-'+frag[4]+'-'+frag[5]+'-'+frag[8]+'.fits'
	finalname='Calib-'+ frag[1] +'-'+frag[3]+'-'+frag[4]+'-'+frag[5]+'-'+frag[6]+'-'+frag[9]+'.fits'
	os.system('mv '+ii+' '+'/'.join(ii.split('/')[:-1])+'/'+finalname)

def simplerename(ii,a,b):
	'''
	simplerename(filename, from, to)
	'''
	import os
		#i=ii.split('/')[-1]
	os.system('rename '+a+' '+b+' '+ii)

def oswalknamesep(i):
	filename=i.split('/')[-1]
	head='/'.join(i.split('/')[:-1])+'/'
	return filename, head

###########################################################################

lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits']

fitslists=[s for s in fitslist if '202012' in s] 

for i in fitslists :  
	print (i) 
	ds9list+=i+' ' 

os.system('ds9 '+ds9list+' &')


### region file ###

radius=""" 10" """
color="red"
ra=['12:56:46.737']
dec=['+21:40:55.90']
regname=['AT2018kkt']
filename=regname[0]+'.reg'
f=open(filename,'w')
head1="# Region file format: DS9 version 4.1\n"
head2="""global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"""
head3="fk5\n"
f.write(head1)
f.write(head2)
f.write(head3)

for n in range(len(ra)):
	body="circle("+str(ra[n])+","+str(dec[n])+","+radius+") # color="+color+" text={"+str(regname[n])+"}\n"	
	f.write(body)
f.close()

