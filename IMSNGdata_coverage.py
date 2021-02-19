import pandas as pd
import os
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.io.fits as fits
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
import matplotlib.dates as mdates
##############################################################################
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

'''
['#', '    IMSNGname ', ' separation ', '         snra ',
       '         sndec ', '  Discovery_date ', ' Discovery_JD ',
       '     Max_date ', '     Max_JD ', '         Host ', '        Type ',
       '    Mag ', ' Mag_Max ', '                           Name ',
       '           ra ', '          dec ', '  dist ', '  ebv ', ' muvmag ',
       ' maxaxis ', ' minaxis ', ' priority ', '         note ',
       'Unnamed: 23']
'''

snlist=pd.read_csv('/data7/cschoi/sngal/recent-sne-check/rochester-list/matched_list_2018_2020.csv') 
telescopes=['DOAO','LOAO','CCA250','LSGT','MCD30INCH','SOAO','MAO','KCT','SAO','ELSAUSE','SQUEAN','UKIRT']

data_directory='/data7/cschoi/IMSNGgalaxies/'
snhost=snlist['    IMSNGname ']   
snname=snlist['                           Name '] 
sndiscovery=snlist['  Discovery_date ']
sntype=snlist['        Type ']
snsep=snlist[' separation ']







os.chdir(data_directory+str.strip(snhost[i]))

lines= oswalkfunc() 
lines.sort() 
fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits' and s.split('/')[-1][:6] =='Calib-']  

datelist=[]
obslist=[]
filterlist=[]
mjd=[]
for n in fitslist : 
	datelist.append(n.split('/')[-1].split('-')[3])
	if  n.split('/')[-1].split('-')[3][:2] != '20' : 
		print( n)
	obslist.append(n.split('/')[-1].split('-')[1]) 
	filterlist.append(n.split('/')[-1].split('-')[5])

#	mjd.append(fits.getheader(n)['MJD'])
#   !rename SAO-STX16803 SAO_STX16803 SAO/STX16803/*/*.fits 
#   !rename T21-ceouobs.changsu T21_ceouobs.changsu itelescope/*/*.fits 
#   !rename NGC5353-190 NGC5353-20190 UKIRT/*.fits
set(obslist) 
set(filterlist)


datelist_value = [datetime.datetime.strptime(d,'%Y%m%d').date() for d in datelist]

newobslist=[]
for obs in obslist : 
	if obs[:3] == 'SOA' : newobslist.append('SOAO')
	elif obs[:3] == 'DOA' : newobslist.append('DOAO') 
	elif obs[:3] == 'MAI' : newobslist.append('MAO')
	elif obs[:3] == 'MAO' : newobslist.append('MAO')
	elif obs[:3] == 'Mai' : newobslist.append('MAO') 
	elif obs[:3] == 'SAO'  : newobslist.append('SAO')
	elif obs[:3] == 'KCT'  : newobslist.append('KCT') 
	elif obs[:3] == 'MCD'  : newobslist.append('MCD30INCH') 
	else : newobslist.append(obs)

fig, ax = plt.subplots()
y_pos = np.arange(len(telescopes))
x = 1920 / 2 / fig.dpi
y = 1080 / 2 / fig.dpi
fig.set_figwidth(x)
fig.set_figheight(y)

ax.set_yticks(y_pos)
ax.set_yticklabels(telescopes)
ax.set_xlabel('DATE')
title=str.strip(snhost[i])+' '+str.strip(snname[i])+' '+str.strip(sntype[i])+' '+str.strip(sndiscovery[i])+' '+str.strip(str(snsep[i]))+' arcmin'
ax.set_title(title)

ax=plt.gca()
formatter = mdates.DateFormatter("%Y-%m-%d")
ax.xaxis.set_major_formatter(formatter)
locator = mdates.DayLocator(interval=10)
ax.xaxis.set_major_locator(locator)
colors1 = ['C{}'.format(k) for k in range(len(telescopes))]



for m,obs in enumerate(newobslist):
	for s, tel in enumerate(telescopes):
		if obs==tel:
			ax.plot(datelist_value[m], s,'o', color=colors1[s], alpha=0.5)


discoverydate=datetime.datetime.strptime(str.strip(sndiscovery[i])[:-4],'%Y/%m/%d').date()

ax.axvline(discoverydate, label='Discovery')
ax.set_xlim(datetime.datetime.strptime('20181201','%Y%m%d').date(),datetime.datetime.strptime('20201231','%Y%m%d').date())
ax.xaxis.set_tick_params(rotation=90, labelsize=10)
plt.tight_layout()
plt.savefig('/data7/cschoi/sngal/'+str.strip(snhost[i])+'_'+str.strip(snname[i])+'.png',dpi=500,overwrite=True)
plt.close()



