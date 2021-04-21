# plot data coverage

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
telescopes=['DOAO','LOAO','CCA250','LSGT','MCD30INCH',
			'SOAO','MAO','KCT','SAO','ELSAUSE',
			'SQUEAN','UKIRT','itelescope']

data_directory='/data7/cschoi/IMSNGgalaxies/'
data_directory='/data7/cschoi/sngal/NGC3367'
snhost=snlist['    IMSNGname ']
snname=snlist['                           Name ']
sndiscovery=snlist['  Discovery_date ']
sntype=snlist['        Type ']
snsep=snlist[' separation ']

os.chdir(data_directory+str.strip(snhost[i]))

lines= oswalkfunc()
lines.sort()
fitslist= [s for s in lines if s.split('/')[-1][:7]=='saCalib' and s.split('/')[-1][-5:]=='.fits']
#fitslist= [s for s in lines if s.split('/')[-1][-5:]=='.fits' and s.split('/')[-1][:6] =='Calib-']

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

# DATE-OBS keyword check
flist=[]
for im in fitslist:
	print(im)
	h=fits.getheader(im)
	if ('DATE-OBS' in list(h.keys())) and (len(h['DATE-OBS'])>=19): pass
	else : flist.append(im)

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




# filename -> MJD
m.split('/')[-1].split('-')[3],m.split('/')[-1].split('-')[4]
dateobs=fits.getheader(m)['DATE-OBS']
t = Time(dateobs, format='isot', scale='utc').mjd
discovery_mjd=Time('2018-01-24T08:42:06',format='isot', scale='utc').mjd

mjds=[]
fitsset=[]
for k,s in enumerate(fitslist):
	dateobs=fits.getheader(s)['DATE-OBS']
	t = Time(dateobs, format='isot', scale='utc').mjd
	mjds.append([k,s,t])
	if (t > 58127) and (t < discovery_mjd):
		print(s,t)
		fitsset.append(s)
