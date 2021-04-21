#!/data7/cschoi/anaconda3/bin/python

# to fine newly discoverd sne from http://www.rochesterastronomy.org/snimages/


import requests
import re
from urllib.request import urlopen
from bs4 import BeautifulSoup
import pandas as pd
from html_table_parser import parser_functions as parser


import astropy.io.ascii as ascii
import os
import sys
from astropy.table import Table, Column
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.coordinates import ICRS
from astropy.coordinates import SkyCoord
import pandas as pd
from datetime import datetime



os.chdir('/data7/cschoi/sngal/recent-sne-check/rochester-list')
#datetime.today().strftime("%Y%m%d%H%M%S")    # YYYYmmddHHMMSS 형태의 시간 출력
#datetime.today().strftime("%Y/%m/%d %H:%M:%S")  # YYYY/mm/dd HH:MM:SS 형태의 시간 출력
today=datetime.today()
today=datetime.today().strftime("%Y%m%d %H:%M:%S")[:8]
print(today)


radius=10.0  # 30 arcmin = 0.5 deg

print ('Radius '+str(radius)+' arcmin')
# print ('Reading recentsnelist.txt file ...')
# colnames=['Ra','Dec','EarliestObs','Host','Type','Last','Max','Link','Discoverer']
# latestsnelist=pd.read_table('recentsnelist.txt')
# latestsnelist=pd.read_table('recentlist.txt')	#,names=colnames,data_start=1,guess='False')
# latestsnelist=ascii.read('recentsnelist.txt',delimiter='\t')	#,names=colnames,data_start=1,guess='False')
imsnglist=ascii.read('/data7/cschoi/IMSNG/target/alltarget.dat')



urlall="http://www.RochesterAstronomy.org/snimages/sndateall.html" # sn date all
url='http://www.rochesterastronomy.org/snimages/sndate.html' # sndate

print ('getting table data from web page from',url)
response=requests.get(url)
print('Done, table data is obtained')
soup = BeautifulSoup(response.content, 'html.parser')
tbl=soup.find_all('table')

soup.find_all('table')[1].find_all('th')
html_table = parser.make2d(tbl[1])
df=pd.DataFrame(html_table[1:], columns=html_table[0])

latestsnelist = df

print ('getting table data from web page from',urlall)
responseall=requests.get(urlall)
print('Done, table data is obtained')
soupall = BeautifulSoup(responseall.content, 'html.parser')
tblall=soupall.find_all('table')

soupall.find_all('table')[1].find_all('th')
html_tableall = parser.make2d(tblall[1])
dfall=pd.DataFrame(html_tableall[1:], columns=html_tableall[0])

latestsnelistall = dfall



if urlall=='http://www.RochesterAstronomy.org/snimages/sndateall.html' :
	print ('sneall.list file first ...')

	f=open('sneall.list','w')
	f.write('#index	ra	dec\n')
	for i in range(len(latestsnelistall)) :
		if i==0 or i==1 : pass
		else :
		#snera=latestsnelist[i][0][:11]
			snera=latestsnelistall['R.A.'][i]
			#snedec=latestsnelist[i][0][13:25]
			snedec=latestsnelistall['Decl.'][i]
			tcoord=SkyCoord(snera,snedec,frame='icrs',unit=(u.hourangle,u.deg))
			f.write(str(i)+'\t'+str(tcoord.ra.deg)+'\t'+str(tcoord.dec.deg)+'\n')

	f.close()
	print('date up to',latestsnelistall['Discovery date'][1])
	print ('sneall.list is created')

if url=='http://www.rochesterastronomy.org/snimages/sndate.html' :
	print ('sne.list file first ...')

	f=open('sne.list','w')
	f.write('#index	ra	dec\n')
	for i in range(len(latestsnelist)) :
			#snera=latestsnelist[i][0][:11]
			snera=latestsnelist['R.A.'][i]
			#snedec=latestsnelist[i][0][13:25]
			snedec=latestsnelist['Decl.'][i]
			tcoord=SkyCoord(snera,snedec,frame='icrs',unit=(u.hourangle,u.deg))
			f.write(str(i)+'\t'+str(tcoord.ra.deg)+'\t'+str(tcoord.dec.deg)+'\n')

	f.close()
	print('date up to',latestsnelist['Discovery date'][1])
	print ('sne.list is created')


'''
f=open('imsng.list','w')
f.write('#index	ra	dec\n')
for i in range(len(imsnglist)) :
	snera=imsnglist['ra'][i]
	snedec=imsnglist['dec'][i]
	imsngname=imsnglist['obj'][i]
	tcoord=SkyCoord(snera,snedec,frame='icrs',unit=(u.hourangle,u.deg))
	f.write(imsngname+'\t'+str(tcoord.ra.deg)+'\t'+str(tcoord.dec.deg)+'\n')

f.close()
print ('IMSNG list is created')
'''
print ('.\n.\n.\n.\n.\n')


print ('Matching starts...')
print ('.\n.\n.\n.\n.\n')

# matchcom="stilts tskymatch2 ifmt1=ascii find=all ifmt2=ascii in1=imsng.list in2=sne.list out=imsng-sne-matched.list ra1=ra dec1=dec ra2=ra dec2=dec error="+str(radius*60.)+" join=1and2 ofmt=ascii omode=out"
matchcom = 'java -jar /data7/cschoi/util/stilts.jar tskymatch2 ifmt1=ascii find=all ifmt2=ascii in1=imsng.list in2=sne.list out=imsng-sne-matched.list ra1=ra dec1=dec ra2=ra dec2=dec error='+str(radius*60.)+' join=1and2 ofmt=ascii omode=out'

matchcomall = 'java -jar /data7/cschoi/util/stilts.jar tskymatch2 ifmt1=ascii find=all ifmt2=ascii in1=imsng.list in2=sneall.list out=imsng-sneall-matched.list ra1=ra dec1=dec ra2=ra dec2=dec error='+str(radius*60.)+' join=1and2 ofmt=ascii omode=out'

print (matchcom)
os.system(matchcom)

print (matchcomall)
os.system(matchcomall)

matchdat=ascii.read('imsng-sne-matched.list')
matchdat.sort('index_1')

matchdatall=ascii.read('imsng-sneall-matched.list')
matchdatall.sort('index_1')

f=open('matched.list','w')
#f.write('IMSNGname	separation	snra	sndec	Discovery_date	Discovery_JD	Max_date	Max_JD	Host	Type	Mag	Mag_Max	Name	Alternate_names  \n')
f.write('IMSNGname	separation	snra	sndec	Discovery_date	Discovery_JD	Max_date	Max_JD	Host	Type	Mag	Mag_Max	Name	ra	dec	dist	ebv	muvmag	maxaxis	minaxis	priority	note\n')
for i in matchdat:
	#comment= i['index_1']+'\t\t\t'+ '% .2f'% (i['Separation']/60.)+'\t'+latestsnelist.loc[i['index_2']].values
	k=latestsnelist.loc[i['index_2']].values
	k=list(k)
	strcomment=''
	line_imsng = imsnglist[np.where(imsnglist['obj']==i['index_1'])]

	for ii in range(len(k)): k[ii] = str(k[ii])
	strcomment1='{:12s}	{:12s}	{:12s}	{:12s}	{:12s}	{:10s}	{:12s}	{:10s}	{:3s}	{:3s}	{:30s}'.format(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8],k[9],k[10])
	strcomment2='{:12s}	{:12s}	{:>5.1f}	{:4.2f}	{:2.2f}	{:4.1f}	{:4.1f}	{:2.1f}	{:12s}'\
				.format(line_imsng['ra'][0], line_imsng['dec'][0], line_imsng['dist'][0], line_imsng['ebv'][0],\
				line_imsng['muvmag'][0], line_imsng['maxaxis'][0], line_imsng['minaxis'][0], line_imsng['priority'][0],\
				line_imsng['note'][0])

	comment= '{:12s}'.format(i['index_1'])+ '\t'+ '{:04.1f}'.format(i['Separation']/60.)+'\t'+strcomment1+'\t'+strcomment2
	#vis= '{:12s}	{:12s}	{:12s}	{:5s}	{:5s}	{:5s}	{:3s}	{:2s}'.format( objname[n],ra[n],dec[n],rtp,ttp,stp,str(int(msep[n])),str(prior[n]))+'\n'
	f.write(comment+'\n')

f.close()

f=open('matchedall.list','w')
#f.write('IMSNGname	separation	snra	sndec	Discovery_date	Discovery_JD	Max_date	Max_JD	Host	Type	Mag	Mag_Max	Name	Alternate_names  \n')
f.write('IMSNGname	separation	snra	sndec	Discovery_date	Discovery_JD	Max_date	Max_JD	Host	Type	Mag	Mag_Max	Name	ra	dec	dist	ebv	muvmag	maxaxis	minaxis	priority	note\n')
for i in matchdatall:
	#comment= i['index_1']+'\t\t\t'+ '% .2f'% (i['Separation']/60.)+'\t'+latestsnelist.loc[i['index_2']].values
	k=latestsnelistall.loc[i['index_2']].values
	k=list(k)
	strcomment=''
	line_imsng = imsnglist[np.where(imsnglist['obj']==i['index_1'])]

	for ii in range(len(k)): k[ii] = str(k[ii])
	strcomment1='{:12s}	{:12s}	{:12s}	{:12s}	{:12s}	{:10s}	{:12s}	{:10s}	{:3s}	{:3s}	{:30s}'.format(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8],k[9],k[10])
	strcomment2='{:12s}	{:12s}	{:>5.1f}	{:4.2f}	{:2.2f}	{:4.1f}	{:4.1f}	{:2.1f}	{:12s}'\
				.format(line_imsng['ra'][0], line_imsng['dec'][0], line_imsng['dist'][0], line_imsng['ebv'][0],\
				line_imsng['muvmag'][0], line_imsng['maxaxis'][0], line_imsng['minaxis'][0], line_imsng['priority'][0],\
				line_imsng['note'][0])

	comment= '{:12s}'.format(i['index_1'])+ '\t'+ '{:04.1f}'.format(i['Separation']/60.)+'\t'+strcomment1+'\t'+strcomment2
	#vis= '{:12s}	{:12s}	{:12s}	{:5s}	{:5s}	{:5s}	{:3s}	{:2s}'.format( objname[n],ra[n],dec[n],rtp,ttp,stp,str(int(msep[n])),str(prior[n]))+'\n'
	f.write(comment+'\n')

f.close()


mat=ascii.read('matched.list')
mat.sort('Discovery_JD')
delim=','
mat.write('matched_order_'+today+'.list',format='ascii.fixed_width',delimiter=delim,overwrite=True)
#mat.write('matched_order_'+today+'.list.csv',format='ascii',delimiter=delim,overwrite=True)
#mat.write('matched_order.list',format='ascii.fixed_width')
print ('imsng-sne-matched_'+today+'.list & matched_'+today+'.list files are created')
print ('\a')

os.system('pluma matched_order_'+today+'.list&')

matall=ascii.read('matchedall.list')
matall.sort('Discovery_JD')
delim=','
matall.write('matched_order_'+today+'_all.list',format='ascii.fixed_width',delimiter=delim,overwrite=True)
#matall.write('matched_order_'+today+'_all.list.csv',format='ascii',delimiter=delim,overwrite=True)
#mat.write('matched_order.list',format='ascii.fixed_width')
print ('imsng-sne-matched_'+today+'_all.list & matched_'+today+'_all.list files are created')
print ('\a')

os.system('pluma matched_order_'+today+'_all.list &')



os.chdir('/data7/cschoi/')



# soup content

'''
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">

<html>
<head>
<title>Bright Supernova pages - Sorted by Date</title>
<meta content="supernova discoveries sorted by Date" name="description"/>
<meta content="supernova, date" name="keywords"/>
<meta content="dbishopx@gmail.com" name="author"/>
<meta content="http://www.RochesterAstronomy.org/snimages/sndate.html" name="identifier-url"/>
<meta content="Global" name="distribution"/>
<meta content="General" name="rating"/>
<meta content="" name="copyright"/>
<meta content="1 day" name="revisit-after"/>
<meta content="document" name="resource-type"/>
<meta content="dbishopx@gmail.com" name="reply-to"/>
<meta content="all" name="robots"/>
<meta content="en" name="language"/>
<meta charset="utf-8" content="text/html" http-equiv="CONTENT-type"/>
<link href="http://www.RochesterAstronomy.org/snimages/sndate.html" rel="canonical"/>
</head>
<p align="right"><font color="#cecece" face="Arial" size="6">
<b><i>List of supernovae sorted by Date</i></b></font></p>
<p>
This list of objects going back 18 months.
Coordinates are in J2000.0 format, the date is the
Discover date is actually the "Earliest Observation" date.  The the "Last Mag" is the last
 observed magnitude. * = 30 days old, ** = 60 days old or older.
Max Mag is the maximum observed magnitude.
Last obs is the date the object was last observed (current year assumed).
Max Date is the date of maximum observation.
For more information on any object, please click on the
 name of the object.  This will bring you directly to
 the main entry for this object.

<br/>Page re-generated at:  2020/04/16.120  UT
</p>
<table><td>
<ul><label>This list was generated from</label></ul></td></table></html>
<li><a href="../sn2018/novae.html">Nova Archive for 2018</a></li>
<li><a href="../sn2018/index.html">Supernova Archive for 2018</a></li>
<li><a href="../sn2019/novae.html">Nova Archive for 2019</a></li>
<li><a href="../sn2019/index.html">Supernova Archive for 2019</a></li>
<li><a href="../novae.html">Latest Extragalactic Novae</a></li>
<li><a href="../supernova.html">Latest Supernovae</a></li>

<td>
<ul><label>Other version of this list</label></ul></td>
<li><a href="snlocations.html">Sorted by Location (RA)</a></li>
<li><a href="sndec.html">Sorted by Location (Decl)</a></li>
<li><a href="snmag.html">Sorted by Magnitude</a></li>
<li><a href="snredshift.html">Sorted by Red Shift</a></li>
<li><a href="snhname.html">Sorted by Host name</a></li>
<li><a href="snname.html">Sorted by Name</a></li>
<li><a href="snstats.html">Statistics</a></li>



<table style='font-family:"Courier New", Courier, monospace;font-size:80%'>
<tr><th>R.A.</th><th>Decl.</th><th>Discovery date</th><th>Discovery JD</th><th>Max date</th><th>Max JD</th><th>Host</th><th>Type</th><th>Mag</th><th>Mag Max</th><th>Name</th><th>Alternate names</th></tr>
<tr><td>18:19:28.382</td><td>+55:08:53.54</td><td>2020/04/15.509</td><td>2458955.009</td><td>2020/04/15</td><td>2458954</td><td>anonymous</td><td>unk</td><td>20.1</td><td>20.1</td><td><a href="../supernova.html#2020hkp" target="_self">AT2020hkp</a></td><td>ZTF20aauxqag</td></tr>
<tr><td>16:52:15.851</td><td>+55:34:04.36</td><td>2020/04/15.465</td><td>2458954.965</td><td>2020/04/15</td><td>2458954</td><td>anonymous</td><td>unk</td><td>19.9</td><td>19.9</td><td><a href="../supernova.html#2020hko" target="_self">AT2020hko</a></td><td>ZTF19abebqgv</td></tr>
<tr><td>14:58:15.808</td><td>+33:06:51.07</td><td>2020/04/15.463</td><td>2458954.963</td><td>2020/04/15</td><td>2458954</td><td>anonymous</td><td>unk</td><td>19.7</td><td>19.7</td><td><a href="../supernova.html#2020hhh" target="_self">AT2020hhh</a></td><td>ZTF20aauvmix</td></tr>


...
...
...

<tr><td>12:35:52.30</td><td>+27:55:55.9</td><td>2012/01/11.250</td><td>2455937.750</td><td>2012/01/11</td><td>2455937</td><td>NGC 4559</td><td>LBV</td><td>16.9*</td><td>16.7</td><td><a href="../supernova.html#2016blu" target="_self">AT2016blu</a></td><td>Gaia16ada, ZTF17aaapufz, XM13MZ, PSN J12355230+2755559, TCP J12355227+2755555, TCP J12355223+2755557, TCP J12355223+2755559</td></tr>
</table>
<p>
Back to <a href="../snimages/archives.html">Archives</a>
<br/>
Back to <a href="../supernova.html">Latest Supernovae</a>
</p>





'''
