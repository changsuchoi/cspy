import astropy.io.ascii as ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.table import Table, vstack
import glob
import os,sys

datdirectory='/data7/cschoi/sngal/NGC3367/phot/'
datfiles=['SN2018kp-DOAO-R.dat']
tbl=ascii.read(datdirectory+datfiles[0])
t = Time('2018-01-24T21:57:12', format='isot', scale='utc')
discoverymjd=t.mjd
# tbl.colnames
colnames= ['FILE', 'MJD', 'OBSERVATORY', 'FILTER', 'NUMBER',
 'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', 'DELTA_J2000',
 'MAG_AUTO', 'MAGERR_AUTO', 'FLUX_AUTO', 'FLUXERR_AUTO',
 'MAG_PSF', 'MAGERR_PSF', 'FLUX_PSF', 'FLUXERR_PSF',
 'MAG_APER', 'MAG_APER_1', 'MAG_APER_2',
 'MAG_APER_3', 'MAG_APER_4', 'MAG_APER_5', 'MAG_APER_6', 'MAG_APER_7',
 'MAG_APER_8',
 'MAGERR_APER', 'MAGERR_APER_1', 'MAGERR_APER_2',
 'MAGERR_APER_3', 'MAGERR_APER_4', 'MAGERR_APER_5', 'MAGERR_APER_6', 'MAGERR_APER_7',
 'MAGERR_APER_8',
 'FLUX_APER', 'FLUX_APER_1', 'FLUX_APER_2',
 'FLUX_APER_3', 'FLUX_APER_4', 'FLUX_APER_5', 'FLUX_APER_6', 'FLUX_APER_7',
 'FLUX_APER_8',
 'FLUXERR_APER', 'FLUXERR_APER_1', 'FLUXERR_APER_2',
 'FLUXERR_APER_3', 'FLUXERR_APER_4', 'FLUXERR_APER_5', 'FLUXERR_APER_6', 'FLUXERR_APER_7',
 'FLUXERR_APER_8',
 'SNR_WIN', 'BACKGROUND', 'THRESHOLD', 'FLUX_MAX', 'MU_THRESHOLD', 'MU_MAX',
 'FLAGS', 'FWHM_IMAGE', 'FWHM_WORLD', 'ELONGATION', 'ELLIPTICITY', 'CLASS_STAR',
 'PSCALE', 'FWHM_PIX', 'SKYVAL', 'SKYSIG', 'OPT_AP', 'MATCHNUM', 'REFCAT',
 'LOWMAG', 'HIGHMAG',
 'NUM_PSF', 'ZP_PSF', 'ZPE_PSF', 'UL5_PSF',
 'NUM_AUTO', 'ZP_AUTO', 'ZPE_AUTO', 'UL5_AUTO',
 'NUM_AP3', 'ZP_AP3', 'ZPE_AP3', 'UL5_AP3',
 'NUM_AP5', 'ZP_AP5', 'ZPE_AP5', 'UL5_AP5',
 'NUM_AP7', 'ZP_AP7', 'ZPE_AP7', 'UL5_AP7',
 'NUM_F10', 'ZP_F10', 'ZPE_F10', 'UL5_F10',
 'NUM_F15', 'ZP_F15', 'ZPE_F15', 'UL5_F15',
 'NUM_F20', 'ZP_F20', 'ZPE_F20', 'UL5_F20',
 'NUM_F25', 'ZP_F25', 'ZPE_F25', 'UL5_F25',
 'NUM_F30', 'ZP_F30', 'ZPE_F30', 'UL5_F30',
 'NUM_OPTA', 'ZP_OPTA', 'ZPE_OPTA', 'UL5_OPTA',
 'sep', 'detect']


'''
fn=tbl['FILE']
# list
comcut=[s for s in list(fn) if 'com' in s and 'cut' in s]
com=[s for s in list(fn) if 'com' in s and 'cut' not in s]
single=[s for s in list(fn) if 'com' not in s and 'cut' not in s]
cut=[s for s in list(fn) if 'com' not in s and 'cut' in s]

# INDEX
comid    = np.reshape([np.where(fn==i) for i in com],-1)
comcutid = np.reshape([np.where(fn==i) for i in comcut],-1)
cutid    = np.reshape([np.where(fn==i) for i in cut],-1)
singleid = np.reshape([np.where(fn==i) for i in single],-1)

comtbl=tbl[comid]
cuttbl=tbl[cutid]
comcuttbl=tbl[comcutid]
othertbl=tbl[singleid]
'''
cctbl=vstack(comtbl,comcuttbl)
cctbl.sort()
ascii.write(cctbl, 'SN2018kp-LOAO-B.dat.com',format='commented_header',overwrite=True)

ddd=[ 'SN2018kp-SOAO-B.dat',
 'SN2018kp-LOAO-B.dat',
 'SN2018kp-DOAO-B.dat',
 'SN2018kp-MCD30INCH-B.dat',
 'SN2018kp-MAO_SNUCAM-B.dat',
 'SN2018kp-MAO_FLI-B.dat',

 'SN2018kp-MAO_FLI-R.dat',
 'SN2018kp-SOAO-R.dat',
 'SN2018kp-MAO_SNUCAM-R.dat',
 'SN2018kp-LOAO-R.dat',
 'SN2018kp-DOAO-R.dat'

 'SN2018kp-DOAO-V.dat',
 'SN2018kp-LOAO-V.dat',
 'SN2018kp-MAO_FLI-V.dat',
 'SN2018kp-CCA250-V.dat',
 'SN2018kp-MCD30INCH-V.dat',

 'SN2018kp-MCD30INCH-I.dat',
 'SN2018kp-LOAO-I.dat',
 'SN2018kp-MAO_FLI-I.dat',

 'SN2018kp-LSGT-g.dat',
 'SN2018kp-LSGT-r.dat',
 'SN2018kp-LSGT-i.dat',
]

# read cut data file make clean data file, then save to phot/cut/*.dat
# cut_SN2018kp-LOAO-R.dat
ds=glob.glob('cut*.dat')

d='cut_SN2018kp-MAO_FLI-B.dat'
magtype='MAG_APER_4'
magerrtype=magerrtype=magtype[:3]+'ERR'+magtype[3:]
zptype='ZP_F15'
zperrtype='ZPE_F15'
ultype='UL5_F15'
def dat_conversion(d,magtype='MAG_APER_4',zptype='ZP_F15',zperrtype='ZPE_F15',ultype='UL5_F15'):
    print(d,magtype, magerrtype, zptype, zperrtype, ultype)
    dd=ascii.read(d)
    t=Table()
    t.add_column(dd['MJD'])
    t.add_column(dd[magtype],name='MAG')
    t.add_column(dd[magerrtype],name='MAGERR')
    t.add_column(dd[zptype],name='ZP')
    t.add_column(dd[zperrtype],name='ZPERR')
    t.add_column(dd[ultype],name='UL5')
    t.add_column(dd['OBSERVATORY'],name='OBS')
    t.add_column(dd['FILTER'],name='FILTER')
    t.add_column(dd['detect'],name='DET')
    t.add_column(dd['FILE'],name='FILE')
    t.write('cut/'+d[4:],format='ascii.commented_header',overwrite=True)
# fil=['R']*175
# t.add_column(fil,name='FILTER')

# each filter
t=Table()
for c in ddd:
    cc=ascii.read(c)
    if '-B.dat' in c:
        t=vstack([t,cc])
t.sort('MJD','OBS')
t.write('SN2018kp-B.dat',format='ascii.commented_header',overwrite=True)
# all to one dat file
t=Table()
for c in ddd:
    cc=ascii.read(c)
    t=vstack([t,cc])
t.sort(['MJD','OBS'])
t.write('SN2018kp-all_in_one.dat',format='ascii.commented_header',overwrite=True)

dpath='/data7/cschoi/sngal/NGC3367/phot/cut/'
ddd=glob.glob(dpath+'SN2018kp*dat')
dd=ddd[0]

# marker
mk={'LOAO':'o' , 'DOAO':'^',
    'MCD30INCH':'s' ,'MAO_FLI':'P',
    'MAO_SNUCAM':'P', 'CCA250':'p',
    'LSGT':'*' , 'SOAO':'X' }
    # color
co={'B':'blue','V':'green','R':'orange','I':'red',
    'g':'lime','r':'darkorange','i':'deeppink'}

## PLOT part ##
def LCPLOT(dd,val=0):
    d=ascii.read(dd)
    print(dd)
    # plot LC
    obs=os.path.splitext(dd)[0].split('-')[-2]
    fil=os.path.splitext(dd)[0].split('-')[-1]
    plt.errorbar(d['MJD'], d['MAG']+val, yerr=d['MAGERR']/2,
                fmt=mk[obs],mec=co[fil],ecolor=co[fil],fillstyle='none',
                label=obs+' '+fil, capsize=2)
    #plot upper limit
    if 'N' in d['DET']:
        idx=np.where((d['DET']=='N') & (d['MJD']<discoverymjd))
        plt.errorbar(d[idx]['MJD'],d[idx]['UL5'],yerr=0.2,lolims=True,
                fmt='_', mec=co[fil],ecolor= co[fil])

plt.close('all')
plt.style.use('seaborn-ticks')
plt.xlabel('MJD')
plt.ylabel('MAG (AB)')
plt.title('SN 2018kp LIGHT CURVE')
plt.xlim(discoverymjd-10,58300)
plt.ylim(21.5,14)
plt.vlines(discoverymjd,14,22,linestyle='--',label='DISCOVERY',color='blue')

for dd in ddd:
    if '-B.dat' in dd: LCPLOT(dd,val=1)
    if '-V.dat' in dd: LCPLOT(dd,val=0.5)
    if '-I.dat' in dd: LCPLOT(dd,val=-1)
    if '-g.dat' in dd: LCPLOT(dd,val=0.5)
    if '-i.dat' in dd: LCPLOT(dd,val=-1)
    if '-R.dat' in dd: LCPLOT(dd,val=0)
    if '-r.dat' in dd: LCPLOT(dd,val=0)
#plt.legend()
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
#legend_elements = [Line2D([0], [0], color='b', lw=4, label='Line'),
#                   Line2D([0], [0], marker='o', color='w', label='Scatter',
#                          markerfacecolor='g', markersize=15),
#                   Patch(facecolor='orange', edgecolor='r',label='Color Patch')]
legend_elements = [Patch(facecolor='blue', edgecolor='none',label='B +1'),
                   Patch(facecolor='green', edgecolor='none',label='V +0.5'),
                   Patch(facecolor='orange', edgecolor='none',label='R'),
                   Patch(facecolor='red', edgecolor='none',label='I -1'),
                   Patch(facecolor='lime', edgecolor='none',label='g +0.5'),
                   Patch(facecolor='darkorange', edgecolor='none',label='r'),
                   Patch(facecolor='deeppink', edgecolor='none',label='i -1'),
                   Line2D([0], [0], marker='o', color='w', label='LOAO',
                        markerfacecolor='g',mec='k',fillstyle='none', markersize=10),
                   Line2D([0], [0], marker='^', color='w', label='DOAO',
                          markerfacecolor='g', mec='k',fillstyle='none', markersize=10),
                   Line2D([0], [0], marker='X', color='w', label='SOAO',
                          markerfacecolor='g', mec='k',fillstyle='none', markersize=10),
                   Line2D([0], [0], marker='s', color='w', label='30INCH',
                           markerfacecolor='g', mec='k',fillstyle='none', markersize=10),
                   Line2D([0], [0], marker='P', color='w', label='MAO',
                        markerfacecolor='g', mec='k',fillstyle='none', markersize=10),
                   Line2D([0], [0], marker='p', color='w', label='CCA250',
                        markerfacecolor='g', mec='k',fillstyle='none', markersize=10),
                   Line2D([0], [0], marker='*', color='w', label='LSGT',
                        markerfacecolor='g', mec='k',fillstyle='none', markersize=10),
                   #Line2D([0], [0], linestyle='--', color='blue', label='Discovery')
                    ]
# Create the figure
#fig, ax = plt.subplots()
#ax.legend(handles=legend_elements, loc='center')
plt.legend(handles=legend_elements,ncol=2)
plt.savefig('SN2018kp_LC_final.eps')

#========================================================================
# color plot

dp='/data7/cschoi/sngal/NGC3367/phot/cut/filters/'
dset=['SN2018kp-all_in_one.dat',  'SN2018kp-g.dat',  'SN2018kp-I.dat',  'SN2018kp-R.dat',
'SN2018kp-B.dat', 'SN2018kp-i.dat',  'SN2018kp-r.dat',  'SN2018kp-V.dat']
# Galactic Extinction
'''
#Bandpass	Central Wavelength	The Galactic extinction	Refcode of the publications
SDSS_g	0.47	0.095	2011ApJ...737..103S
SDSS_r	0.62	0.066	2011ApJ...737..103S
SDSS_i	0.75	0.049	2011ApJ...737..103S
Landolt_B	0.43	0.104	2011ApJ...737..103S
Landolt_V	0.54	0.079	2011ApJ...737..103S
Landolt_R	0.64	0.062	2011ApJ...737..103S
Landolt_I	0.80	0.043	2011ApJ...737..103S
'''
GE={'g':0.47,'r':0.62,'i':0.75,'B':0.43,'V':0.54,'R':0.64,'I':0.80}

#Galatic extinction correction, add column MAG_GEC
dat=ascii.read('SN2018kp-i.dat')
dat.add_column(dat['MAG']-GE['i'],name='MAG_GEC')
dat.write('GE_'+'SN2018kp-i.dat',format='ascii.commented_header')


# B-V color
d1=dset[4]
d2=dset[7]
df1=d1.split('.')[0].split('-')[1]
df2=d2.split('.')[0].split('-')[1]
dd1=ascii.read(dp+d1)
dd2=ascii.read(dp+d2)
dd1.add_column(dd1['MAG']-GE[df1],name='MAG_GEC')
dd2.add_column(dd2['MAG']-GE[df2],name='MAG_GEC')
vsd=vstack([dd1,dd2])
vsd.sort('MJD')
vsd0=vsd[vsd['DET']=='Y']

def near_filter_diff(vsd, row_index,std='B'):
    vsd.sort('MJD')
    vtt=vsd['MJD_GEC']
    vft=vsd['FILTER']
    tt=vsd[row_index]['MJD']
    ft=vsd[row_index]['FILTER']
    mt=vsd[row_index]['MAG_GEC']
    met=vsd[row_index]['MAGERR']
    print('This row has filter of',ft)
    #g1=vsd[vsd['FILTER']==ft]
    g2=vsd[vsd['FILTER']!=ft]
    idv=np.where(abs(g2['MJD']-tt)==np.min(abs(g2['MJD']-tt)) )
    idv0=np.where(vsd==g2[idv])
    fs=g2[idv]['FILTER']
    ts=g2[idv]['MJD']
    ms=g2[idv]['MAG_GEC']
    mes=g2[idv]['MAGERR']
    print('Filter',ft.item(),fs.item())
    if ft==std:
        mjddiff=ts-tt
        cval=mt-ms
        cerr=np.sqrt(met**2+mes**2)
    else :
        print(fs.item(),'-',ft.item())
        mjddiff=ts-tt
        cval=ms-mt
        cerr=np.sqrt(met**2+mes**2)
    return idv0[0].item(),round(((ts+tt)/2).item(),7), round(mjddiff.item()*24,7), round(cval.item(),3), round(cerr.item(),3)


cmjd,cmjddiff,cmag,cerr=[],[],[],[]
for i,n in enumerate(vsd0):
    print(i,n['MJD'],n['FILTER'])
    cres=near_filter_diff(vsd0, i,std='B')
    cmjd.append(cres[1])
    cmjddiff.append(cres[2])
    cmag.append(cres[3])
    cerr.append(cres[4])

t1=Table()
t1.add_column(cmjd,name='MJD')
t1.add_column(cmjddiff,name='MJD_DIFF')
t1.add_column(cmag,name='MAG')
t1.add_column(cerr,name='MAGERR')
t1.sort('MJD')
t1.write('B_V_color.dat',format='ascii.commented_header',overwrite=True)
