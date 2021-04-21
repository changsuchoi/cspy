import astropy.io.ascii as ascii
import numpy as np
import matplotlib.pyplot as plt


datdirectory='/data7/cschoi/sngal/NGC3367/phot/'
datfiles=[
'SN2018kp-CCA250-V.dat','SN2018kp-LOAO-V.dat',        'SN2018kp-MCD30INCH-B.dat',
'SN2018kp-DOAO-B.dat'  ,'SN2018kp-LSGT-g.dat',        'SN2018kp-MCD30INCH-I.dat',
'SN2018kp-DOAO-R.dat'  ,'SN2018kp-LSGT-i.dat',        'SN2018kp-MCD30INCH-R.dat',
'SN2018kp-DOAO-V.dat'  ,'SN2018kp-LSGT-r.dat',        'SN2018kp-MCD30INCH-V.dat',
'SN2018kp-LOAO-B.dat'  ,'SN2018kp-MAO_FLI-R.dat',     'SN2018kp-SOAO-B.dat',
'SN2018kp-LOAO-I.dat'  ,'SN2018kp-MAO_SNUCAM-B.dat',  'SN2018kp-SOAO-R.dat',
'SN2018kp-LOAO-R.dat'  ,'SN2018kp-MAO_SNUCAM-R.dat']

# tbl.colnames
colnames= ['FILE', 'MJD', 'OBSERVATORY', 'FILTER', 'NUMBER',
 'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', 'DELTA_J2000',
 'MAG_AUTO', 'MAGERR_AUTO', 'FLUX_AUTO', 'FLUXERR_AUTO',
 'MAG_PSF', 'MAGERR_PSF', 'FLUX_PSF', 'FLUXERR_PSF',
 'MAG_APER', 'MAG_APER_1', 'MAG_APER_2',
 'MAG_APER_3', 'MAG_APER_4', 'MAG_APER_5', 'MAG_APER_6', 'MAG_APER_7', 'MAG_APER_8',
 'MAGERR_APER', 'MAGERR_APER_1', 'MAGERR_APER_2',
 'MAGERR_APER_3', 'MAGERR_APER_4', 'MAGERR_APER_5', 'MAGERR_APER_6', 'MAGERR_APER_7', 'MAGERR_APER_8',
 'FLUX_APER', 'FLUX_APER_1', 'FLUX_APER_2',
 'FLUX_APER_3', 'FLUX_APER_4', 'FLUX_APER_5', 'FLUX_APER_6', 'FLUX_APER_7', 'FLUX_APER_8',
 'FLUXERR_APER', 'FLUXERR_APER_1', 'FLUXERR_APER_2',
 'FLUXERR_APER_3', 'FLUXERR_APER_4', 'FLUXERR_APER_5', 'FLUXERR_APER_6', 'FLUXERR_APER_7', 'FLUXERR_APER_8',
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
magtypes=['MAG_AUTO', 'MAG_PSF','MAG_APER','MAG_APER_1','MAG_APER_2',
	'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7','MAG_APER_8']
zptypes  = ['ZP_PSF','ZP_AUTO','ZP_AP3','ZP_AP5','ZP_AP7',
	'ZP_F10','ZP_F15','ZP_F20','ZP_F25','ZP_F30','ZP_OPTA']
zpetypes = ['ZPE_PSF','ZPE_AUTO','ZPE_AP3','ZPE_AP5','ZPE_AP7',
	'ZPE_F10','ZPE_F15','ZPE_F20','ZPE_F25','ZPE_F30','ZPE_OPTA']
ul5types = ['UL5_PSF','UL5_AUTO','UL5_AP3','UL5_AP5','UL5_AP7',
	'UL5_F10','UL5_F15','UL5_F20','UL5_F25','UL5_F30','UL5_OPTA']
#def lcplot(w):

datfile=datfiles[w]
frag=os.path.splitext(datfile)[0].split('-')
tbl=ascii.read(datdirectory+datfile)
fn=tbl['FILE']
comcut=[s for s in list(fn) if 'com' in s and 'cut' in s]
com=[s for s in list(fn) if 'com' in s and 'cut' not in s]
single=[s for s in list(fn) if 'com' not in s and 'cut' not in s]
cut=[s for s in list(fn) if 'com' not in s and 'cut' in s]
comid=[np.where(fn==i) for i in com]
comcutid=[np.where(fn==i) for i in comcut]
cutid=[np.where(fn==i) for i in cut]
singleid=[np.where(fn==i) for i in single]
detyid=np.where(tbl['detect']=='Y')[0]
detnid=np.where(tbl['detect']=='N')[0]
comtbl=tbl[comid]
cuttbl=tbl[cutid]
comcuttbl=tbl[comcutid]
othertbl=tbl[singleid]
detytbl=tbl[detyid]
detntbl=tbl[detnid]

magtype='MAG_APER_5'
magerrtype=magtype[:3]+'ERR'+magtype[3:]
fig, ax = plt.subplots()
ax.set_xlim(58120,58350)
ax.set_ylim(22,14)
ax.set_xlabel('MJD')
ax.set_ylabel('MAG (AB)')
ax.errorbar(comtbl['MJD'],comtbl[magtype],yerr=comtbl[magerrtype].reshape(-1),
                fmt='ro',label=frag[-2]+' : '+frag[-1])
ax.errorbar(cuttbl['MJD'],cuttbl[magtype],yerr=cuttbl[magerrtype].reshape(-1),
                fmt='bo',label=frag[-2]+' : '+frag[-1])
ax.errorbar(comcuttbl['MJD'],comcuttbl[magtype],yerr=comcuttbl[magerrtype].reshape(-1),
                fmt='go',label=frag[-2]+' : '+frag[-1])
ax.errorbar(othertbl['MJD'],othertbl[magtype],yerr=othertbl[magerrtype].reshape(-1),
                fmt='yo',label=frag[-2]+' : '+frag[-1])
ax.errorbar(detntbl['MJD'],detntbl['UL5_F15'],yerr=0.1,lolims=detntbl['ZPE_F15'],fmt='ro',fillstyle='none')

ax.vlines(58142.3625,ymin=12,ymax=23,label='Discovery',linestyle='--')
