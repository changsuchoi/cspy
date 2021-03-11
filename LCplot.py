import astropy.io.ascii as ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.table import Table, vstack

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

cctbl=vstack(comtbl,comcuttbl)
cctbl.sort()
ascii.write(cctbl, 'SN2018kp-LOAO-B.dat.com',format='commented_header',overwrite=True)



plt.errorbar(comtbl['MJD'], comtbl['MAG_APER_8'], yerr=comtbl['MAGERR_APER_8'].reshape(-1),fmt='bo')
plt.errorbar(cuttbl['MJD'], cuttbl['MAG_APER_8'], yerr=cuttbl['MAGERR_APER_8'].reshape(-1),fmt='ro')
plt.errorbar(comcuttbl['MJD'], comcuttbl['MAG_APER_8'], yerr=comcuttbl['MAGERR_APER_8'].reshape(-1),fmt='yo')
plt.errorbar(othertbl['MJD'], othertbl['MAG_APER_8'], yerr=othertbl['MAGERR_APER_8'].reshape(-1),fmt='go')

plt.xlabel(MJD)
plt.ylable('MAG(AB)')
plt.title(datfiles)
plt.xlim(discoverymjd-20,discoverymjd+100)
plt.ylim(22,14)
