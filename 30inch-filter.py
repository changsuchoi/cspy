
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clipping
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.io.ascii as ascii
from scipy.interpolate import UnivariateSpline
from multiprocessing import Process,Pool
import matplotlib
matplotlib.font_manager._rebuild()


seconfigdir ='/data7/cschoi/code/cspy/sex.config/'
seconfig    ='se1.sex'
separam     ='se1.param'
separam_noPSF = 'se1_noPSF.param'
growthparam = 'growth.param'
seconv      ='default.conv'
sennw       ='default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH  = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)
lowmag=13
highmag=19
#filname,filerr='R','Rerr'
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
#magtype=magtypes[0]
#refcat='../../ps1-Tonry-NGC3367.cat'
# source extractor command
psf=True
def matching(intbl, reftbl, inra, indec, refra, refdec, sep=2.0):
    """
    MATCHING TWO CATALOG WITH RA, Dec COORD. WITH python
    INPUT   :   SE catalog, SDSS catalog file name, sepertation [arcsec]
    OUTPUT  :   MATCED CATALOG FILE & TABLE
    """
    incoord     = SkyCoord(inra, indec, unit=(u.deg, u.deg))
    refcoord    = SkyCoord(refra, refdec, unit=(u.deg, u.deg))
    #   INDEX FOR REF.TABLE
    indx, d2d, d3d  = incoord.match_to_catalog_sky(refcoord)
    mreftbl         = reftbl[indx]
    mreftbl['sep']  = d2d
    mergetbl        = intbl
    for col in mreftbl.colnames:
        mergetbl[col]   = mreftbl[col]
    indx_sep        = np.where(mergetbl['sep']*3600.<sep)
    mtbl            = mergetbl[indx_sep]
    #mtbl.write(mergename, format='ascii', overwrite=True)
    return mtbl

def starcut(mtbl,lowmag=lowmag,highmag=highmag,filname=filname,magtype=magtype):
	idx=np.where( (mtbl['SNR_WIN'] >10) &
				(mtbl['FLAGS'] == 0) &
				(mtbl[filname] < highmag) &
				(mtbl[filname] > lowmag) &
				(mtbl[magtype[:3]+'ERR'+magtype[3:]]<0.5)		)
	return mtbl[idx]

magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']

filset=[ ['B','Berr'],['R','Rerr'],['V','Verr'],['I','Ierr'] ]

def zpcal(mtbl1,filname, magtype):
    zp=mtbl1[filname]-mtbl1[magtype]
    #zp3=sigma_clipped_stats(zp, sigma=3, maxiters=10)
    zp2=sigma_clipped_stats(zp, sigma=2, maxiters=10)
    print ('zp ',zp2[0], 'zp err',zp2[2])
    filtered_data=sigma_clip(zp,sigma=2,maxiters=10)
    selected, nonselected= ~filtered_data.mask, filtered_data.mask
    zperrp=np.sqrt( np.sum(mtbl1[filerr][selected]**2 + \
						mtbl1[magtype[:3]+'ERR'+magtype[3:]][selected]**2)\
						/ len(mtbl1) )
    print(magtype, 'zp', '{},'.format(round(zp2[0],3)),
	 	'zperr', '{},'.format(round(zperrp,3)),
		len(mtbl1[selected]),'stars from',len(mtbl1))
    return zp2, selected,zperrp

def zp_plot(mtbl1, zp2, selected, im, magtype='MAG_APER_3', filname=filname, filerr=filerr):
	fn=os.path.splitext(im)[0]
	zp=mtbl1[filname]-mtbl1[magtype]
	xr=np.linspace(np.min(mtbl1[filname]), np.max(mtbl1[filname]), len(zp))
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	zperrp=np.sqrt( np.sum(mtbl1[filerr][selected]**2 + \
						mtbl1[magtype[:3]+'ERR'+magtype[3:]][selected]**2)\
						/ len(mtbl1) )
#	plt.plot(mtbl1[filname],zp,'o',markersize=12)
	plt.errorbar(mtbl1[filname],zp,yerr=zperrp,fmt='bo',markersize=12)
	plt.ylim(zp2[0]-1,zp2[0]+1)
	plt.xlim(np.min(mtbl1[filname]),np.max(mtbl1[filname]))
	#plt.hlines(zp3[0],xmin=12,xmax=20,color='b')
	plt.hlines(zp2[0],xmin=12,xmax=20,color='r')
	sig2=np.ones(len(mtbl1))*zp2[2]
	zp2a=np.ones(len(mtbl1))*zp2[0]
	#plt.fill_between(xr,zp3a+sig3,zp3a-sig3,color='b',alpha=0.5)
	plt.fill_between(xr,zp2a+sig2,zp2a-sig2,color='yellow',alpha=0.5)
	plt.plot(mtbl1[filname][~selected],zp[~selected],'ro')
	plt.plot(mtbl1[filname][selected],zp[selected],'ko')
	plt.text(18,zp2[0]+0.8,magtype+' ZP '+str(round(zp2[0],3))+' err '+str(round(zp2[2],3)))
	plt.title(fn+', '+filname+' '+magtype)
	plt.ylabel('Zeropoint (AB Mag)')
	plt.xlabel(filname +' Reference Mag (AB)')
	plt.savefig(fn+'_'+filname+'_'+magtype+'_zp.png')
	plt.close()



se1='saCalib-MCD30INCH-NGC3367-20160305-041500-R-300_com.se1'
ref='../../ps1-Tonry-NGC3367.cat'
def filter_30inch(se1, ref):
	se1cat=ascii.read(se1)
	refcat=ascii.read(ref)
	mtbl=matching(se1cat, refcat, se1cat['ALPHA_J2000'],se1cat['DELTA_J2000'],refcat['ra'],refcat['dec'])
	mtbl=matching(se1cat, refcat, se1cat['ALPHA_J2000'],se1cat['DELTA_J2000'],refcat['ra'],refcat['dec'])
	filset=[ ['B','Berr'],['R','Rerr'],['V','Verr'],['I','Ierr']]
	filname, filerr=filset[0]
	magtype='MAG_APER_2'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	lowmag=13
	highmag=19
	mtblB=starcut(mtbl, filname='B', magtype=magtype)
	mtblV=starcut(mtbl, filname='V', magtype=magtype)
	mtblR=starcut(mtbl, filname='R', magtype=magtype)
	mtblI=starcut(mtbl, filname='I', magtype=magtype)
	Rzp=zpcal(mtblR,'R','MAG_APER_2')
	zpB=mtbl['B']-mtbl[magtype]
	plt.plot(mtbl['B'], zpB,'b.',label='B')
	zpV=mtbl['V']-mtbl[magtype]
	plt.plot(mtbl['V'], zpV,'g.',label='V')
	zpR=mtbl['R']-mtbl[magtype]
	plt.plot(mtbl['R'], zpR,'r.',label='R')
	zpI=mtbl['I']-mtbl[magtype]
	plt.plot(mtbl['I'], zpI,'y.',label='I')
	plt.xlabel('REF MAG')
	plt.ylabel('SE MAG')
	plt.title(os.path.splitext(se1)[0])
	plt.legend(loc='lower left')
	plt.xlim(13,20)
	plt.ylim(Rzp[0][0]-3,Rzp[0][0]+3)
	plt.text(13.5,27.5,'B'+'\t'+str(round( np.std(mtblB['B']-mtblB[magtype]),2)) )
	plt.text(13.5,27.0,'V'+'\t'+str(round( np.std(mtblV['V']-mtblV[magtype]),2)) )
	plt.text(13.5,26.5,'R'+'\t'+str(round( np.std(mtblR['R']-mtblR[magtype]),2)) )
	plt.text(13.5,26.0,'I'+'\t'+str(round( np.std(mtblI['I']-mtblI[magtype]),2)) )
	plt.savefig(os.path.splitext(se1)[0]+'_chfil.jpg')
	plt.close()
	print('B sigma',np.std(mtblB['B']-mtblB[magtype]))
	print('V sigma',np.std(mtblV['V']-mtblV[magtype]))
	print('R sigma',np.std(mtblR['R']-mtblR[magtype]))
	print('I sigma',np.std(mtblI['I']-mtblI[magtype]))
	return [np.std(mtblB['B']-mtblB[magtype]),np.std(mtblV['V']-mtblV[magtype]),
			np.std(mtblR['R']-mtblR[magtype]),np.std(mtblI['I']-mtblI[magtype])]


# run
selist=glob.glob('saCalib*se1')
selist.sort()
plt.ioff()
stdfil=[]
for n,se1 in enumerate(selist):
	print(n+1,'of',len(selist),'\n', os.path.splitext(se1)[0],'\n')
	stdfil.append(filter_30inch(se1,ref))

# check
fils=['B','V','R,','I']
for n,i in enumerate(stdfil):
	if list(i)[3] != np.min(list(i)):
		#list(i).index(np.min(list(i)))
		print (selist[n],'\t', fils[list(i).index(np.min(list(i)))],
				round(i[0],2),round(i[1],2),round(i[2],2),round(i[3],2)
				)
