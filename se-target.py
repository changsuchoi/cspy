#source extractor on subtracted images

import os
import glob
import astropy.io.fits as fits
import numpy as np
import subprocess
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clipping
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.io.ascii as ascii
from scipy.interpolate import UnivariateSpline
from multiprocessing import Process,Pool
from astropy.time import Time
from astropy.wcs import WCS
import astropy.wcs.utils as wcsutils
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord


seconfigdir ='/data7/cschoi/code/cspy/sex.config/'
seconfig    ='se1.sex'
separam     ='se1.param'
separam_noPSF = 'se1_noPSF.param'
#growthparam = 'growth.param'
seconv      ='default.conv'
sennw       ='default.nnw'
DETECT_MINAREA = str(5)
DETECT_THRESH  = str(3)
DEBLEND_NTHRESH = str(32)
DEBLEND_MINCONT = str(0.005)
# lowmag=13
# highmag=19
# filname,filerr='R','Rerr'
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
#magtype=magtypes[0]

#tra,tdec=[]
def targetfind(tra, tdec, refra, refdec, sep=2.0):
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	targ_coord	= SkyCoord(tra, tdec, unit=(u.deg, u.deg))
	phot_coord	= SkyCoord(refra, refdec, unit=(u.deg, u.deg))
	indx, d2d, d3d	= targ_coord.match_to_catalog_sky(phot_coord)
	return indx.item(), round(d2d.arcsec.item(),2)

# isolated field object, no need to subtraciotn
psf=True
def setarget(im,psf=psf):
	DETECT_MINAREA = str(3)
	DETECT_THRESH  = str(1.5)
	DEBLEND_NTHRESH = str(32)
	DEBLEND_MINCONT = str(0.005)
	PSCALE=fits.getheader(im)['PSCALE']
	fwhm_pix=fits.getheader(im)['FWHM_PIX']
	opt_ap=fits.getheader(im)['OPT_AP']
	aper_list=[3,5,7]
	aper_list_fwhm=[1.0,1.5,2.0,2.5,3.0]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	#aper_input = aper_input[:-1]
	for i in aper_list_fwhm: aper_input += '{},'.format(round(i*fwhm_pix,1))
	aper_input = aper_input+str(opt_ap)
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.set'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	opt5a=' -CHECKIMAGE_TYPE NONE '
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	opt9=' -SEEING_FWHM '+str(round(fwhm_pix*PSCALE))+ ' '
	if psf==True:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 + im
	else:
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 + im
	os.system(secommand)

# for subtracted images
psf=True
# se_sub(im,psf=psf)
def se_sub(subim,psf=psf,dual=False,det=''):
	subfn=os.path.splitext(subim)[0]
	psfname=subfn[6:]+'.psf '
	if 'cut' in subfn : psfname=subfn[6:][:-10]+'.psf'
	PSCALE=fits.getheader(subim)['PSCALE']
	fwhm_pix=fits.getheader(subim)['FWHM_PIX']
	opt_ap=fits.getheader(subim)['OPT_AP']
	DETECT_MINAREA = str(3)
	DETECT_THRESH  = str(1.5)
	DEBLEND_NTHRESH = str(32)
	DEBLEND_MINCONT = str(0.005)
	aper_list=[3,5,7]
	aper_list_fwhm=[1.0,1.5,2.0,2.5,3.0]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	#aper_input = aper_input[:-1]
	for i in aper_list_fwhm: aper_input += '{},'.format(round(i*fwhm_pix,1))
	aper_input = aper_input+str(opt_ap)
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ subfn+'.sub'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+subfn+'_seg.fits'+','+subfn+'_ap.fits'
	opt5a=' -CHECKIMAGE_TYPE  SEGMENTATION'
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+psfname
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	opt9=' -SEEING_FWHM '+str(round(fwhm_pix*PSCALE))+ ' '
	if (psf==True) & (dual==True):
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 +det+' '+ subim
	elif (psf==True) & (dual==False):
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 + subim
	elif (psf==False) & (dual==True):
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 +det+' '+ subim
	else:
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 + subim
	os.system(secommand)


def se_zp(im,ext='.set'):
	from astropy.table import Table
	from astropy import units as u
	from astropy.table import Column
	fn = os.path.splitext(im)[0]
	hdr=fits.getheader(im)
	sefcat=ascii.read(fn+ext)
	#print('Making final catalog',fn+'.dat')
	sefcat['MAG_AUTO']=sefcat['MAG_AUTO']+hdr['ZP_AUTO']
	sefcat['MAGERR_AUTO']=np.sqrt(sefcat['MAGERR_AUTO']**2 + hdr['ZPE_AUTO']**2)
	sefcat['MAG_APER']=sefcat['MAG_APER']+hdr['ZP_AP3']
	sefcat['MAGERR_APER']=np.sqrt(sefcat['MAGERR_APER']**2 + hdr['ZPE_AP3']**2)
	sefcat['MAG_APER_1']=sefcat['MAG_APER_1']+hdr['ZP_AP5']
	sefcat['MAGERR_APER_1']=np.sqrt(sefcat['MAGERR_APER_1']**2 + hdr['ZPE_AP5']**2)
	sefcat['MAG_APER_2']=sefcat['MAG_APER_2']+hdr['ZP_AP7']
	sefcat['MAGERR_APER_2']=np.sqrt(sefcat['MAGERR_APER_2']**2 + hdr['ZPE_AP7']**2)
	sefcat['MAG_APER_3']=sefcat['MAG_APER_3']+hdr['ZP_F10']
	sefcat['MAGERR_APER_3']=np.sqrt(sefcat['MAGERR_APER_3']**2 + hdr['ZPE_F10']**2)
	sefcat['MAG_APER_4']=sefcat['MAG_APER_4']+hdr['ZP_F15']
	sefcat['MAGERR_APER_4']=np.sqrt(sefcat['MAGERR_APER_4']**2 + hdr['ZPE_F15']**2)
	sefcat['MAG_APER_5']=sefcat['MAG_APER_5']+hdr['ZP_F20']
	sefcat['MAGERR_APER_5']=np.sqrt(sefcat['MAGERR_APER_5']**2 + hdr['ZPE_F20']**2)
	sefcat['MAG_APER_6']=sefcat['MAG_APER_6']+hdr['ZP_F25']
	sefcat['MAGERR_APER_6']=np.sqrt(sefcat['MAGERR_APER_6']**2 + hdr['ZPE_F25']**2)
	sefcat['MAG_APER_7']=sefcat['MAG_APER_7']+hdr['ZP_F30']
	sefcat['MAGERR_APER_7']=np.sqrt(sefcat['MAGERR_APER_7']**2 + hdr['ZPE_F30']**2)
	sefcat['MAG_APER_8']=sefcat['MAG_APER_8']+hdr['ZP_OPTA']
	sefcat['MAGERR_APER_8']=np.sqrt(sefcat['MAGERR_APER_8']**2 + hdr['ZPE_OPTA']**2)
	if 'MAG_PSF' in sefcat.colnames:
		sefcat['MAG_PSF']=sefcat['MAG_PSF']+hdr['ZP_PSF']
		sefcat['MAGERR_PSF']=np.sqrt(sefcat['MAGERR_PSF']**2 + hdr['ZPE_PSF']**2)
	sefcat.write(fn+'.sep',format='ascii.commented_header',overwrite=True)
	return sefcat

def se_zp_sub(im,ext='.sub'):
	from astropy.table import Table
	from astropy import units as u
	from astropy.table import Column
	fn = os.path.splitext(im)[0]
	hdr=fits.getheader(im)
	sefcat=ascii.read(fn+ext)
	#print('Making final catalog',fn+'.dat')
	sefcat['MAG_AUTO']=sefcat['MAG_AUTO']+hdr['ZP_AUTO']
	sefcat['MAGERR_AUTO']=np.sqrt(sefcat['MAGERR_AUTO']**2 + hdr['ZPE_AUTO']**2)
	sefcat['MAG_APER']=sefcat['MAG_APER']+hdr['ZP_AP3']
	sefcat['MAGERR_APER']=np.sqrt(sefcat['MAGERR_APER']**2 + hdr['ZPE_AP3']**2)
	sefcat['MAG_APER_1']=sefcat['MAG_APER_1']+hdr['ZP_AP5']
	sefcat['MAGERR_APER_1']=np.sqrt(sefcat['MAGERR_APER_1']**2 + hdr['ZPE_AP5']**2)
	sefcat['MAG_APER_2']=sefcat['MAG_APER_2']+hdr['ZP_AP7']
	sefcat['MAGERR_APER_2']=np.sqrt(sefcat['MAGERR_APER_2']**2 + hdr['ZPE_AP7']**2)
	sefcat['MAG_APER_3']=sefcat['MAG_APER_3']+hdr['ZP_F10']
	sefcat['MAGERR_APER_3']=np.sqrt(sefcat['MAGERR_APER_3']**2 + hdr['ZPE_F10']**2)
	sefcat['MAG_APER_4']=sefcat['MAG_APER_4']+hdr['ZP_F15']
	sefcat['MAGERR_APER_4']=np.sqrt(sefcat['MAGERR_APER_4']**2 + hdr['ZPE_F15']**2)
	sefcat['MAG_APER_5']=sefcat['MAG_APER_5']+hdr['ZP_F20']
	sefcat['MAGERR_APER_5']=np.sqrt(sefcat['MAGERR_APER_5']**2 + hdr['ZPE_F20']**2)
	sefcat['MAG_APER_6']=sefcat['MAG_APER_6']+hdr['ZP_F25']
	sefcat['MAGERR_APER_6']=np.sqrt(sefcat['MAGERR_APER_6']**2 + hdr['ZPE_F25']**2)
	sefcat['MAG_APER_7']=sefcat['MAG_APER_7']+hdr['ZP_F30']
	sefcat['MAGERR_APER_7']=np.sqrt(sefcat['MAGERR_APER_7']**2 + hdr['ZPE_F30']**2)
	sefcat['MAG_APER_8']=sefcat['MAG_APER_8']+hdr['ZP_OPTA']
	sefcat['MAGERR_APER_8']=np.sqrt(sefcat['MAGERR_APER_8']**2 + hdr['ZPE_OPTA']**2)
	if 'MAG_PSF' in sefcat.colnames:
		sefcat['MAG_PSF']=sefcat['MAG_PSF']+hdr['ZP_PSF']
		sefcat['MAGERR_PSF']=np.sqrt(sefcat['MAGERR_PSF']**2 + hdr['ZPE_PSF']**2)
	sefcat.write(fn+'.sep',format='ascii.commented_header',overwrite=True)
	return sefcat

tblcols=['NUMBER','X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000',
'MAG_AUTO','MAGERR_AUTO','FLUX_AUTO','FLUXERR_AUTO',
'MAG_PSF','MAGERR_PSF','FLUX_PSF','FLUXERR_PSF',
'MAG_APER','MAG_APER_1','MAG_APER_2',
'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7','MAG_APER_8',
'MAGERR_APER','MAGERR_APER_1','MAGERR_APER_2',
'MAGERR_APER_3','MAGERR_APER_4','MAGERR_APER_5','MAGERR_APER_6','MAGERR_APER_7','MAGERR_APER_8',
'FLUX_APER','FLUX_APER_1','FLUX_APER_2',
'FLUX_APER_3','FLUX_APER_4','FLUX_APER_5','FLUX_APER_6','FLUX_APER_7','FLUX_APER_8',
'FLUXERR_APER','FLUXERR_APER_1','FLUXERR_APER_2',
'FLUXERR_APER_3','FLUXERR_APER_4','FLUXERR_APER_5','FLUXERR_APER_6','FLUXERR_APER_7','FLUXERR_APER_8',
'SNR_WIN','BACKGROUND','THRESHOLD','FLUX_MAX',
'MU_THRESHOLD','MU_MAX','FLAGS',
'FWHM_IMAGE','FWHM_WORLD','ELONGATION','ELLIPTICITY','CLASS_STAR']
hdrkey=['PSCALE', 'FWHM_PIX','SKYVAL',  'SKYSIG', 'OPT_AP',  'MATCHNUM','REFCAT',  'LOWMAG',
'HIGHMAG', 'NUM_PSF' ,'ZP_PSF',  'ZPE_PSF','UL5_PSF', 'NUM_AUTO','ZP_AUTO', 'ZPE_AUTO',
'UL5_AUTO','NUM_AP3' ,'ZP_AP3',  'ZPE_AP3','UL5_AP3', 'NUM_AP5' ,'ZP_AP5',  'ZPE_AP5',
'UL5_AP5', 'NUM_AP7' ,'ZP_AP7',  'ZPE_AP7','UL5_AP7', 'NUM_F10' ,'ZP_F10',  'ZPE_F10',
'UL5_F10', 'NUM_F15' ,'ZP_F15',  'ZPE_F15','UL5_F15', 'NUM_F20' ,'ZP_F20',  'ZPE_F20',
'UL5_F20', 'NUM_F25' ,'ZP_F25',  'ZPE_F25','UL5_F25', 'NUM_F30' ,'ZP_F30',  'ZPE_F30',
'UL5_F30', 'NUM_OPTA','ZP_OPTA', 'ZPE_OPTA','UL5_OPTA']


# tra,tdec=161.6715913,13.7165212
def field_target_phot(im, tra, tdec, sep=3.0):
	hdr=fits.getheader(im)
	setarget(im)
	tbl=se_zp(im)
	t = Time(hdr['DATE-OBS'], format='isot', scale='utc').mjd
	trow, tsep = targetfind(tra, tdec,
				refra=tbl['ALPHA_J2000'], refdec=tbl['DELTA_J2000'], sep=sep)
	if tsep <= 2.0:
		detect='Y'
		body1=im+'\t'+str(round(t,7))+'\t'+im.split('-')[1]+'\t'+ im.split('-')[5]+'\t'
		body2=''
		for j in tblcols: body2+=str(tbl[trow][j])+'\t'
		body3=''
		for h in hdrkey: body3+=str(hdr[h])+'\t'
		body=body1+body2+body3+str(tsep)+'\t'+detect+'\n'
	else :
		detect='N'
		print(tsep, 'Not detected')
		body1=im+'\t'+str(round(t,7))+'\t'+im.split('-')[1]+'\t'+ im.split('-')[5]+'\t'
		body2=''
		for j in tblcols: body2+=str(-99)+'\t'
		body3=''
		for h in hdrkey: body3+=str(hdr[h])+'\t'
		body=body1+body2+body3+str(tsep)+'\t'+detect+'\n'
	return body,detect

def sub_target_phot(im, tra, tdec,sep=3.0):
	# im='hdreg_Calib-MAO_FLI-NGC3367-20170124-211937-R-60_4_com.fits'
	#hdr=fits.getheader(im[6:])
	hdr=fits.getheader(im)
	t = Time(hdr['DATE-OBS'], format='isot', scale='utc').mjd
	se_sub(im)
	tbl=se_zp_sub(im)
	trow, tsep = targetfind(tra, tdec,
				refra=tbl['ALPHA_J2000'], refdec=tbl['DELTA_J2000'], sep=sep)
	if tsep <= 2.0:
		detect='Y'
		body1=im+'\t'+str(round(t,7))+'\t'+im.split('-')[1]+'\t'+ im.split('-')[5]+'\t'
		body2=''
		for j in tblcols: body2+=str(tbl[trow][j])+'\t'
		body3=''
		for h in hdrkey: body3+=str(hdr[h])+'\t'
		body=body1+body2+body3+str(tsep)+'\t'+detect+'\n'
	else :
		detect='N'
		print(tsep, 'Not detected')
		body1=im+'\t'+str(round(t,7))+'\t'+im.split('-')[1]+'\t'+ im.split('-')[5]+'\t'
		body2=''
		for j in tblcols: body2+=str(-99)+'\t'
		body3=''
		for h in hdrkey: body3+=str(hdr[h])+'\t'
		body=body1+body2+body3+str(tsep)+'\t'+detect+'\n'
	return body,detect

pos=(161.637792 +13.741869)
tra,tdec=161.637792, 13.741869
pos=(tra,tdec)
sz=(5,5)
def trimstamp(inim, positions=pos, sizes=sz):
	'''
	positions=(ra,dec) in deg unit
	sizes=(px,py) in pixel unit
	'''
	hdu = fits.open(inim)[0]
	hdr = hdu.header
	w = WCS(hdu.header)
	outim=os.path.splitext(inim)[0]+'_'+'{0:02d}'.format(sizes[0])+'min_stamp.fits'
	print('Input center position(deg)',positions)
	xpscale,ypscale=wcsutils.proj_plane_pixel_scales(w)*3600
	pixscale=(xpscale+ypscale)/2.
	size=u.Quantity(sizes,u.arcmin)
	print('sizes',size)
	positions=SkyCoord(positions[0]*u.deg, positions[1]*u.deg)
	cutout = Cutout2D(hdu.data, position=positions, size=size, wcs=w,
					mode='trim',fill_value=1.0e-30)
	hdu.data = cutout.data
	hdu.header.update(cutout.wcs.to_header())
	hdu.writeto(outim, overwrite=True)
	return 'Done'

def fitplot_target_stamp(inim,tra,tdec,sizes):
	im=os.path.splitext(inim)[0]+'_'+'{0:02d}'.format(sizes[0])+'min_stamp.fits'
	import matplotlib.pyplot as plt
	from astropy.wcs import WCS
	from astropy.io import fits
	from astropy.visualization import MinMaxInterval,ZScaleInterval,PercentileInterval
	from astropy.visualization import SqrtStretch,LinearStretch
	from astropy.visualization import ImageNormalize
	imdata,imhdr=fits.getdata(im,header=True)
	norm = ImageNormalize(imdata, interval=ZScaleInterval(), stretch=LinearStretch() )
	wcs=WCS(imhdr)
	fig,ax=plt.subplots()
	ax=plt.subplot(projection=wcs)
	ax.set_xlabel('RA')
	ax.set_ylabel('DEC')
	#ax.invert_xaxis()
	ax.set_title(im)
	ax.scatter(tra,tdec, transform=ax.get_transform('fk5'),
				s=100, edgecolor='yellow',facecolor='none')
	#ax.scatter(mtbl1[~selected]['ALPHA_J2000'],mtbl1[~selected]['DELTA_J2000'],
	#	transform=ax.get_transform('fk5'),s=20, edgecolor='red',facecolor='none')
	img=ax.imshow(imdata,cmap='gray',norm=norm,origin='lower')
	ax.invert_yaxis()
	fig.colorbar(img)
	plt.savefig(os.path.splitext(inim)[0]+'_'+'{0:02d}'.format(sizes[0])+'min_stamp.png')
	plt.close()


# for field target: use below
'''
targetname='test'
tra,tdec=161.6715913,13.7165212
pos=(tra,tdec)
sizes=(5,5)
sz=sizes
thead=''
for j in tblcols: thead+=j+'\t'
for h in hdrkey: thead+=h+'\t'
resulthead= '#FILE'+'\t'+'MJD'+ '\t'+'OBSERVATORY'+'\t'+'FILTER'+'\t'+thead+'sep'+'\t'+'detect'+'\n'
tra,tdec=161.6715913,13.7165212
f=open(targetname+'-'+os.getcwd().split('/')[-2]+'-'+os.getcwd().split('/')[-1]+'.dat','w')
f.write(resulthead)
for im in imlist:
	body,detect=field_target_phot(im, tra, tdec)
	print( im,detect,'\n',body)
	f.write(body)
	#trimstamp(im, positions=pos, sizes=sz)
	#fitplot_target_stamp(im,tra,tdec,sizes)
f.close()
'''
# for sub target: use below
'''
targetname='SN2018kp'
tra,tdec=161.637792, 13.741869
pos=(tra,tdec)
sizes=(5,5)
sz=sizes
thead=''
for j in tblcols: thead+=j+'\t'
for h in hdrkey: thead+=h+'\t'
resulthead= '#FILE'+'\t'+'MJD'+'\t'+'OBSERVATORY'+'\t'+'FILTER'+'\t'+thead+'\t'+'sep'+'\t'+'detect'+'\n'
f=open(targetname+'-'+os.getcwd().split('/')[-2]+'-'+os.getcwd().split('/')[-1]+'.dat','w')
f.write(resulthead)
for im in hdlist:
	body,detect=sub_target_phot(im, tra, tdec)
	print( im, detect,'\n',body)
	f.write(body)
	trimstamp(im, positions=pos, sizes=sz)
	fitplot_target_stamp(im,tra,tdec,sizes)
f.close()
'''











'''
('NUMBER','X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000',
'MAG_AUTO','MAGERR_AUTO','FLUX_AUTO','FLUXERR_AUTO',
'MAG_PSF','MAGERR_PSF','FLUX_PSF','FLUXERR_PSF',
'MAG_APER','MAG_APER_1','MAG_APER_2',
'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7','MAG_APER_8',
'MAGERR_APER','MAGERR_APER_1','MAGERR_APER_2',
'MAGERR_APER_3','MAGERR_APER_4','MAGERR_APER_5','MAGERR_APER_6','MAGERR_APER_7','MAGERR_APER_8',
'FLUX_APER','FLUX_APER_1','FLUX_APER_2',
'FLUX_APER_3','FLUX_APER_4','FLUX_APER_5','FLUX_APER_6','FLUX_APER_7','FLUX_APER_8',
'FLUXERR_APER','FLUXERR_APER_1','FLUXERR_APER_2',
'FLUXERR_APER_3','FLUXERR_APER_4','FLUXERR_APER_5','FLUXERR_APER_6','FLUXERR_APER_7','FLUXERR_APER_8',
'SNR_WIN','BACKGROUND','THRESHOLD','FLUX_MAX',
'MU_THRESHOLD','MU_MAX','FLAGS',
'FWHM_IMAGE','FWHM_WORLD','ELONGATION','ELLIPTICITY','CLASS_STAR')
'''

'''
PSCALE  =     0.91990318614984
FWHM_PIX=                4.011 / FWHM PIXEL VALUE
SKYVAL  =              122.625 / sky median value form sextractor
SKYSIG  =              4.94231 / sky sigma value form sextractor
OPT_AP  =                 5.71 / Optimal aperture size
MATCHNUM=                   19 / Number of matched stars
REFCAT  = 'ps1-Tonry-NGC3367.cat' / Referenc Catalog used
LOWMAG  =                   13 / Low MAG CUT for ZP Calculation
HIGHMAG =                   19 / HIGH MAG CUT for ZP Calculation
NUM_PSF =                   11 / number of stars used for zp calculation
ZP_PSF  =               24.141 / MAG_PSF ZERO POINT(AB)
ZPE_PSF =                 0.03 / MAG_PSF ZERO POINT ERROR
UL5_PSF =               18.781 / MAG_PSF 5 sigma upperlimit
NUM_AUTO=                    9 / number of stars used for zp calculation
ZP_AUTO =                24.15 / MAG_AUTO ZERO POINT(AB)
ZPE_AUTO=                0.045 / MAG_AUTO ZERO POINT ERROR
UL5_AUTO=                18.71 / MAG_AUTO 5 sigma upperlimit
NUM_AP3 =                   11 / number of stars used for zp calculation
ZP_AP3  =               22.927 / 3 arcsec aperture ZERO POINT(AB)
ZPE_AP3 =                0.039 / 3 arcsec aperture ZERO POINT ERROR
UL5_AP3 =               18.383 / 3 arcsec aperture 5 sigma upperlimit
NUM_AP5 =                   10 / number of stars used for zp calculation
ZP_AP5  =               23.633 / 5 arcsec aperture ZERO POINT(AB)
ZPE_AP5 =                0.032 / 5 arcsec aperture ZERO POINT ERROR
UL5_AP5 =               19.089 / 5 arcsec aperture 5 sigma upperlimit
NUM_AP7 =                   10 / number of stars used for zp calculation
ZP_AP7  =               23.922 / 7 arcsec aperture ZERO POINT(AB)
ZPE_AP7 =                0.032 / 7 arcsec aperture ZERO POINT ERROR
UL5_AP7 =               19.378 / 7 arcsec aperture 5 sigma upperlimit
NUM_F10 =                   11 / number of stars used for zp calculation
ZP_F10  =               23.238 / 1.0 FWHM APERTURE ZERO POINT(AB)
ZPE_F10 =                0.035 / 1.0 FWHM APERTURE ZERO POINT ERROR
UL5_F10 =               18.694 / 1.0 FWHM APERTURE 5 sigma upperlimit
NUM_F15 =                   10 / number of stars used for zp calculation
ZP_F15  =               23.746 / 1.5 FWHM APERTURE ZERO POINT(AB)
ZPE_F15 =                0.032 / 1.5 FWHM APERTURE ZERO POINT ERROR
UL5_F15 =               19.203 / 1.5 FWHM APERTURE 5 sigma upperlimit
NUM_F20 =                   10 / number of stars used for zp calculation
ZP_F20  =               23.955 / 2.0 FWHM APERTURE ZERO POINT(AB)
ZPE_F20 =                0.033 / 2.0 FWHM APERTURE ZERO POINT ERROR
UL5_F20 =               19.411 / 2.0 FWHM APERTURE 5 sigma upperlimit
NUM_F25 =                    9 / number of stars used for zp calculation
ZP_F25  =               24.086 / 2.5 FWHM APERTURE ZERO POINT(AB)
ZPE_F25 =                0.033 / 2.5 FWHM APERTURE ZERO POINT ERROR
UL5_F25 =               19.542 / 2.5 FWHM APERTURE 5 sigma upperlimit
NUM_F30 =                   10 / number of stars used for zp calculation
ZP_F30  =               24.125 / 3.0 FWHM APERTURE ZERO POINT(AB)
ZPE_F30 =                0.041 / 3.0 FWHM APERTURE ZERO POINT ERROR
UL5_F30 =               19.581 / 3.0 FWHM APERTURE 5 sigma upperlimit
NUM_OPTA=                   10 / number of stars used for zp calculation
ZP_OPTA =               23.692 / Optimal APERTURE ZERO POINT(AB)
ZPE_OPTA=                0.032 / Optimal APERTURE ZERO POINT ERROR
UL5_OPTA=               19.149 / Optimal APERTURE 5 sigma upperlimit
'''
