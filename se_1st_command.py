# im='Calib-LOAO-NGC3367-20131023-123229-R-60.fits'
'''
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
lowmag=14
highmag=19
filname,filerr='R','Rerr'
magtypes=['MAG_AUTO', 'MAG_PSF',
		'MAG_APER','MAG_APER_1','MAG_APER_2',
		'MAG_APER_3','MAG_APER_4','MAG_APER_5','MAG_APER_6','MAG_APER_7',
		'MAG_APER_8']
'''



magtype_dict={'MAG_PSF'    : ['NUM_PSF',  'ZP_PSF',  'ZPE_PSF',  'UL5_PSF'  ,'MAG_PSF'],
			'MAG_AUTO' : ['NUM_AUTO', 'ZP_AUTO', 'ZPE_AUTO', 'UL5_AUTO' ,'MAG_AUTO'],
			'MAG_APER'   : ['NUM_AP3',  'ZP_AP3',  'ZPE_AP3',  'UL5_AP3'  ,'3 arcsec aperture'],
			'MAG_APER_1' : ['NUM_AP5',  'ZP_AP5',  'ZPE_AP5',  'UL5_AP5'  ,'5 arcsec aperture'],
			'MAG_APER_2' : ['NUM_AP7',  'ZP_AP7',  'ZPE_AP7',  'UL5_AP7'  ,'7 arcsec aperture'],
			'MAG_APER_3' : ['NUM_F10',  'ZP_F10',  'ZPE_F10',  'UL5_F10'  ,'1.0 FWHM APERTURE'],
			'MAG_APER_4' : ['NUM_F15',  'ZP_F15',  'ZPE_F15',  'UL5_F15'  ,'1.5 FWHM APERTURE'],
			'MAG_APER_5' : ['NUM_F20',  'ZP_F20',  'ZPE_F20',  'UL5_F20'  ,'2.0 FWHM APERTURE'],
			'MAG_APER_6' : ['NUM_F25',  'ZP_F25',  'ZPE_F25',  'UL5_F25'  ,'2.5 FWHM APERTURE'],
			'MAG_APER_7' : ['NUM_F30',  'ZP_F30',  'ZPE_F30',  'UL5_F30'  ,'3.0 FWHM APERTURE'],
			'MAG_APER_8' : ['NUM_OPTA', 'ZP_OPTA', 'ZPE_OPTA', 'UL5_OPTA' ,'Optimal APERTURE'],
			}


def magtype_zp(im, setbl, mtbl, magtype, filname=filname, filerr=filerr, ploton=True, skysig=''):
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	kk=magtype_dict[magtype]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	starnum=len(mtbl1[selected])
	puthdr(im,kk[0],starnum,hdrcomment='number of stars used for zp calculation')
	print(kk[0],starnum)
	puthdr(im, kk[1], round(zp2[0],3), hdrcomment=kk[4]+' '+'ZERO POINT(AB)' )
	print(kk[1],round(zp2[0],3))
	puthdr(im, kk[2], round(zperr,3), hdrcomment=kk[4]+' '+'ZERO POINT ERROR')
	print(kk[2], round(zperr,3))
	if (magtype == 'MAG_AUTO') or (magtype=='MAG_PSF'):
		ul5 = UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	else: ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, kk[3], ul5, hdrcomment=kk[4]+' '+'5 sigma upperlimit')
	print(kk[3], ul5)
	if ploton==True:
		zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
		fitplot(im,mtbl1,magtype,selected)
	else : print('plot off')


# refcat='../../ps1-Tonry-NGC3367.cat'
psf=True
def se1st(im,psf=psf):
	fn=os.path.splitext(im)[0]
	# se 1st
	skyval, skysig=secom(im,psf=psf)
	setbl=ascii.read(fn+'.se1')
	reftbl=ascii.read(refcat)
	mtbl=matching(setbl, reftbl, setbl['ALPHA_J2000'],setbl['DELTA_J2000'],reftbl['ra'],reftbl['dec'])
	puthdr(im, 'MATCHNUM', len(mtbl), hdrcomment='Number of matched stars')
	if len(mtbl)<5 : return None
	# magtype='MAG_PSF'
	mtbl1=starcut(mtbl,lowmag=13,highmag=19,filname=filname,magtype='MAG_AUTO')
	if len(mtbl1)<3 : return None
	#fwhm_img1=fwhm_img(im,mtbl1)
	puthdr(im,'REFCAT',refcat.split('/')[-1],hdrcomment='Referenc Catalog used')
	puthdr(im,'LOWMAG',lowmag,hdrcomment='Low MAG CUT for ZP Calculation')
	puthdr(im,'HIGHMAG',highmag,hdrcomment='HIGH MAG CUT for ZP Calculation')
	for magtype in magtype_dict :
		magtype_zp(im, setbl, mtbl, magtype,
					filname=filname, filerr=filerr,
					ploton= True, skysig=skysig)
	return 'Done'
'''
	#MAG_PSF
	if psf==True:
		magtype='MAG_PSF'
		magerrtype=magtype[:3]+'ERR'+magtype[3:]
		mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
		zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
		zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
		fitplot(im,mtbl1,magtype,selected)
		starnum=len(mtbl1[selected])
		puthdr(im,'NUM_PSF',starnum,hdrcomment='Number of stars used for zp calculation')
		puthdr(im, 'ZP_PSF',round(zp2[0],3), hdrcomment='MAG_PSF'+' '+'ZERO POINT(AB)' )
		puthdr(im, 'ZPE_PSF',round(zperr,3), hdrcomment='MAG_PSF'+' '+'ZERO POINT ERROR')
		sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
		puthdr(im, 'UL5_PSF', sig5ul, hdrcomment='MAG_PSF'+' '+'5 sigma upperlimit')
	else: pass
	#MAG_AUTO
	# magerr=magtype[:3]+'ERR'+magtype[3:]
	#for magtype in magtypes:
	magtype='MAG_AUTO'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_AUTO',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_AUTO',round(zp2[0],3), hdrcomment='MAG_AUTO'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AUTO',round(zperr,3), hdrcomment='MAG_AUTO'+' '+'ZERO POINT ERROR')
	sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	puthdr(im, 'UL5_AUTO', sig5ul,
		hdrcomment='MAG_AUTO'+' '+'5 sigma upper limit')
	#MAG_APER 3"
	magtype='MAG_APER'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_AP3',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_AP3', round(zp2[0],3), hdrcomment='3 arcsec aperture'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AP3', round(zperr,3), hdrcomment='3 arcsec aperture'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_AP3', ul5, hdrcomment='AP3'+' '+'5 sigma upperlimit')
	#MAG_APER_1 5"
	magtype='MAG_APER_1'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_AP5',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_AP5', round(zp2[0],3), hdrcomment='5 arcsec aperture'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AP5',round(zperr,3), hdrcomment='5 arcsec aperture'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_AP5', ul5,	hdrcomment='AP5'+' '+'5 sigma upperlimit')
	#MAG_APER_2 7"
	magtype='MAG_APER_2'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_AP7',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_AP7', round(zp2[0],3), hdrcomment='7 arcsec aperture'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AP7',round(zperr,3), hdrcomment='7 arcsec aperture'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_AP7', ul5,	hdrcomment='AP7'+' '+'5 sigma upperlimit')
	#MAG_APER_3, 1.0 FWHM aperture size
	magtype='MAG_APER_3'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_F10',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_F10', round(zp2[0],3), hdrcomment='1.0 FWHM APERTURE'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_F10',round(zperr,3), hdrcomment='1.0 FWHM APERTURE'+' '+'ZERO POINT ERROR')
	sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	# ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_F10', ul5,	hdrcomment='1.0 FWHM APERTURE'+' '+'5 sigma upperlimit')
	#MAG_APER_4, 1.5 FWHM aperture size
	magtype='MAG_APER_4'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_F15',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_F15', round(zp2[0],3), hdrcomment='1.5 FWHM APERTURE'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_F15',round(zperr,3), hdrcomment='1.5 FWHM APERTURE'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_F15', ul5,	hdrcomment='1.5 FWHM APERTURE'+' '+'5 sigma upperlimit')
	#MAG_APER_5, 2.0 FWHM aperture size
	magtype='MAG_APER_5'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_F20',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_F20', round(zp2[0],3), hdrcomment='2.0 FWHM APERTURE'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_F20',round(zperr,3), hdrcomment='2.0 FWHM APERTURE'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_F20', ul5,	hdrcomment='2.0 FWHM APERTURE'+' '+'5 sigma upperlimit')
	#MAG_APER_6, 2.5 FWHM aperture size
	magtype='MAG_APER_6'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_F25',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_F25', round(zp2[0],3), hdrcomment='2.5 FWHM APERTURE'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_F25',round(zperr,3), hdrcomment='2.5 FWHM APERTURE'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_F25', ul5,	hdrcomment='2.5 FWHM APERTURE'+' '+'5 sigma upperlimit')
	#MAG_APER_7, 3.0 FWHM aperture size
	magtype='MAG_APER_7'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_F30',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_F30', round(zp2[0],3), hdrcomment='3.0 FWHM APERTURE'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_F30',round(zperr,3), hdrcomment='3.0 FWHM APERTURE'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_F30', ul5,	hdrcomment='3.0 FWHM APERTURE'+' '+'5 sigma upperlimit')
	# MAG_APER_8, optimal aperture size
	magtype='MAG_APER_8'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	starnum=len(mtbl1[selected])
	puthdr(im,'NUM_OPTA',starnum,hdrcomment='number of stars used for zp calculation')
	puthdr(im, 'ZP_OPTA', round(zp2[0],3), hdrcomment='Optimal APERTURE'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_OPTA',round(zperr,3), hdrcomment='Optimal APERTURE'+' '+'ZERO POINT ERROR')
	# sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_OPTA', ul5,	hdrcomment='Optimal APERTURE'+' '+'5 sigma upperlimit')
	os.system('rm '+fn+'.se1')
	return 'Done'
'''

def secom_final(im,psf=psf):
	DETECT_THRESH='1.5'
	hdr=fits.getheader(im)
	PSCALE=hdr['PSCALE']
	skyval, skysig,fwhm,opt_ap=hdr['SKYVAL'],hdr['SKYSIG'],hdr['FWHM_PIX'],hdr['OPT_AP']
	aper_list,aper_list2=[3,5,7],[1.0*fwhm, 1.5*fwhm, 2.0*fwhm, 2.5*fwhm, 3.0*fwhm, opt_ap]
	aper_input = ''
	for i in aper_list: aper_input += '{},'.format(round(i/PSCALE,1))
	for i in aper_list2: aper_input += '{},'.format(round(i,1))
	aper_input = aper_input[:-1]
	fn = os.path.splitext(im)[0]
	opt1= seconfigdir+seconfig+' -CATALOG_TYPE ASCII_HEAD -CATALOG_NAME '+ fn+'.sef'
	opt2a=' -PARAMETERS_NAME '+seconfigdir+separam
	opt2b= ' -PARAMETERS_NAME '+seconfigdir+separam_noPSF
	opt2=' -FILTER_NAME '+seconfigdir+seconv +' -STARNNW_NAME '+seconfigdir+sennw
	opt3=' -DETECT_MINAREA '+ DETECT_MINAREA + ' -DETECT_THRESH '+DETECT_THRESH
	opt4=' -DEBLEND_NTHRESH '+ DEBLEND_NTHRESH +' -DEBLEND_MINCONT '+ DEBLEND_MINCONT
	opt5=' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES ' +\
	 		' -CHECKIMAGE_NAME '+fn+'_seg.fits'+','+fn+'_ap.fits'
	# opt5=' -CHECKIMAGE_TYPE NONE '
	opt6=' -PHOT_APERTURES '+aper_input+' '
	opt7=' -PSF_NAME '+fn+'.psf '
	opt8=' -PIXEL_SCALE '+str(PSCALE)+' '
	opt9=' -SEEING_FWHM '+str(round(PSCALE*fwhm,3))+' '
	if psf==True:
		secommand= 'sex -c '+opt1+opt2+opt2a+opt3+opt4+opt5+opt6+opt7+opt8+opt9 +im
	else:
		secommand= 'sex -c '+opt1+opt2+opt2b+opt3+opt4+opt5+opt6+opt8+opt9 +im
	print(secommand)
	#sexout = subprocess.getoutput(secommand)
	#line = [s for s in sexout.split('\n') if 'RMS' in s]
	#skymed, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	os.system(secommand)


'''
imlist=glob.glob('Calib*.fits')+\
		glob.glob('saCalib*.fits')+\
		glob.glob('Calib*com.fits')

imlist.sort()
badlist=[]
for nn,im in enumerate(imlist) :
	print('=' *60,'\n')
	print(nn+1, 'of ',len(imlist),im)
	if fits.getheader(im).get('UL5_OPTA',default=True):
		try : se1st(im)
		except:
			print(im,'has a problem, add it to badlist')
			badlist.append(im)
	else: pass
print(badlist)
'''
