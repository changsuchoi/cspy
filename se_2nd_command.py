# im='Calib-LOAO-NGC3367-20131023-123229-R-60.fits'

filname,filerr='R','Rerr'

def se1st(im):
	fn=os.path.splitext(im)[0]
	# se 1st
	skyval, skysig=secom(im,psf=True)
	puthdr(im, 'SKYVAL', skyval,
		hdrcomment='sky median value form sextractor')
	puthdr(im, 'SKYSIG', skysig,
		hdrcomment='sky sigma value form sextractor')
	setbl=ascii.read(fn+'.se1')
	refcat='../../ps1-Tonry-NGC3367.cat'
	reftbl=ascii.read(refcat)
	magtypes=['MAG_AUTO', 'MAG_PSF', 'MAG_APER','MAG_APER_1','MAG_APER_2','MAG_APER_3']
	magname =[]
	mtbl=matching(setbl, pstbl, setbl['ALPHA_J2000'],setbl['DELTA_J2000'],pstbl['ra'],pstbl['dec'])
	# magtype='MAG_PSF'
	mtbl1=starcut(mtbl,lowmag=14,highmag=19,filname=filname,magtype='MAG_AUTO')
	fwhm_img1=fwhm_img(im,mtbl1)

	#MAG_AUTO
	# magerr=magtype[:3]+'ERR'+magtype[3:]
	#for magtype in magtypes:
	magtype='MAG_AUTO'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=magtype)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	puthdr(im, 'ZP_AUTO',round(zp2[0],3), hdrcomment='MAG_AUTO'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AUTO',round(zperr,3), hdrcomment='MAG_AUTO'+' '+'ZERO POINT ERROR')
	sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	puthdr(im, 'UL5_AUTO', sig5ul,
		hdrcomment='MAG_AUTO'+' '+'5 sigma upper limit')

	#MAG_PSF
	magtype='MAG_PSF'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=mag)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	puthdr(im, 'ZP_PSF',round(zp2[0],3), hdrcomment='MAG_PSF'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_PSF',round(zperr,3), hdrcomment='MAG_PSF'+' '+'ZERO POINT ERROR')
	sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	puthdr(im, 'UL5_PSF', sig5ul, hdrcomment='MAG_PSF'+' '+'5 sigma upperlimit')

	#MAG_APER 3"
	magtype='MAG_APER'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=mag)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	puthdr(im, 'ZP_AP3', round(zp2[0],3), hdrcomment='3 arcsec aperture'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AP3', round(zperr,3), hdrcomment='3 arcsec aperture'+' '+'ZERO POINT ERROR')
	sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_AP3', ul5, hdrcomment='AP3'+' '+'5 sigma upperlimit')

	#MAG_APER_1 5"
	magtype='MAG_APER_1'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=mag)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	puthdr(im, 'ZP_AP5', round(zp2[0],3), hdrcomment='5 arcsec aperture'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AP5',round(zperr,3), hdrcomment='5 arcsec aperture'+' '+'ZERO POINT ERROR')
	sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_AP5', ul5,	hdrcomment='AP5'+' '+'5 sigma upperlimit')

	#MAG_APER_2 7"
	magtype='MAG_APER_2'
	magerrtype=magtype[:3]+'ERR'+magtype[3:]
	mtbl1=starcut(mtbl,filname=filname,magtype=mag)
	zp2,selected, zperr=zpcal(mtbl1,filname,magtype)
	zp_plot(mtbl1,zp2,selected,magtype,im,filname=filname,filerr=filerr)
	fitplot(im,mtbl1,magtype,selected)
	puthdr(im, 'ZP_AP7', round(zp2[0],3), hdrcomment='7 arcsec aperture'+' '+'ZERO POINT(AB)' )
	puthdr(im, 'ZPE_AP7',round(zperr,3), hdrcomment='7 arcsec aperture'+' '+'ZERO POINT ERROR')
	sig5ul=	UL_5sig_err(im,setbl,mtbl,mtbl1,magtype,zp2)
	ul5=limitmag(5, zp2[0], 3, skysig)
	puthdr(im, 'UL5_AP5', ul5,	hdrcomment='AP5'+' '+'5 sigma upperlimit')



imlist=glob.glob('Calib*.fits')+\
		glob.glob('saCalib*.fits')+\
		glob.glob('Calib*com.fits')

imlist.sort()
for nn,im in enumerate(imlist) :
	print('=' *60,'\n')
	print(nn+1, 'of ',len(imlist),im)
	se1st(im)
