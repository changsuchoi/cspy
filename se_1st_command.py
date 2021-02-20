# im='Calib-LOAO-NGC3367-20131023-123229-R-60.fits'
fn=os.path.splitext(im)[0]
skyval, skysig=secom(im,psf=True)
setbl=ascii.read(fn+'.se1')
refcat='../../ps1-Tonry-NGC3367.cat'
reftbl=ascii.read(refcat)
magtypes=['MAG_AUTO', 'MAG_PSF', 'MAG_APER','MAG_APER_1','MAG_APER_2','MAG_APER_3']
mtbl=matching(setbl, pstbl, setbl['ALPHA_J2000'],setbl['DELTA_J2000'],pstbl['ra'],pstbl['dec'])
# magtype='MAG_PSF'
mtbl1=starcut(mtbl,lowmag=14,highmag=19,filname='R',magtype=magtype)

filname,filerr='R','Rerr'

# magerr=magtype[:3]+'ERR'+magtype[3:]
for magtype in magtypes:
magerrtype=mag[:3]+'ERR'+mag[3:]
mtbl1=starcut(mtbl,filname='R',magtype=mag)
zp2,selected, fwhm_img,zperr=zpcal(mtbl1,filname,mag)
zp_plot(mtbl1,zp2,selected,mag,im,filname='R',filerr='Rerr')
fitplot(im,mtbl1,mag,selected)
