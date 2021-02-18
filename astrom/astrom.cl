
wcsreset fr_fdobj.M82.20140228.0032.fits physical
wcsreset fr_fdobj.M82.20140228.0032.fits world

hedit fr_fdobj.M82.20140228.0032.fits WAT0_001 'system=image' ver-
hedit fr_fdobj.M82.20140228.0032.fits WAT1_001 'wtype=tan axtype=ra' ver-
hedit fr_fdobj.M82.20140228.0032.fits WAT2_001 'wtype=tan axtype=dec' ver-

hedit fr_fdobj.M82.20140228.0032.fits RADECSYS 'FK5'   add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits EQUINOX 2000. add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CTYPE1 'RA---TAN' add+ ver- 
hedit fr_fdobj.M82.20140228.0032.fits CTYPE2 'DEC--TAN' add+ ver-

hedit fr_fdobj.M82.20140228.0032.fits CRVAL1 148.89732 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CRVAL2 69.64831 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CRPIX1 911 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CRPIX2 1149 add+ ver-

hedit fr_fdobj.M82.20140228.0032.fits CD1_1 2.2E-4 add+ ver- 
hedit fr_fdobj.M82.20140228.0032.fits CD1_2 1.9E-5 add+ ver- 
hedit fr_fdobj.M82.20140228.0032.fits CD2_1 1.9E-5 add+ ver-
hedit fr_fdobj.M82.20140228.0032.fits CD2_2 -2.2E-4 add+ ver-

!sex fr_fdobj.M82.20140228.0032.fits -c astrom.sex -CATALOG_NAME=fr_fdobj.M82.20140228.0032.cat 
!scamp fr_fdobj.M82.20140228.0032.cat -c astrom.scamp 
!swarp fr_fdobj.M82.20140228.0032.fits -c astrom.swarp -IMAGEOUT_NAME=afr_fdobj.M82.20140228.0032.fits

!ds9 afr_fdobj.M82.20140228.0032.fits &


#sex test.fits -c default.sex
#scamp test.cat -c default.scamp
#swarp test.fits -c default.swarp -IMAGEOUT_NAME=test.fits



#   $rmag<18          #ds9 catalog tool filter option
