# Calculate MJD and modifiy its header to have MJD information
# The 1st code from IMSNG team collaboration.
# Usage : putMJD(inim,'C2') ; for LOAO
# Usage : putMJD(inim,'C5') ; for 30INCH
# last modified : 2018/10/28


import astropy.io.fits as fits
import astropy.io.ascii as ascii
from astropy.time import Time

# obsname=ascii.read('/home/changsu/Documents/IMSNG/IMSNG_obs_DB.xlsx')

def putMJD(inim,case):
	hdr=fits.getheader(inim)
	if   case == 'C1' : 
		mjdval=Time(hdr['JD'],format='jd',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
		print ('header updated', fits.getheader(inim)['MJD'])
	elif case == 'C2' : 
		mjdval=Time(hdr['J_DATE'],format='jd',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
		print ('header updated', fits.getheader(inim)['MJD'])
	elif case == 'C3' :
		mjdval=Time(hdr['DATE-OBS'],format='isot',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
		print ('header updated', fits.getheader(inim)['MJD'])
	elif case == 'C4' : 
		datetimestr=hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
		mjdval = Time(datetimestr,format='isot',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
	elif case =='C5' : 
		datetimestr=hdr['DATE-OBS']+'T'+hdr['UT']
		mjdval = Time(datetimestr,format='isot',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	

	else : print ('Not in known cases.')	
	
'''
CASE
C1		JD							SOAO, DOAO, BOAO, SAO1m, CCA250, T52(SNUCAMII)
C2		J_DATE						LOAO
C3		DATE-OBS					MAO(FLI), iTel(except T52)
C4		'DATE-OBS'+'T'+'TIME-OBS'	MAO(SNUCAM)			
C5		'DATE-OBS'+'T'+'UT'			30inch		
'''


def putmjd(i):
	try : 
		hdr=fits.getheader(i)
		mjdval=Time(hdr['DATE-OBS'],format='isot',scale='utc').mjd
		fits.setval(i,'MJD',value=mjdval,comment='MJD appended')	
		print ('header updated', fits.getheader(i)['MJD'])

	except : print(i,'DATE-OBS keyword is missing')
