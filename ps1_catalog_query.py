def ps1_query(name, radeg, dedeg, path, radius=1.0):
	"""
	#	SELECT STARS FROM STARS & GALAXIES (iPSF - iKron <= 0.05)
	https://outerspace.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies#
	"""
	from astroquery.vizier import Vizier
	from astropy.coordinates import Angle
	from astropy.table import Table
	from astroquery.vizier import Vizier
	import astropy.units as u
	import astropy.coordinates as coord
	import numpy as np
	comment = 'NAME'+'\t'+': '+name+'\n' \
			+ 'RA'+'\t'+': '+str(round(radeg, 3))+'\n' \
			+ 'Dec'+'\t'+': '+str(round(dedeg, 3))+'\n' \
			+ 'Radius'+'\t'+': '+str(radius)+' deg'+'\n'*2 \
			+ 'LOADING PS1 Catalog ...'+'\n'
	print(comment)
	outname	= 'ps1-'+name+'.cat'
	#	QUERY PART
	Vizier.ROW_LIMIT    = -1
	query   = Vizier.query_region(coord.SkyCoord(ra=radeg, dec=dedeg, \
								unit=(u.deg, u.deg), frame='icrs'), \
								width=str(radius*60)+'m', catalog=["II/349/ps1"])
	dum0    = query[0]
	colnames= dum0.colnames
	#	REMOVE MASKED VALUE ROW
	for col in colnames:
		indx    = np.where( dum0[col].mask == False )
		dum1    = dum0[indx]
	f_objID_bin		= []
	for i in dum1['f_objID']:
		f_objID_bin.append(bin(i)[2:])
	f_objID_bin		= np.array( f_objID_bin )
	#	SELECT POINT SOURCE & NON-VARIABLE & GOOD QUALITY STARS
	indx_sel		= []
	for j in range(len(f_objID_bin)):
		i	= f_objID_bin[j]
		#	REJECT EXTENDED SOURCES THAT CONFIMED BY PS1 & 2MASS (23, 24)
		#	REJECT QSO, RR Lyra, VARIABLE, TRANSIENT (2, 3, 4, 5, 6, 7, 8)
		#	REJECT POOR-QUALITY STACK OBJECT (30 = 0) -> not applied yet
		try:
			if (i[-23] != '1') and (i[-24] != '1') and (i[-2] != '1') and (i[-3] != '1') and (i[-4] != '1') and (i[-5] != '1') and (i[-6] != '1') and (i[-7] != '1') and (i[-8] != '1'):# and (i[0] != '1'):
				indx_sel.append(j)
		except:
			pass
	dum2	= dum1[indx_sel]
	#	SELECT STARS FROM STARS & GALAXIES (iPSF - iKron <= 0.05)
	indx_stars		= np.where( (dum2['imag'] - dum2['iKmag']) <= 0.05 )
	dum		= dum2[indx_stars]
	#	CHANGE TO GENERTAL COL. NAMES
	querytbl			= Table()
	querytbl['NUMBER']  = dum['objID']
	querytbl['RA_ICRS'] = dum['RAJ2000']
	querytbl['DE_ICRS'] = dum['DEJ2000']
	querytbl['Q']		= dum['Qual']
	querytbl['Numb_obs']= dum['Nd']
	querytbl['Numb_img']= dum['Ns']
	querytbl['gmag']    = dum['gmag']
	querytbl['e_gmag']  = dum['e_gmag']
	querytbl['rmag']    = dum['rmag']
	querytbl['e_rmag']  = dum['e_rmag']
	querytbl['imag']    = dum['imag']
	querytbl['e_imag']  = dum['e_imag']
	querytbl['zmag']    = dum['zmag']
	querytbl['e_zmag']  = dum['e_zmag']
	querytbl['ymag']    = dum['ymag']
	querytbl['e_ymag']  = dum['e_ymag']

	querytbl.write(path+'/'+outname, format='ascii.tab', overwrite=True)
	return querytbl


ps1_query('NGC3938', 178.206042, 44.120722, '.', radius=1.0)
