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

def ps1_Tonry(intbl, name):
	'''
	PS1 -> Johnson/COusins [Vega] -> [AB]	(Tonry+12)
	#   Vega - AB Magnitude Conversion (Blanton+07)
	U       : m_AB - m_Vega = 0.79
	B       : m_AB - m_Vega =-0.09
	V       : m_AB - m_Vega = 0.02
	R       : m_AB - m_Vega = 0.21
	I       : m_AB - m_Vega = 0.45
	J       : m_AB - m_Vega = 0.91
	H       : m_AB - m_Vega = 1.39
	K       : m_AB - m_Vega = 1.85
	'''
	import numpy as np
	from astropy.table import Table
	from astropy.io import ascii
	#	REJECT BAD QUALITY
	intbl	= intbl[	(intbl['gmag']>12)	&
						(intbl['rmag']>12)	&
						(intbl['imag']>12)	&
						(intbl['zmag']>12)	&
						(intbl['ymag']>12)]


	outfile	= 'ps1-Tonry-'+name+'.cat'
	Q		= intbl['Q']
	indx    = np.where(Q < 128)

	intbl	= intbl[indx]
	Q		= intbl['Q']

	name    = intbl['NUMBER']
	ra, de  = intbl['RA_ICRS'],    intbl['DE_ICRS']

	g, ger  = intbl['gmag'],	intbl['e_gmag']
	r, rer  = intbl['rmag'],	intbl['e_rmag']
	i, ier  = intbl['imag'],	intbl['e_imag']
	z, zer  = intbl['zmag'],	intbl['e_zmag']
	y, yer  = intbl['ymag'],	intbl['e_ymag']
	#	TRANSF. ERROR FOR B CONST. TERMS
	Bsig, Vsig, Rsig, Isig	= 0.034, 0.012, 0.01, 0.016
	#	COLOR TERM
	gr		= intbl['gmag']-intbl['rmag']
	grer	= sqsum(intbl['e_gmag'], intbl['e_rmag'])
	#	CONVERT TO B
	B0		= 0.213
	B1		= 0.587
	B		= B0 + B1*gr + intbl['gmag'] - 0.09
	Ber		= sqsum( Bsig, sqsum(B1*grer, intbl['e_gmag']) )
	#	CONVERT TO V
	B0		= 0.006
	B1		= 0.474
	V		= B0 + B1*gr + intbl['rmag'] + 0.02
	Ver	    = sqsum( Bsig, sqsum(B1*grer, intbl['e_rmag']) )
	#	CONVERT TO R
	B0		=-0.138
	B1		=-0.131
	R		= B0 + B1*gr + intbl['rmag'] + 0.21
	Rer		= sqsum( Rsig, sqsum(B1*grer, intbl['e_rmag']) )
	#	CONVERT TO I
	B0		=-0.367
	B1		=-0.149
	I		= B0 + B1*gr + intbl['imag'] + 0.45
	Ier		= sqsum( Isig, sqsum(B1*grer, intbl['e_imag']) )
	outbl	= Table([name, ra, de, g, ger, r, rer, i, ier, z, zer, y, yer, B, Ber, V, Ver, R, Rer, I, Ier, Q],
					names=['name', 'ra', 'dec', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'y', 'yerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'Q'])
	outtbl0	= Table([name, ra, de, g, ger, r, rer, i, ier, z, zer, y, yer, B, Ber, V, Ver, R, Rer, I, Ier, Q],
					names=['#name', 'ra', 'dec', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'y', 'yerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'Q'])
	outtbl.write(outfile, format='ascii.tab', overwrite=True)
	return outbl, outfile
