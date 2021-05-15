import astropy.io.ascii as ascii
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.lines import Line2D

filelist = ['B-phot-clean.csv','V-phot-clean.csv','R-phot-clean.csv']

obsname=['DOAO','LOAO','SOAO','30INCH','CCA250','Maidanak']
filtername=['B','V','R']
aperture=['auto','3"','5"','7"','psf']

vandykdat=ascii.read('vandykphot.dat')

Bdat=ascii.read(filelist[0])
Vdat=ascii.read(filelist[1])
Rdat=ascii.read(filelist[2])

galactic_extinction=[0.077,0.058,0.046] # BVR, I = 0.032 from NED, SF11

'''
 # column names
['filename', 'dateobs', 'jd',
 'zpauto', 'zpautoerr',
 'zpap3', 'zpap3err',
 'zpap5', 'zpap5err',
 'zpap7', 'zpap7err',
 'zppsf', 'zppsferr',
 'limit3', 'limit5', 'limit7',
 'starnum', 'note',
 'automag', 'automagerr',
 'ap3mag', 'ap3magerr',
 'ap5mag', 'ap5magerr',
 'ap7mag', 'ap7magerr',
 'psfmag', 'psfmagerr',
 'obs']

['MJD', 'B', 'Berr',
 'V', 'Verr',
 'R', 'Rerr',
 'unfilterd', 'unfilterderr',
 'I', 'Ierr', 'source']

'''

# Van Dyk light curve
dat=vandykdat # Van Dyk et al 2018
plt.errorbar(dat['MJD'], dat['B'],dat['Berr'],fmt='b.',alpha=0.5)
plt.errorbar(dat['MJD'], dat['V'],dat['Verr'],fmt='g.',alpha=0.5)
plt.errorbar(dat['MJD'], dat['R'],dat['Rerr'],fmt='r.',alpha=0.5)
#plt.xlim()
plt.ylim(21.5,14)
plt.xlim(57880,58150)
plt.xlabel('MJD')
plt.ylabel('Mag (VEGA)')
plt.legend(filtername)# ,'Van Dyk +2018')

#Bdat
indexreal, indexnone = np.where(Bdat['note']=='real'), np.where(Bdat['note']=='none')
Bdatreal, Bdatnone = Bdat[indexreal], Bdat[indexnone]
indexDOAO=np.where(Bdatreal['obs']=='DOAO')
indexLOAO=np.where(Bdatreal['obs']=='LOAO')
indexSOAO=np.where(Bdatreal['obs']=='SOAO')
index30inch=np.where(Bdatreal['obs']=='30inch')
indexCCA250=np.where(Bdatreal['obs']=='CCA250')
indexMaidanak=np.where(Bdatreal['obs']=='Maidanak')

plt.errorbar(Bdatreal['jd'][indexDOAO]-2400000.5,     Bdatreal['psfmag'][indexDOAO]+2,     Bdatreal['psfmagerr'][indexDOAO],fmt='bo', fillstyle='none')
plt.errorbar(Bdatreal['jd'][indexLOAO]-2400000.5,     Bdatreal['psfmag'][indexLOAO]+2,     Bdatreal['psfmagerr'][indexLOAO],fmt='b.')
plt.errorbar(Bdatreal['jd'][indexSOAO]-2400000.5,     Bdatreal['psfmag'][indexSOAO]+2,     Bdatreal['psfmagerr'][indexSOAO],    fmt='b^')
plt.errorbar(Bdatreal['jd'][index30inch]-2400000.5,   Bdatreal['psfmag'][index30inch]+2,   Bdatreal['psfmagerr'][index30inch],  fmt='bs', fillstyle='none')
plt.errorbar(Bdatreal['jd'][indexCCA250]-2400000.5,   Bdatreal['psfmag'][indexCCA250]+2,   Bdatreal['psfmagerr'][indexCCA250],  fmt='b+')
plt.errorbar(Bdatreal['jd'][indexMaidanak]-2400000.5, Bdatreal['psfmag'][indexMaidanak]+2, Bdatreal['psfmagerr'][indexMaidanak],fmt='b*')

#plt.errorbar(Bdat['jd'][indexnone]-2400000.5, Bdat['limit7'][indexnone], Bdat['zpap7err'][indexnone],lolims=Bdat['zpap7err'][indexnone])

### Bolometric Luminosity 
galactic_extinction=[0.077,0.058,0.046]
c0,c1,c2,rms = -0.083,-0.139,-0.691,0.109 # Lyman +2014 B-V 0.0~1.3 range Bolometrice correction term
# BC_B= M_bol - M_B
# BC__B= c0 + c1*(B-V) + c2 * ((B-V)**2)
# (BVcolor-0.019), B-V(galactic extinction corrected) = B-0.077 - V-0.058 = BVcolor-0.019
BolCorB= c0 + c1*((BVcolor['BVcolor']-0.019)) + c2 * (((BVcolor['BVcolor']-0.019))**2)





#Vdat
indexreal, indexnone = np.where(Vdat['note']=='real'), np.where(Vdat['note']=='none')
Vdatreal, Vdatnone = Vdat[indexreal], Vdat[indexnone]

indexDOAO=np.where(Vdatreal['obs']=='DOAO')
indexLOAO=np.where(Vdatreal['obs']=='LOAO')
indexSOAO=np.where(Vdatreal['obs']=='SOAO')
index30inch=np.where(Vdatreal['obs']=='30inch')
indexCCA250=np.where(Vdatreal['obs']=='CCA250')
indexMaidanak=np.where(Vdatreal['obs']=='Maidanak')

plt.errorbar(Vdatreal['jd'][indexDOAO]-2400000.5,     Vdatreal['psfmag'][indexDOAO]+1,     Vdatreal['psfmagerr'][indexDOAO],    fmt='go', fillstyle='none')
plt.errorbar(Vdatreal['jd'][indexLOAO]-2400000.5,     Vdatreal['psfmag'][indexLOAO]+1,     Vdatreal['psfmagerr'][indexLOAO],    fmt='g.')
plt.errorbar(Vdatreal['jd'][indexSOAO]-2400000.5,     Vdatreal['psfmag'][indexSOAO]+1,     Vdatreal['psfmagerr'][indexSOAO],    fmt='g^')
plt.errorbar(Vdatreal['jd'][index30inch]-2400000.5,   Vdatreal['psfmag'][index30inch]+1,   Vdatreal['psfmagerr'][index30inch],  fmt='gs', fillstyle='none')
plt.errorbar(Vdatreal['jd'][indexCCA250]-2400000.5,   Vdatreal['psfmag'][indexCCA250]+1,   Vdatreal['psfmagerr'][indexCCA250],  fmt='g+')
plt.errorbar(Vdatreal['jd'][indexMaidanak]-2400000.5, Vdatreal['psfmag'][indexMaidanak]+1, Vdatreal['psfmagerr'][indexMaidanak],fmt='g*')

#Rdat
indexreal, indexnone = np.where(Rdat['note']=='real'), np.where(Rdat['note']=='none')
Rdatreal, Rdatnone = Rdat[indexreal], Rdat[indexnone]

indexDOAO=np.where(Rdatreal['obs']=='DOAO')
indexLOAO=np.where(Rdatreal['obs']=='LOAO')
indexSOAO=np.where(Rdatreal['obs']=='SOAO')
index30inch=np.where(Rdatreal['obs']=='30inch')
indexCCA250=np.where(Rdatreal['obs']=='CCA250')
indexMaidanak=np.where(Rdatreal['obs']=='Maidanak')

plt.errorbar(Rdatreal['jd'][indexDOAO]-2400000.5,     Rdatreal['psfmag'][indexDOAO],     Rdatreal['psfmagerr'][indexDOAO],    fmt='ro', fillstyle='none')
plt.errorbar(Rdatreal['jd'][indexLOAO]-2400000.5,     Rdatreal['psfmag'][indexLOAO],     Rdatreal['psfmagerr'][indexLOAO],    fmt='r.')
plt.errorbar(Rdatreal['jd'][indexSOAO]-2400000.5,     Rdatreal['psfmag'][indexSOAO],     Rdatreal['psfmagerr'][indexSOAO],    fmt='r^')
plt.errorbar(Rdatreal['jd'][index30inch]-2400000.5,   Rdatreal['psfmag'][index30inch],   Rdatreal['psfmagerr'][index30inch],  fmt='rs', fillstyle='none')
plt.errorbar(Rdatreal['jd'][indexCCA250]-2400000.5,   Rdatreal['psfmag'][indexCCA250],   Rdatreal['psfmagerr'][indexCCA250],  fmt='r+')
plt.errorbar(Rdatreal['jd'][indexMaidanak]-2400000.5, Rdatreal['psfmag'][indexMaidanak], Rdatreal['psfmagerr'][indexMaidanak],fmt='r*')




legend_elements = [Line2D([], [], marker='o', color='k', label='DOAO',ls='none', fillstyle='none'),
				   Line2D([], [], marker='.', color='k', label='LOAO',ls='none'),
				   Line2D([], [], marker='^', color='k', label='SOAO',ls='none'),
				   Line2D([], [], marker='s', color='k', label='30INCH',ls='none', fillstyle='none'),
				   Line2D([], [], marker='+', color='k', label='CCA250',ls='none'),
				   Line2D([], [], marker='*', color='k', label='Maidanak',ls='none'),
				  ]
plt.legend(handles=legend_elements)


#plt.errorbar(Bdat['jd'][indexnone]-2400000.5, Bdat['limit7'][indexnone], Bdat['zpap7err'][indexnone],lolims=Bdat['zpap7err'][indexnone])

plt.text(58140,18,'B+2',color='b',weight='bold',fontsize=12)
plt.text(58140,18.5,'V+1',color='g',weight='bold',fontsize=12)
plt.text(58140,19,'R',color='r',weight='bold',fontsize=12)

#plt.vlines(57898.99,14,22) # discovery date 57898.99 by Ron Arbour
plt.xlim(57880,58190)
plt.ylim(23,14)
plt.ylabel('Mag (Vega)')
plt.xlabel('MJD')
plt.title('SN 2017ein Light Curve')
plt.savefig('sn2017ein-LC.png')
plt.close()




#Bdat
#Vdat
#Rdat


BVcolor=ascii.read('BV-phot-clean-color.csv')
VRcolor=ascii.read('VR-phot-clean-color.csv')
BRcolor=ascii.read('BR-phot-clean-color.csv')

Vmax=2457913.1 # Van Dyk 2018
plt.subplots_adjust(hspace=0)
plt.subplot(311)
plt.errorbar(BVcolor['jd_1']-Vmax ,BVcolor['BVcolor'],BVcolor['BVcolorerr'], marker='o',fillstyle='none',color='k',ls='')
plt.ylabel('B-V')

plt.subplot(312)
plt.errorbar(VRcolor['jd_1']-Vmax ,VRcolor['VRcolor'],VRcolor['VRcolorerr'], marker='s',fillstyle='none',color='k',ls='')
plt.ylabel('V-R')

plt.subplot(313)
plt.errorbar(BRcolor['jd_1']-Vmax ,BRcolor['BRcolor'],BRcolor['BRcolorerr'], marker='^',fillstyle='none',color='k',ls='')
plt.ylabel('B-R')
plt.xlabel('Day since V Maximum')

plt.savefig('SN2017ein-color.png')

plt.close()


galactic_extinction=[0.077,0.058,0.046]
c0,c1,c2,rms = -0.083,-0.139,-0.691,0.109 # Lyman +2014 B-V 0.0~1.3 range Bolometrice correction term
# BC_B= M_bol - M_B
# BC__B= c0 + c1*(B-V) + c2 * ((B-V)**2)
# (BVcolor-0.019), B-V(galactic extinction corrected) = B-0.077 - V-0.058 = BVcolor-0.019
BolCorB= c0 + c1*((BVcolor['BVcolor']-0.019)) + c2 * (((BVcolor['BVcolor']-0.019))**2)


M_B=BVcolor['psfmag_1']-31.75-0.077
plt.plot( BVcolor['jd_1']-Vmax, Mbol_B,'bo')
plt.title('Bolometric Light Curve')
plt.xlabel('Day MJD')
plt.ylable('')
plt.show()
