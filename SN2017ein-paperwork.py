import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
from astropy.table import Table, Column, MaskedColumn
import numpy as np


Bphot, Vphot, Rphot = ascii.read('B-phot-clean.csv'), ascii.read('V-phot-clean.csv'), ascii.read('R-phot-clean.csv')
Bphot, Vphot, Rphot = ascii.read('B-phot.csv'), ascii.read('V-phot.csv'), ascii.read('R-phot.csv')
## new table for paper in latex format

phot=[Bphot, Vphot, Rphot]


mjd=phot[i]['jd']-2400000.5
mag=phot[i]['psfmag']
magerr=phot[i]['psfmagerr']
obs=phot[i]['obs']
limitmag=phot[i]['limit5']
note=phot[i]['note']

obslist=['30inch', 'DOAO', 'LOAO', 'Maidanak', 'SOAO','CCA250']

## each Tables B
realid=np.where(phot[0]['note']=='real')
Bmjd=Column(phot[0]['jd'][realid]-2400000.5,format='5.2f')
Bmag=Column(phot[0]['psfmag'][realid],format='2.3f')
Bmagerr=Column(phot[0]['psfmagerr'][realid],format='1.3f')
Bobs=Column(phot[0]['obs'][realid])
Blimitmag=Column(phot[0]['limit5'][realid],format='2.3f')
Blimitmagerr=Column(phot[0]['zpap5err'][realid],format='2.3f')
Bnote=Column(phot[0]['note'][realid])

Btable=Table([Bmjd,Bmag,Bmagerr,Bobs,Blimitmag,Blimitmagerr,Bnote],names=['MJD','Mag','Error','Observatory','LimitMag','LimitMagErr','Note'])
Btable.sort('MJD')
ascii.write(Btable,'Btable.txt',format='latex')
ascii.write(Btable,'Btbl.dat')
## each Tables V
realid=np.where(phot[1]['note']=='real')
Vmjd=Column(phot[1]['jd'][realid]-2400000.5,format='5.2f')
Vmag=Column(phot[1]['psfmag'][realid],format='2.3f')
Vmagerr=Column(phot[1]['psfmagerr'][realid],format='1.3f')
Vobs=Column(phot[1]['obs'][realid])
Vlimitmag=Column(phot[1]['limit5'][realid],format='2.3f')
Vlimitmagerr=Column(phot[1]['zpap5err'][realid],format='2.3f')
Vnote=Column(phot[1]['note'][realid])

Vtable=Table([Vmjd,Vmag,Vmagerr,Vobs,Vlimitmag,Vlimitmagerr,Vnote],names=['MJD','Mag','Error','Observatory','LimitMag','LimitMagErr','Note'])
Vtable.sort('MJD')
ascii.write(Vtable,'Vtable.txt',format='latex')
ascii.write(Vtable,'Vtbl.dat')


## each Tables R
realid=np.where(phot[2]['note']=='real')
Rmjd=Column(phot[2]['jd'][realid]-2400000.5,format='5.2f')
Rmag=Column(phot[2]['psfmag'][realid],format='2.3f')
Rmagerr=Column(phot[2]['psfmagerr'][realid],format='1.3f')
Robs=Column(phot[2]['obs'][realid])
Rlimitmag=Column(phot[2]['limit5'][realid],format='2.3f')
Rlimitmagerr=Column(phot[2]['zpap5err'][realid],format='2.3f')
Rnote=Column(phot[2]['note'][realid])

Rtable=Table([Rmjd,Rmag,Rmagerr,Robs,Rlimitmag,Rlimitmagerr,Rnote],names=['MJD','Mag','Error','Observatory','LimitMag','LimitMagErr','Note'])
Rtable.sort('MJD')
ascii.write(Rtable,'Rtable.txt',format='latex')
ascii.write(Rtable,'Rtbl.dat')
## All combined tables BVR mag and error only

bvr=ascii.read('BVRrealphot.txt')
ascii.write(bvr,'ptable.txt',format='latex')

## plot
#early
plt.clf()
plt.subplot(1,2,1)


#late
plt.subplot(1,2,2)

plt.tight_layout()
plt.title()

plt.savefig('SN2017einLC-paper.eps')

plt.close()


tbl=ascii.read('bvr-tbl.txt')
tbl.filled('\nodata')
mjd=Column(tbl['MJD'],format='5.3f')
tbl.remove_column('MJD')
tbl.add_column(mjd,0)


ascii.write(tbl,'tbl-latex.txt',format='latex',overwrite=True)
