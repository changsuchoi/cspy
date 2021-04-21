import astropy.io.ascii as ascii
from astropy.table import Table, hstack
import pandas as pd
from astropy.table import Table, join

Bdata=ascii.read('SN2018kp-B.dat')
Rdata=ascii.read('SN2018kp-R.dat')
Vdata=ascii.read('SN2018kp-V.dat')
Idata=ascii.read('SN2018kp-I.dat')
gdata=ascii.read('SN2018kp-g.dat')
rdata=ascii.read('SN2018kp-r.dat')
idata=ascii.read('SN2018kp-i.dat')

Bmjd,Vmjd,Rmjd=[],[],[]
for n in Bdata['MJD']:
	Bmjd.append(round(n,2))
	#Bmjd.append(float(format(n,'.2f')))
for n in Vdata['MJD']:
	Vmjd.append(round(n,2))
for n in Rdata['MJD']:
	Rmjd.append(round(n,2))

Bdata['MJD']=Bmjd
Vdata['MJD']=Vmjd
Rdata['MJD']=Rmjd
Bpd=Bdata.to_pandas()
Vpd=Vdata.to_pandas()
Rpd=Rdata.to_pandas()

#join(Bdata,Vdata,Rdata,keys='MJD',join_type='outer')

BV=join(Bdata,Vdata,keys='MJD',join_type='outer',table_names=['B', 'V'],uniq_col_name='{table_name}_{col_name}')
BVR=join(BV,Rdata,keys='MJD',join_type='outer',table_names=['', 'R'],uniq_col_name='{table_name}_{col_name}')
BVRI=join(BVR,Idata,keys='MJD',join_type='outer',table_names=['R','I'],uniq_col_name='{table_name}_{col_name}')
BVR.write('BVR.dat',format='ascii.commented_header',overwrite=True)
