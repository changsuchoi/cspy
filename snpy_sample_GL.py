import os, sys
import numpy as np
import snpy
​
#datpath = '/data1/SN2019ein/work/lc/'
# Define filter BVRIJHK
'''
B = snpy.filters.filter('B',file='/data1/SN2019ein/work/lc_v2/filters/Generic_Johnson.B.dat', zp=6.293e-9, comment='Johnson/Cousin B, Vega, erg/cm2/s/A')
V = snpy.filters.filter('V',file='/data1/SN2019ein/work/lc_v2/filters/Generic_Johnson.V.dat', zp=3.575e-9, comment='Johnson/Cousin V, Vega, erg/cm2/s/A')
R = snpy.filters.filter('R',file='/data1/SN2019ein/work/lc_v2/filters/Generic_Johnson.R.dat', zp=1.882e-9, comment='Johnson/Cousin R, Vega, erg/cm2/s/A')
I = snpy.filters.filter('I',file='/data1/SN2019ein/work/lc_v2/filters/Generic_Johnson.I.dat', zp=9.329e-10, comment='Johnson/Cousin I, Vega, erg/cm2/s/A')
J = snpy.filters.filter('J',file='/data1/SN2019ein/work/lc_v2/filters/UKIRT_WFCAM.J_filter.dat', zp=2.960e-10, comment='WFCAM J, Vega, erg/cm2/s/A')
H = snpy.filters.filter('H',file='/data1/SN2019ein/work/lc_v2/filters/UKIRT_WFCAM.H_filter.dat', zp=1.142e-10, comment='WFCAM H, Vega, erg/cm2/s/A')
K = snpy.filters.filter('K',file='/data1/SN2019ein/work/lc_v2/filters/UKIRT_WFCAM.K_filter.dat', zp=3.808e-11, comment='WFCAM K, Vega, erg/cm2/s/A')
'''

Bd=ascii.read('SN2018kp-B.dat')
Vd=ascii.read('SN2018kp-V.dat')
Rd=ascii.read('SN2018kp-R.dat')
Id=ascii.read('SN2018kp-I.dat')
gd=ascii.read('SN2018kp-g.dat')
rd=ascii.read('SN2018kp-r.dat')
iid=ascii.read('SN2018kp-i.dat')

datpath = '/data1/SN2019ein/work/lc_v2/'
s = snpy.get_sn(datpath+'SN2019ein_snpy.txt', sql=None)
​
# 2nd try
model = 'EBV_model'
s.choose_model(model, stype='dm15')
​
s.restbands['BANDI']  = 'Bi'
s.restbands['VANDI']  = 'Vi'
s.restbands['RANDI']  = 'Ri'
s.restbands['IANDI']  = 'Ii'
s.restbands['WFCAMJ'] = 'WFCAMJ'
s.restbands['WFCAMH'] = 'WFCAMH'
s.restbands['WFCAMK'] = 'WFCAMK'
​
#s.fit( mangle=True, kcorr=True, k_stretch=True) # 20201009 -> success!
#DM = 33.003  +/-  0.012  +/- 0.121 (sys)
#dm15 = 1.430  +/-  0.012  +/- 0.060 (sys)
#EBVhost = 0.114  +/-  0.009  +/- 0.060 (sys)
#Tmax = 58619.315  +/-  0.001  +/- 0.340 (sys)
​
s.fit( mangle=False, kcorr=True, k_stretch=True) # 20201009_2 -> success!
#DM = 33.041  +/-  0.014  +/- 0.119 (sys)
#dm15 = 1.412  +/-  0.015  +/- 0.060 (sys)
#EBVhost = 0.098  +/-  0.011  +/- 0.060 (sys)
#Tmax = 58619.387  +/-  0.057  +/- 0.340 (sys)
s.get_max(['Bi','Vi','Ri','Ii'], deredden=False)
​
​
​
band1 = ['BANDI','VANDI','RANDI', 'IANDI', 'WFCAMJ', 'WFCAMH', 'WFCAMK']
band2 = ['Bs','V','Rs', 'Is','WFCAMJ', 'WFCAMH', 'WFCAMK']
band3 = ['Bs','Vs','Rs','Is']
band4 = ['Bs', 'Vs']
band5= ['Bi','Vi','Ri','Ii','WFCAMJ', 'WFCAMH', 'WFCAMK']
band6= ['WFCAMJ', 'WFCAMH', 'WFCAMK']
​
s.fit( mangle=True, kcorr=True, k_stretch=True)
s.fit(band6)
s.getEBVgal()
print("MW extinction (E(B-V)_MW) =", s.EBVgal)
​
#s.plot( yrange=(19.5, 7), title='SN 2019ein', single=True, offset=True, legend=True, fsize=15, linewidth=1.5, symbols={'BANDI':'s', 'VANDI':'s', 'RANDI':'s', 'IANDI':'s', 'WFCAMJ':'s', 'WFCAMH':'s', 'WFCAMK':'s'}, colors={'BANDI':'royalblue', 'VANDI':'seagreen', 'RANDI':'orange', 'IANDI':'orangered', 'WFCAMJ':'tomato', 'WFCAMH':'red', 'WFCAMK':'brown'}, flux=False, epoch=True, outfile='SN2019ein.snpy.lc.pdf')
​
​
# 1st try
#s.fit(bands=['WFCAMK', 'IANDI', 'RANDI', 'WFCAMK', 'WFCAMJ', 'VANDI', 'BANDI'], kcor=True, reset_kcorrs=False, mangle=True )
#s.fit(bands=['VANDI', 'BANDI'], kcor=True, reset_kcorrs=True, mangle=True )
s.choose_model('EBV_model', stype='dm15')
s.fit(['BANDI','VANDI'])
print("The value of dm15 is",s.dm15,"+/-",s.e_dm15)
s.plot(outfile='SN2019ein_EBV.pdf')
s.summary()
print("Ia distance:",s.DM,"+/-",s.e_DM,"\nHubble distance",s.get_distmod())
'''
s.choose_model('EBV_model', stype='dm15')
s.fit(['VANDI','RANDI'])
print("The value of dm15 is",s.dm15,"+/-",s.e_dm15)
s.plot(outfile='SN2019ein_EVR.pdf')
s.summary()
print("Ia distance:",s.DM,"+/-",s.e_DM,"\nHubble distance",s.get_distmod())
'''
​
​
# B-V color
​
#BV_MJD, BV, BVe, BVflag= s.get_color('BANDI','VANDI')
#plt.scatter(BV_MJD, BV)
#plt.show()
'''
s.choose_model('max_model', stype='dm15')
s.fit(['BANDI','VANDI'])
s.plot(outfile='SN2019ein_max.pdf')
s.summary()
​
s.save('SN2019ein_maxmodel.snpy')
​
print("B-band data covers MJD=", s.B.MJD.min(), "to MJD=", s.B.MJD.max())
'''
