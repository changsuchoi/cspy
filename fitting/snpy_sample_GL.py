# SNPY fitting
# https://users.obs.carnegiescience.edu/cburns/SNooPyDocs/html/index.html
# https://csp.obs.carnegiescience.edu/data/snpy/documentation

import os, sys
import numpy as np
import snpy
​from astropy.table import Table, vstack

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
dlist=[Bd, Vd, Rd, Id, gd, rd, iid]

Bd0=Bd[Bd['DET']=='Y']
Bd0.sort('MJD')
Vd0=Vd[Vd['DET']=='Y']
Vd0.sort('MJD')
Rd0=Rd[Rd['DET']=='Y']
Rd0.sort('MJD')
Id0=Id[Id['DET']=='Y']
Id0.sort('MJD')
gd0=gd[gd['DET']=='Y']
gd0.sort('MJD')
rd0=rd[rd['DET']=='Y']
rd0.sort('MJD')
iid0=iid[iid['DET']=='Y']
iid0.sort('MJD')

vsd=vstack([Bd0['MJD','MAG','MAGERR'], Vd0['MJD','MAG','MAGERR'],
        Rd0['MJD','MAG','MAGERR'], Id0['MJD','MAG','MAGERR'],
        gd0['MJD','MAG','MAGERR'], rd0['MJD','MAG','MAGERR'], iid0['MJD','MAG','MAGERR']])
vsd.write('SN2018kp_snpy.dat',format='ascii.commented_header', overwrite=True)
# edit file incuding filter names

# filter set check
snpy.fset.list_filters()
'''
APO
	SDSS
		'u_s':  sloan u at APO
		'g_s':  sloan g at APO
		'r_s':  sloan r at APO
		'i_s':  sloan i at APO
		'z_s':  sloan z at APO
CTIO
	0.9meter
		'Uctio':  CTIO 0.9 m U band
		'Bctio':  CTIO 0.9 m B band
		'Vctio':  CTIO 0.9 m V band
		'Rctio':  CTIO 0.9 m R band
		'Ictio':  CTIO 0.9 m I band
	ANDICAM
		'BANDI':  CTIO 1.3 m ANDICAM B band
		'VANDI':  CTIO 1.3 m ANDICAM V band
		'RANDI':  CTIO 1.3 m ANDICAM R band
		'IANDI':  CTIO 1.3 m ANDICAM I band
		'JANDI':  CTIO 1.3 m ANDICAM J band
		'HANDI':  CTIO 1.3 m ANDICAM H band
		'KANDI':  CTIO 1.3 m ANDICAM K band

standard
	Persson
		'J_K':  J-band for Swope at LCO
		'H_K':  H-band for Swope at LCO
		'K_K':  K-band for Swope at LCO
	stritzinger
		'Us':  Kron-Cousins U filter based on Stritzinger et al. 2005
		'Bs':  Kron-Cousins B filter based on Stritzinger et al. 2005
		'Vs':  Kron-Cousins B filter based on Stritzinger et al. 2005
		'Rs':  Kron-Cousins R filter based on Stritzinger et al. 2005
		'Is':  Kron-Cousins I filter based on Stritzinger et al. 2005
'''

datpath = '/data1/SN2019ein/work/lc_v2/'
datpath = '/data7/cschoi/sngal/NGC3367/phot/cut/fitting/'
# open new terminal and type 'snpy' starting with new environmnet, to avoid sql connetion error

s = get_sn('SN2018kp_snpy.txt', sql=None)
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
​s.get_max(['Bs','Vs','Rs','Is','g_s','r_s','i_s'])
​
​
band1 = ['BANDI','VANDI','RANDI', 'IANDI', 'WFCAMJ', 'WFCAMH', 'WFCAMK']
band2 = ['Bs','V','Rs', 'Is','WFCAMJ', 'WFCAMH', 'WFCAMK']
band3 = ['Bs','Vs','Rs','Is']
band4 = ['Bs', 'Vs']
band5= ['Bi','Vi','Ri','Ii','WFCAMJ', 'WFCAMH', 'WFCAMK']
band6= ['WFCAMJ', 'WFCAMH', 'WFCAMK']
​band7=['Bs','Vs','Rs']
bands=['Bs','Vs','Rs']
ss.fit( mangle=1, kcorr=0, k_stretch=0)
ss.fit(band)
s.getEBVgal()
print("MW extinction (E(B-V)_MW) =", s.EBVgal)

'''​
s.plot( yrange=(19.5, 7), title='SN 2019ein', single=True, offset=True,
        legend=True, fsize=15, linewidth=1.5,
        symbols={'BANDI':'s', 'VANDI':'s', 'RANDI':'s', 'IANDI':'s', 'WFCAMJ':'s', 'WFCAMH':'s', 'WFCAMK':'s'},
        colors={'BANDI':'royalblue', 'VANDI':'seagreen', 'RANDI':'orange', 'IANDI':'orangered',
                'WFCAMJ':'tomato', 'WFCAMH':'red', 'WFCAMK':'brown'},
        flux=False, epoch=True, outfile='SN2019ein.snpy.lc.pdf')
​'''
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




# param print
for param in s.parameters:
    print("{} = {} +/- {}".format(param, s.parameters[param], s.errors[param]))

from snpy  import get_sn
s = get_sn('https://sne.space/astrocats/astrocats/supernovae/output/json/SN2006ax.json')
bands = ['u','g','r','i','B','V','Y','J','H']
bands = ['g_s','r_s','i_s','Bs','Vs','Rs','Is']
s.filter_order = bands   # Only plot these filters
s.plot()
for band in bands:
    s.data[band].curve_fit(method='spline')    # often the defaults are fine
s.plot()                      # Check that the curves look good

s.choose_model('EBV_model2')
s.fit(bands)
DM = s.get_distmod(H0=70.0)
print("E(B-V) = {}, distance = {}".format(s.EBVhost, DM))


'''
Step 3: Fit Bolometric Light Curve (Direct Method)
The direct method simply takes the magnitude in each filter (interpolating if necessary) and converts to flux.
These fluxes are then integrated over wavelength from the bluest to reddest filter.
The reddening (EBVhost) is removed from the flux assuming a reddening law (Rv) and the distance (DM) is used to convert to luminosity.
A good value of Rv is 2.0 for SNe Ia, but only matters for relatively large E(B−V).
Lastly, we tell the function to interpolate using the spline (well, GP, actually, but spline what generically call the interpolating function).
'''
t,Lbol,filts,limits = s.bolometric(bands, method='direct', DM=DM,
                                  EBVhost=s.EBVhost, Rv=3.1,
                                  interpolate='spline', refband="V", interp_all=False)
print(limits)

from matplotlib import pyplot as plt
fig,ax = plt.subplots()
ax.plot(t, Lbol, 'o')
ax.set_xlabel('Epoch (days)')
ax.set_ylabel('Bolometric Luminosity (erg/s)')

'''
Step 3: Fit Bolometric Light Curve (SED Method)
If the object in question is a type Ia supernova,
then another approach is to fit a model of the SN spectral energy distribution (SED)
to the observed photometry and integrate that over wavelength.
For the model, SNooPy will use the Hsiao et al. (2007) SED and
use the same color-matching mechanics used to compute k-corrections to warp it to fit the observed photometry.
'''
t2,Lbol2,filts2,limits2 = s.bolometric(bands, method='SED', DM=DM,
                                      EBVhost=s.EBVhost, Rv=3.1,
                                      interpolate='spline', refband='Vs', interp_all=False,
                                      lam1=3584.0, lam2=15863.0)
fig,ax = plt.subplots()
ax.plot(t, Lbol, 'o', label='direct')
ax.plot(t2,Lbol2,'o', label='SED')
ax.set_xlabel('Epoch (days)')
ax.set_ylabel('Bolometric Luminosity (erg/s)')


#s = get_sn('SN1981D.txt')
s.summary()
s.plot()
s.restbands
s.restbands['B'] = 'Bs'
s.restbands['Vs'] = 'V'
s.restbands['Rs'] = 'R'
s.restbands['Is'] = 'I'
s.restbands['g_s'] = 'g'
s.restbands['r_s'] = 'r'
s.restbands['i_s'] = 'i'
s.fit(['Bs','Vs','Rs'])
s.fit([’g_m’,’r_m’,’i_m’,’Jc’], dm15=s.dm15, Tmax=s.Tmax)
s.fit(['Bs','Vs','Rs'], dm15=s.dm15, Tmax=s.Tmax)
s.B.template(method=’spline’, task=0, s=len(s.B.mag))
s.B.template(method=’spline’, task=1, s=len(s.B.mag) - sqrt(2*len(s.B.mag)))
s.B.template(method=’chebyshev’, n=5)
s.summary()
s.save('SN1981D.snpy')
s.dump_lc()

'''
3.8.5    Interactive Fitting
A new feature is the ability to interactively fit the data using the matplotlib library.
Simply specify interactive=True when you call the template()member function.
A window that looks like figure XXX will pop up, showing you the light-curve, the current fit,
and the residuals from the fit.
You will then have access to certain keystrokes depending on the fitting method.In all cases,
you can use the following keys:
•’x’:  The point closest to the mouse pointer will be masked
(red X will cover it and itwon’t contribute to the fit).
If you use ’x’ again near this point, you un-mask the data.
•’r’:  re-fit the plot (if necessary) and re-draw the plot.
•’c’:  (re)compute the light-curve parameters and plot them on the light-curve.
•’q’:  quit the interactive plotting and close the graph.
•’ ?’:  quick help on what keystrokes are availableWhen fitting Gaussian Processes,
the following keys can be used:
•’s’ and ’S’: decrease and increase the scale by 10%, respectively.
•’a’ and ’A’: decrease and increase the amplitude by 10%, respectively.
•’d’ and ’D’: decrease and increase the degree of differentiability by 1
When fitting with Dierckx splines, the following keys can be used:
•’a’: Add a knot point at the cursor position.  Note: this will change the task parameter to -1.
•’d’:  Delete knot point closest to the cursor position.
Note:  this will change thetaskparameter to -1.
•’m’:  Move  the  knot  point  closest  to  the  cursor  position  to  a  new  position  (hit  ’m’ again).
Note:  this will change the task parameter to -1.
When fitting with polynomials, the following keys can be used:
•’n’ and ’N’: decrease and increase the order of the polynomial by 1
•’m’:  specify range over which to fit (press ’m’ at beginning and again at end).
Pressing ’m’ twice in the same location will reset to default range.
'''




#====================================================

# model EBVmodel
# model EBVmodel
# model maxmodel
# model colormodel

models = ['max_model', 'EBV_model', 'EBV_model2', 'color_model']
models= [ 'EBV_model', 'EBV_model2',
        'max_model', 'max_model2',
        'Rv_model',
        'color_model',
        'SALT_model',
        'MLCS_model']
model=models[0]
s = get_sn('SN2018kp_snpy.txt', sql=None)
s.choose_model(model, stype='dm15')

'''
s.restbands['Bs'] = 'B'
s.restbands['Vs'] = 'V'
s.restbands['Rs'] = 'R'
s.restbands['Is'] = 'I'
s.restbands['g_s'] = 'g'
s.restbands['r_s'] = 'r'
s.restbands['i_s'] = 'i'
'''
bands=['Bs','Vs','Rs','Is','g_s','r_s','i_s']
bands=['Bs','Vs','Rs']

s.filter_order = bands   # Only plot these filters
s.plot()
for band in bands:
    s.data[band].curve_fit(method='spline')    # often the defaults are fine
s.plot()

s.fit( mangle=False, dokcorr=0, k_stretch=True)
s.get_max(['Bs','Vs','Rs','Is','g_s','r_s','i_s'])
# bands=['Bs','Vs','Rs']
s.fit( mangle=True, dokcorr=0, k_stretch=True)
s.fit(bands)
# LC fit band ['Bs', 'Vs', 'Rs']
# dm15 DM Tmax, MaxMag
# plot : spline, LC fit, color, Bolometric
#

s.choose_model('EBV_model2',stype='dm15')
s.fit(bands,mangle=True, dokcorr=0, k_stretch=True)
DM = s.get_distmod(H0=70.0)
print("E(B-V) = {}, distance = {}".format(s.EBVhost, DM))


t,Lbol,filts,limits = s.bolometric(['Bs','Vs','Rs'], method='direct', DM=DM,
                                  EBVhost=s.EBVhost, Rv=1.766 ,
                                  interpolate='spline', refband="Vs", interp_all=False)
print(limits)

from matplotlib import pyplot as plt
fig,ax = plt.subplots()
ax.plot(t, Lbol, 'o')
ax.set_xlabel('Epoch (days)')
ax.set_ylabel('Bolometric Luminosity (erg/s)')

'''
Step 3: Fit Bolometric Light Curve (SED Method)
If the object in question is a type Ia supernova,
then another approach is to fit a model of the SN spectral energy distribution (SED)
to the observed photometry and integrate that over wavelength.
For the model, SNooPy will use the Hsiao et al. (2007) SED and
use the same color-matching mechanics used to compute k-corrections to warp it to fit the observed photometry.
'''
t2,Lbol2,filts2,limits2 = s.bolometric(bands, method='SED', DM=DM,
                                      EBVhost=s.EBVhost, Rv=1.766,
                                      interpolate='spline', refband='Vs', interp_all=False,
                                      lam1=3584.0, lam2=15863.0)
fig,ax = plt.subplots()
ax.plot(t, Lbol, 'o', label='direct')
ax.plot(t2,Lbol2,'o', label='SED')
ax.set_xlabel('Epoch (days)')
ax.set_ylabel('Bolometric Luminosity (erg/s)')



s.choose_model('EBVmodel', stype='dm15')
s.choose_model('EBVmodel2', stype='dm15')
s.choose_model('max_model', stype='dm15')
s.choose_model('color_model', stype='st') #color_model 'st' only

s.fit( mangle=True, dokcorr=True, k_stretch=True)
s.summary()
for param in s.parameters:
    print("{} = {} +/- {}".format(param, s.parameters[param], s.errors[param]))

s.fit( bands, mangle=True, dokcorr=True, k_stretch=True)
s.summary()
for param in s.parameters:
    print("{} = {} +/- {}".format(param, s.parameters[param], s.errors[param]))


# from online documents
# https://users.obs.carnegiescience.edu/cburns/SNooPyDocs/html/overview.html#light-curves
# spline fit

# color
s.plot_color('Bs','Vs')
s.plot_color('Vs','Rs')
s.plot_color('Bs','Rs')

# MCMC
s.fitMCMC(bands, R_V="N,2.3,0.9")

#lira
s.lira('Bs','Vs',dokcorr=0,plot=1,interpolate=1)
'''
Warning:  can't k-correct spline fits, you're getting observed maxima!
Vmax occurred at 58160.726562
Slope of (B-V) vs. t-Tvmax was -0.014341(0.026924)
median E(B-V) = 0.415282    1.49*mad(E(B-V)) = 0.021416
Out[31]:
(0.41528224676904296,
 0.021416154286220346,
 -0.014340638785003572,
 0.02692357509401204)
'''
# plot kcorrection
s.plot_kcorrs()

# data table and plot saving

# module check
