#code for plot of SN2017ein Light Curve


import astropy.io.ascii as ascii
import matplotlib.pyplot as plt 
import numpy as np 
import matplotlib.patches as mpatches


#file read
B1data=ascii.read('B-phot-cut.dat')
V1data=ascii.read('V-phot-cut.dat')
R1data=ascii.read('R-phot-cut.dat')

firstdetect= 2457899.27211
obsnames=['30inch', 'DOAO', 'LOAO', 'Maidanak', 'SOAO','CCA250']
obsmarkers=['+','8','o','x','^','D']
colnames=['filename', 'dateobs', 'jd',
 'zpauto', 'zpautoerr', 'zpap3', 'zpap3err',
 'zpap5', 'zpap5err', 'zpap7', 'zpap7err',
 'zppsf', 'zppsferr',
 'limit3', 'limit5', 'limit7',
 'starnum', 'note',
 'automag', 'automagerr', 'ap3mag', 'ap3magerr',
 'ap5mag', 'ap5magerr', 'ap7mag', 'ap7magerr',
 'psfmag', 'psfmagerr',
 'obs']

# Rdata
# Rdata[Rdata['psfmagerr'] < 0.2]

Bdata= B1data[(B1data['psfmagerr'] < 0.2) &  (B1data['starnum'] >= 5)] - 0.09
Vdata= V1data[(V1data['psfmagerr'] < 0.2) &  (V1data['starnum'] >= 5)] + 0.02
Rdata= R1data[(R1data['psfmagerr'] < 0.2) &  (R1data['starnum'] >= 5)] + 0.21

# ABoffset
#filter cen_wl  mab-mv  Msun_ab Msun_v
#    B	4344	-0.09	5.36	5.45
#    V	5456	0.02	4.80	4.78
#    R	6442	0.21	4.61	4.41


# Distance= 31.17 +/- 0.10 Tully+2009
# extinction 
# E(B-V)host= 0.40+/-0.06 Xiang+2019
# Rv=3.1
# ABhost,AVhost,ARhost
# AVHost= 1.29+/-0.18, Rv = 3.1
# Cardelli 1989

# From SK Lee analysis, E(B-V) = 0.225 Rv=4.05 following Calzetti (2000) Laww
ABhost, AVhost, ARhost = 1.72, 1.29, 0.97  # Xiang
ABhost, AVhost, ARhost = 1.38, 0.91, 0.85  # SED fitting SK Lee

Vmax=Vmaxdate =  57909.86+3.6 # (Xiang+2019)
firstdetect= 2457899.27211
explosiontime=57897.0 # +/-0.3 (Xiang+2019)
ABMW,AVMW,ARMW=0.077,0.058,0.046
Bext,Vext,Rext=ABMW+ABhost, AVMW+AVhost, ARMW+ARhost


plt.rc('text',usetex=True)
plt.rc('font', family='serif')


filtername=['B','V','R']

#plt.errorbar(Bdata['jd']-2400000.5,Bdata['psfmag'],Bdata['psfmagerr'],marker='o',fillstyle='none',color='blue',ls='')

#plt.errorbar(Bdata['jd'][Bdata['obs']=obsname]-2400000.5, Bdata['psfmag'][obsname[i]], Bdata['psfmagerr'][obsname[i]], marker=obsmarker[i],fillstyle='none',color='blue',ls='')


#plt.errorbar(Rdata['jd']-2400000.5,Rdata['psfmag'],Vdata['psfmagerr'],marker='o',fillstyle='none',color='green',ls='')

#plt.errorbar(Rdata['jd']-2400000.5,Rdata['psfmag'],Rdata['psfmagerr'],marker='o',fillstyle='none',color='red',ls='')

#plt.legend()


plt.plot(-99,-99,label='B + 0.5',color='blue')
plt.plot(-99,-99,label='V',color='green')
plt.plot(-99,-99,label='R - 0.5',color='red')
plt.legend(loc=2,frameon=False)

obsnames=['30inch', 'DOAO', 'LOAO', 'Maidanak', 'SOAO','CCA250']
obsmarkers=['+','s','o','x','^','D']
for i  in range(len(obsnames)):
    plt.scatter([],[],label=obsnames[i],marker=obsmarkers[i],color='k')

plt.legend(loc=1,frameon=False)

plt.xlabel('MJD',fontsize=14)
plt.ylabel('Magnitude (Vega)',fontsize=14)
plt.title('Light Curve of SN 2017ein',fontsize=16)
plt.xlim(57890,57990)
#plt.legend()
plt.ylim(20,13)



for i in range(len(obsnames)) :
    plt.errorbar(Bdata[Bdata['obs']==obsnames[i]]['jd']-2400000.5, Bdata[Bdata['obs']==obsnames[i]]['psfmag']+0.5, Bdata[Bdata['obs']==obsnames[i]]['psfmagerr'], \
        marker=obsmarkers[i],fillstyle='none',color='blue',ls='',legend=None)

for i in range(len(obsnames)) :
    plt.errorbar(Vdata[Vdata['obs']==obsnames[i]]['jd']-2400000.5, Vdata[Vdata['obs']==obsnames[i]]['psfmag'], Vdata[Vdata['obs']==obsnames[i]]['psfmagerr'], \
        marker=obsmarkers[i],fillstyle='none',color='green',ls='',legend=None)

for i in range(len(obsnames)) :
    plt.errorbar(Rdata[Rdata['obs']==obsnames[i]]['jd']-2400000.5, Rdata[Rdata['obs']==obsnames[i]]['psfmag']-0.5, Rdata[Rdata['obs']==obsnames[i]]['psfmagerr'], \
        marker=obsmarkers[i],fillstyle='none',color='red',ls='',legend=None)


plt.savefig('SN2017ein_LC_BVR.eps',format='eps',dpi=1200)


plt.close()

Rdata.rename_column('filename','#filename')
ascii.write(Rdata,'R-phot-final.dat',overwrite=True)
Vdata.rename_column('filename','#filename')
ascii.write(Vdata,'V-phot-final.dat',overwrite=True)
Bdata.rename_column('filename','#filename')
ascii.write(Bdata,'B-phot-final.dat',overwrite=True)



######## Color - Color
Vmax=Vmaxdate =  57909.86+3.6 # (Xiang+2019)
firstdetect= 2457899.27211
# ABoffset

ABMW,AVMW,ARMW=0.077,0.058,0.046
#ABhost,AVhost,ARhost

# Distance= 31.17 +/- 0.10 Tully+2009
# extinction 
# E(B-V)host= 0.40+/-0.06 Xiang+2019
# Rv=3.1
# ABhost,AVhost,ARhost
# AVHost= 1.29+/-0.18, Rv = 3.1
# Cardelli 1989
ABhost,AVhost,ARhost = 1.72,1.29,0.97

Bext,Vext,Rext=ABMW+ABhost, AVMW+AVhost, ARMW+ARhost

BVdata = ascii.read('BV-color-final.dat')
VRdata = ascii.read('VR-color-final.dat')
BRdata = ascii.read('BR-color-final.dat')


# B-V

plt.subplots(3,1,sharex=True)
plt.subplots_adjust(hspace=0.2)
plt.subplot(311,xlim=(-5,70))
bv=(BVdata['psfmag_1']-Bext)-(BVdata['psfmag_2']-Vext)
bverr=np.sqrt(BVdata['psfmagerr_1']**2 + BVdata['psfmagerr_2']**2)
plt.errorbar(BVdata['jd_1'][bverr<0.2]-firstdetect ,bv[bverr<0.2],bverr[bverr<0.2], marker='o',fillstyle='none',color='k',ls='')
plt.ylabel('B-V',fontsize=12)
plt.vlines(Vmax-firstdetect+2400000.5,ymin=-2,ymax=2,linestyles='--')
#plt.xlim(-5,70)
plt.ylim(-0.5,1.5)
plt.title ('Color evolution')
# V-R
plt.subplot(312,xlim=(-5,70))
vr=(VRdata['psfmag_1']-Vext)-(VRdata['psfmag_2']-Rext)
vrerr=np.sqrt(VRdata['psfmagerr_1']**2 + VRdata['psfmagerr_2']**2)
plt.errorbar(VRdata['jd_1'][vrerr<0.2]-firstdetect ,vr[vrerr<0.2],vrerr[vrerr<0.2], marker='o',fillstyle='none',color='k',ls='')
plt.ylabel('V-R',fontsize=12)
plt.vlines(Vmax-firstdetect+2400000.5,ymin=-2,ymax=2,linestyles='--')
#plt.xlim(-5,70)
plt.ylim(-0.5,1.5)

# B-R 
plt.subplot(313,xlim=(-5,70))
br=(BRdata['psfmag_1']-Bext)-(BRdata['psfmag_2']-Rext)
brerr=np.sqrt(BRdata['psfmagerr_1']**2 +BRdata['psfmagerr_2']**2)
plt.errorbar(BRdata['jd_1'][brerr<0.2]-firstdetect ,br[brerr<0.2],brerr[brerr<0.2], marker='o',fillstyle='none',color='k',ls='')
plt.ylabel('B-R',fontsize=12)
plt.vlines(Vmax-firstdetect+2400000.5,ymin=-2,ymax=2,linestyles='--')
#plt.xlim(-5,70)
plt.ylim(-0.5,1.5)


plt.xlabel('Day since first detection',fontsize=14)

plt.savefig('SN2017ein-color.eps',dpi=1200)

plt.close()



######## discovery fits image plot
# image plot for SN 2017 discovery 
#

import astropy.io.fits as fits
import aplpy 
import matplotlib.pyplot as plt 
from astropy.visualization import MinMaxInterval, SqrtStretch, ImageNormalize                     
from astropy.visualization import ZScaleInterval, LinearStretch         
from astropy.wcs import WCS
from regions import PixCoord, CirclePixelRegion
from matplotlib.patches import Circle



# image,hdr=fits.getdata('gx20190117_218_mos.fit',header=True) 
# wcs=WCS(hdr)

 

# norm1 = ImageNormalize(image, interval=MinMaxInterval(), stretch=SqrtStretch()) 
# norm2 = ImageNormalize(image, interval=ZScaleInterval(), stretch=LinearStretch()) 


# fig = plt.figure() 

# ax = plt.subplot(projection=wcs)

# im = ax.imshow(image, origin='lower', norm=norm2) 

# fig.colorbar(im) 


pathto='/data7/changsu/sngal/NGC3938/MAO/FLI/R/'

inputfits1 = pathto+'Calibrated-Maidanak_FLI-NGC3938-20170520-181903-R-60_3_com_gregister.fits'
inputfits2 = pathto+'Calibrated-Maidanak_FLI-NGC3938-20170523-182640-R-60_3_com_gregister.fits'
inputfits3 = pathto+'Calibrated-Maidanak_FLI-NGC3938-20170525-183150-R-60_3_com_gregister.fits'
inputfits4 = pathto+'Calibrated-Maidanak_FLI-NGC3938-20170527-180840-R-60_4_com_gregister.fits'


indata1,hdr1= fits.getdata(inputfits1,header=True)
indata2,hdr2= fits.getdata(inputfits2,header=True)
indata3,hdr3= fits.getdata(inputfits3,header=True)
indata4,hdr4= fits.getdata(inputfits4,header=True)





xlim1,xlim2 = 1300, 1600
ylim1,ylim2 = 1300, 1900

wcs=WCS(hdr1)
plt.subplots_adjust(hspace=0.01,wspace=0.01)
plt.title("Emergency of SN 2017ein",loc='center',fontsize=14)

ax1=plt.subplot(141)
#ax1=plt.subplot(141,projection=wcs)
norm2 = ImageNormalize(indata1, interval=ZScaleInterval(), stretch=LinearStretch()) 
ax1.imshow(indata1,cmap='gist_heat',norm=norm2,origin='lower')
#plt.axis('off')
#
ax1.set_xlim(xlim1,xlim2)
ax1.set_ylim(ylim1,ylim2)
circ = Circle((1472, 1593), 15, color='white', fill=None, linewidth='1')
ax1.add_patch(circ)
ax1.set_axis_off()
ax1.text(1350,1850, '2017-05-20.76',fontsize=10, weight='bold',color='white')

# plt.arrow(1450,1590,1455,1595,head_width=0.5,head_length=0.5,fc='white',ec='white')


wcs=WCS(hdr2)
ax2=plt.subplot(142)
#ax2=plt.subplot(142,projection=wcs)
norm2 = ImageNormalize(indata2, interval=ZScaleInterval(), stretch=LinearStretch()) 
ax2.imshow(indata2,cmap='gist_heat',norm=norm2,origin='lower')
ax2.set_xlim(xlim1,xlim2)
ax2.set_ylim(ylim1,ylim2)
circ = Circle((1472, 1593), 15, color='white', fill=None, linewidth='1')
ax2.add_patch(circ)
ax2.set_axis_off()
ax2.text(1350,1850, '2017-05-23.77',fontsize=10, weight='bold',color='white')

wcs=WCS(hdr3)
ax3=plt.subplot(143)
#ax2=plt.subplot(142,projection=wcs)
norm2 = ImageNormalize(indata3, interval=ZScaleInterval(), stretch=LinearStretch()) 
ax3.imshow(indata3,cmap='gist_heat',norm=norm2,origin='lower')
ax3.set_xlim(xlim1,xlim2)
ax3.set_ylim(ylim1,ylim2)
circ = Circle((1472, 1593), 15, color='white', fill=None, linewidth='1')
ax3.add_patch(circ)
ax3.set_axis_off()
ax3.text(1350,1850, '2017-05-25.77',fontsize=10, weight='bold',color='white')
ax3.text(1350,1800, 'First Detection',fontsize=10, weight='bold',color='white')

wcs=WCS(hdr4)
ax4=plt.subplot(144)
#ax2=plt.subplot(142,projection=wcs)
norm2 = ImageNormalize(indata4, interval=ZScaleInterval(), stretch=LinearStretch()) 
ax4.imshow(indata4,cmap='gist_heat',norm=norm2,origin='lower')
ax4.set_xlim(xlim1,xlim2)
ax4.set_ylim(ylim1,ylim2)
circ = Circle((1472, 1593), 15, color='white', fill=None, linewidth='1')
ax4.add_patch(circ)
ax4.set_axis_off()
ax4.text(1350,1850, '2017-05-27.75',fontsize=10, weight='bold',color='white')

plt.savefig('SN2017ein_emergency.eps',format='eps',dpi=1200)
'''

fig=plt.figure()
f1=aplpy.FITSFigure(inputfits1,figure=fig,subplot=[0.05,0.05,0.35,0.8])
f1.show_grayscale()
f2=aplpy.FITSFigure(inputfits2,figure=fig,subplot=[0.25,0.05,0.35,0.8])
f2.show_grayscale()

f2.hide_yaxis_label()
f2.hide_ytick_labels()

fig.canvas.draw()

'''


# Bolometric light curve plot 

BRdata = ascii.read('../BR-color-final.dat')
BRdata.sort('jd_1') 
br=(BRdata['psfmag_1']-Bext)-(BRdata['psfmag_2']-Rext)
brAB=br[brerr<0.2]-0.3
BOLbr = -0.029 -0.302*(brAB)-0.224*(brAB**2)

#plt.plot (BRdata['jd_1'][brerr<0.2]-firstdetect,brAB,'bo') 
MBbol = BRdata['psfmag_1'][brerr<0.2]-Bext+BOLbr 
plt.plot (BRdata['jd_1'][brerr<0.2]-firstdetect,MBbol-31.17,'bo',label='M Bol')
plt.title('Bolometric Light Curve')
plt.xlabel('Date from Fisrt Detection (MJD)')
plt.xlim (-5,70)
plt.ylim(-14,-18)
plt.ylabel('$M_{Bol}$')

spl = UnivariateSpline(BRdata['jd_1'][brerr<0.2], MBbol ,k=1,s=0.0)#,bbox=[BRdata['jd_1'][0],BRdata['jd_1'][-1]])
xs = np.arange(np.min(BRdata['jd_1'])-5, np.max(BRdata['jd_1']+5), 0.001)
# plt.close()

#plt.plot(Rdata['jd']-2457912.23911262,Rdata['psfmag'],'ro') 
plt.plot(xs-firstdetect,spl(xs)-31.17,'k')
#plt.ylim(-13,-18)
#plt.show()

# np.where( spl(xs)==np.min(spl(xs)))
# xs[np.where( spl(xs)==np.min(spl(xs)))]
bolpeak=2457912.23911262

boltemp=ascii.read('../Lyman2016_Bol_LC_template.dat') 
phase=boltemp['col1']
nphase = np.array(phase)

IbL,IbLerr = np.array(boltemp['col2']),np.array(boltemp['col3'])
IcL,IcLerr = np.array(boltemp['col4']),np.array(boltemp['col5'])

# IcL
# Bol_L= L0 *10 ^ (-0.4 * Mbol)
# L0 = 3.828×10^33 erg/s 
Mbol_IcL = np.log10(10**IcL / (3.828*10**33))/-0.4 + 5.5
Mbol_IbL = np.log10(10**IbL / (3.828*10**33))/-0.4 + 5.0
plt.plot(phase + bolpeak-firstdetect, Mbol_IbL, 'r^',label='Ib')
plt.plot(phase + bolpeak-firstdetect, Mbol_IcL, 'gs',label='Ic')
plt.legend()
plt.close()


### physical explosion parameter estimation
'''
Mpeak_BOL  =  np.min(MBbol)-31.17 

M_Ni= 10**(-0.415 * (-17.14)-8.184) # = 0.085

Mej = 1/2 * (b*c/k) * (tm**2)*vsc = 2.0357 M_sun
b=13.7
c=2.99792458*10**10 # cm/s
k=0.07 # opacity
tp=14.139 #explosion to Bol peak
tm=tp #
vph=9300 #km/s Xiang+2019
vsc= vph

Ek = 3/10 * Mej * vph**2 =1.056 * 10^51 erg
'''
### Valenti+2008 Bolometric LC function
from scipy.integrate import quad
from numpy import exp

Msun = 1.98847 *  10**33 # g
M_ej = 2.0357 * Msun
Ek = 1.056 * 10**51 # erg
M_Ni = 0.085 * Msun

k = 0.07 # optical opacity
b = 13.8 # integran constant
c = 2.99792458 * 10**10 # cm/s,  speed of light

tau_m = np.sqrt(k / b * c) * ((10*M_ej) / (3* Ek))**(1/4)
t1=np.linspace(0,60,6001)

tau_Ni, tau_Co = 8.80 * 86400, 111.3 * 86400 # decay time (days)
y = tau_m / (2*tau_Ni)

s = tau_m * (tau_Co - tau_Ni) / (2 * tau_Co * tau_Ni)

e_Ni = 3.90 * 10**10 # erg/s/g,   
e_Co = 6.78 * 10**9 # erg/s/g, 

# integral, Az = 2 * z *np.exp(-2*z*y + z**2)
def integrandAz(z): return  2 * z *np.exp(-2*z*y + z**2)
#integral_Az = quad(integrandAz, 0, x, args=(y))
# integral, Bz = 2 * z *np.exp(-2*z*y + 2*z*s+ z**2)
def integrandBz(z) : return 2 * z *np.exp(-2*z*y + 2*z*s+ z**2)
#integral_Bz = quad(integrandBz, 0, x, args=(y,s))

Lph1=[]
for i in range(len(t1)) :
    t=t1[i] 
    x = t / tau_m
    def integrandAz(z): 
        return  2 * z *np.exp(-2*z*y + z**2)
    def integrandBz(z) : 
        return 2 * z *np.exp(-2*z*y + 2*z*s+ z**2)

    Lph = M_Ni * np.exp(-x**2)  * ( (e_Ni - e_Co)* quad(integrandAz, 0, x)[0] + e_Co * quad(integrandBz, 0, x)[0]   ) 
    Lph1.append(Lph) 
    print (t,np.log10(Lph))

Lph1=[]
for t in t1[:-1]:    
    
    x = t / tau_m
    def integrandAz(z): 
        return  2 * z *np.exp(-2*z*y + z**2)
    def integrandBz(z) : 
        return 2 * z *np.exp(-2*z*y + 2*z*s+ z**2)

    Lph = M_Ni * np.exp(-x**2)  * ( (e_Ni - e_Co)* quad(integrandAz, 0, x)[0] + e_Co * quad(integrandBz, 0, x)[0]   ) 
    Lph1.append(Lph) 
    print (t,np.log10(Lph))


plt.plot(t1[:1900],Lph1[:1900],'ro')

