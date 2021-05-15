# code for Prof Yoon's model check with SN 2017ein data

import astropy.io.ascii as ascii
import numpy as np 
import matplotlib.pyplot as plt 
import astropy.modeling as modelling


B1data=ascii.read('B-phot-cut.dat')
V1data=ascii.read('V-phot-cut.dat')
R1data=ascii.read('R-phot-cut.dat')

Bdata= B1data[(B1data['psfmagerr'] < 0.2) & (B1data['starnum'] >= 5)]
Vdata= V1data[(V1data['psfmagerr'] < 0.2) & (V1data['starnum'] >= 5)]
Rdata= R1data[(R1data['psfmagerr'] < 0.2) & (R1data['starnum'] >= 5)]


Vmax=Vmaxdate =  57909.86+3.6 # (Xiang+2019)
firstdetect= 2457899.27211
explosiontime=57897.0 # +/-0.3 (Xiang+2019)
ABMW,AVMW,ARMW=0.077,0.058,0.046
#ABhost,AVhost,ARhost

# Distance= 31.17 +/- 0.10 Tully+2009
# extinction 
# E(B-V)host= 0.40+/-0.06 Xiang+2019
# Rv=3.1
# ABhost,AVhost,ARhost
# AVHost= 1.29+/-0.18, Rv = 3.1
# Cardelli 1989

# From SK Lee analysis, E(B-V) = 0.225 Rv=4.05 following Calzetti (2000) Laww
ABhost, AVhost, ARhost = 1.72, 1.29, 0.97  # Xiang
# ABhost, AVhost, ARhost = 1.38, 0.91, 0.71  # SED fitting SK Lee, Rv= 4.05
# ABhost, AVhost, ARhost = 0.93, 0.70, 0.52  # SED fitting SK Lee, Rv= 3.1



Bext,Vext,Rext=ABMW+ABhost, AVMW+AVhost, ARMW+ARhost

BVdata = ascii.read('BV-color-final.dat')
VRdata = ascii.read('VR-color-final.dat')
BRdata = ascii.read('BR-color-final.dat')


'''
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.0fm0.15.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.0fm0.30.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.0fm0.50.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.0fm0.90.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.0fm5.00.tt

/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.5fm0.15.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.5fm0.30.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.5fm0.50.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.5fm0.90.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.5fm5.00.tt

/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.8fm0.15.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.8fm0.30.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.8fm0.50.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.8fm0.90.tt
/data7/changsu/sngal/NGC3938/phot/snlc-yoon/CO3.93E1.8fm5.00.tt
'''
infilelist=['CO3.93E1.0fm0.15.tt',
'CO3.93E1.0fm0.30.tt',
'CO3.93E1.0fm0.50.tt',
'CO3.93E1.0fm0.90.tt',
'CO3.93E1.0fm5.00.tt',
'CO3.93E1.5fm0.15.tt',
'CO3.93E1.5fm0.30.tt',
'CO3.93E1.5fm0.50.tt',
'CO3.93E1.5fm0.90.tt',
'CO3.93E1.5fm5.00.tt',
'CO3.93E1.8fm0.15.tt',
'CO3.93E1.8fm0.30.tt',
'CO3.93E1.8fm0.50.tt',
'CO3.93E1.8fm0.90.tt',
'CO3.93E1.8fm5.00.tt']

infilelist=['HE3.87E1.0fm0.15.tt',
'HE3.87E1.0fm0.30.tt',
'HE3.87E1.0fm0.50.tt',
'HE3.87E1.0fm0.90.tt',
'HE3.87E1.0fm5.00.tt',
'HE3.87E1.5fm0.15.tt',
'HE3.87E1.5fm0.30.tt',
'HE3.87E1.5fm0.50.tt',
'HE3.87E1.5fm0.90.tt',
'HE3.87E1.5fm5.00.tt',
'HE3.87E1.8fm0.15.tt',
'HE3.87E1.8fm0.30.tt',
'HE3.87E1.8fm0.50.tt',
'HE3.87E1.8fm0.90.tt',
'HE3.87E1.8fm5.00.tt']

plt.close()
plt.figure(figsize=(20,10))
pathto='/data7/changsu/sngal/NGC3938/phot/snlc-yoon/'
for infile in infilelist:

    colnames=['time','Tbb','rbb','Teff','Rlast_sc','R_tau2_3','Mbol','MU','MB','MV','MI','MR','Mbolavg','gdepos']
    data=ascii.read(pathto+infile,header_start=78,data_start=79)
    const =0
    #plt.plot(Rdata['jd']-2400000.5-57897.0+const ,Rdata['psfmag']-31.17-Rext,'s',fillstyle='none')
    plt.plot(data['time'],data['MB'],'o',label=infile)
    #plt.plot(data['time'],data['MV'],'o',label=infile)
    #plt.plot(data['time'],data['MR'],'o',label=infile)

plt.plot(Bdata['jd']-2400000.5-57897.0+const ,Bdata['psfmag']-31.17-Bext,'ks',fillstyle='none',label='SN 2017ein')
#plt.plot(Vdata['jd']-2400000.5-57897.0+const ,Vdata['psfmag']-31.17-Vext,'ks',fillstyle='none',label='SN 2017ein')
#plt.plot(Rdata['jd']-2400000.5-57897.0+const ,Rdata['psfmag']-31.17-Rext,'ks',fillstyle='none',label='SN 2017ein')

plt.ylim(-12,-18)
plt.xlim=(-5,70)
#plt.title(infile)   
plt.legend(loc=4)
plt.ylabel('Abs Mag (MB)')
plt.xlabel('Time since explosion (day)')
plt.title('HE core MB light curve')


plt.savefig('fit-yoon-model-HE-core-MB.png')
#plt.plot(xnew,ynew,'D')
plt.close()

#############################   COLOR   #####################################
plt.close()
plt.figure(figsize=(20,10))
for infile in infilelist:

    colnames=['time','Tbb','rbb','Teff','Rlast_sc','R_tau2_3','Mbol','MU','MB','MV','MI','MR','Mbolavg','gdepos']
    data=ascii.read(pathto+infile,header_start=78,data_start=79)
    const =0
    #plt.plot(Rdata['jd']-2400000.5-57897.0+const ,Rdata['psfmag']-31.17-Rext,'s',fillstyle='none')
    plt.plot(data['time'],data['MB']-data['MR'],'o',label=infile)
    #plt.plot(data['time'],data['MV'],'o',label=infile)
    #plt.plot(data['time'],data['MR'],'o',label=infile)

#plt.plot(Bdata['jd']-2400000.5-57897.0+const ,Bdata['psfmag']-31.17-Bext,'ks',fillstyle='none',label='SN 2017ein')
#plt.plot(Vdata['jd']-2400000.5-57897.0+const ,Vdata['psfmag']-31.17-Vext,'ks',fillstyle='none',label='SN 2017ein')
#plt.plot(Rdata['jd']-2400000.5-57897.0+const ,Rdata['psfmag']-31.17-Rext,'ks',fillstyle='none',label='SN 2017ein')

plt.ylim(-0.5,2.5)
plt.xlim=(-5,70)
#plt.title(infile)   
plt.legend(loc=4)
plt.ylabel('B-R',fontsize=12)
plt.xlabel('Time since explosion (day)',fontsize=12)
plt.title('B-R value for CO core',fontsize=14)


plt.savefig('fit-Yoon-model-CO-core-B-R-color.png')
#plt.plot(xnew,ynew,'D')
plt.close()



# B-V
plt.subplots_adjust(hspace=0)
plt.subplot(311)
bv=(BVdata['psfmag_1']-Bext)-(BVdata['psfmag_2']-Vext)
bverr=np.sqrt(BVdata['psfmagerr_1']**2 + BVdata['psfmagerr_2']**2)
plt.errorbar(BVdata['jd_1'][bverr<0.2]-firstdetect ,bv[bverr<0.2],bverr[bverr<0.2], marker='o',fillstyle='none',color='k',ls='')
plt.ylabel('B-V',fontsize=12)
plt.vlines(Vmax-firstdetect+2400000.5,ymin=-2,ymax=2,linestyles='--')
plt.xlim(-5,70)
plt.ylim(-0.5,1.5)
plt.title ('Color evolution')
plt.set
# V-R
plt.subplot(312)
vr=(VRdata['psfmag_1']-Vext)-(VRdata['psfmag_2']-Rext)
vrerr=np.sqrt(VRdata['psfmagerr_1']**2 + VRdata['psfmagerr_2']**2)
plt.errorbar(VRdata['jd_1'][vrerr<0.2]-firstdetect ,vr[vrerr<0.2],vrerr[vrerr<0.2], marker='o',fillstyle='none',color='k',ls='')
plt.ylabel('V-R',fontsize=12)
plt.vlines(Vmax-firstdetect+2400000.5,ymin=-2,ymax=2,linestyles='--')
plt.xlim(-5,70)
plt.ylim(-0.5,1.5)

# B-R 
plt.subplot(313)
br=(BRdata['psfmag_1']-Bext)-(BRdata['psfmag_2']-Rext)
brerr=np.sqrt(BRdata['psfmagerr_1']**2 +BRdata['psfmagerr_2']**2)
plt.errorbar(BRdata['jd_1'][brerr<0.2]-firstdetect ,br[brerr<0.2],brerr[brerr<0.2], marker='o',fillstyle='none',color='k',ls='')
plt.ylabel('B-R',fontsize=12)
plt.vlines(Vmax-firstdetect+2400000.5,ymin=-2,ymax=2,linestyles='--')
plt.xlim(-5,70)
plt.ylim(-0.5,1.5)


plt.xlabel('Day since first detection',fontsize=14)

plt.savefig('SN2017ein-color.eps',dpi=1200)

plt.close()




### modelling part 
data=ascii.read(pathto+infile,header_start=78,data_start=79)


### interpolation of SN 2017ein data to function
#from scipy.interpolate import CubicSpline
#cs = CubicSpline(Rdata['jd'], Rdata['psfmag'])
#xs = np.arange(np.min(Rdata['jd'])-5, np.max(Rdata['jd']+5), 0.001)

#from scipy import interpolate
#from scipy.interpolate import BSpline
#from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline
from scipy.interpolate import UnivariateSpline
spl = UnivariateSpline(Rdata['jd'], Rdata['psfmag'],k=1,s=0.2,bbox=[Rdata['jd'][0],Rdata['jd'][-1]])
xs = np.arange(np.min(Rdata['jd'])-3, np.max(Rdata['jd']+3), 0.001)
plt.close()

plt.plot(Rdata['jd'],Rdata['psfmag'],'ro') 
plt.plot(xs,spl(xs),'k')
plt.ylim(14,20)
plt.show()


# fitting model data to polynomial 

jd , mag, magerr = Rdata['jd'], Rdata['psfmag'], Rdata['psfmagerr']

p = np.polyfit(jd,mag,7)
p1=  np.poly1d(p)

plt.plot(jd,p1(jd),'k--')


