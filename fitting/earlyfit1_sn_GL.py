import extinction
import astropy.io.ascii as ascii
import numpy as np
f=extinction.Fitzpatrick99(1.766)
A_lamda=f(np.array([4400,5400,6500,8100]),1.06)
mwdust=[0.123, 0.095, 0.076, 0.055]

bandlist  = ['B', 'V', 'R', 'I']#,'g','r','i']#, 'J', 'H', 'K']

incat1=incat[incat['DET']=='Y'] # use data detected only
Rcat=incat1[incat1['FILTER']=='R']
Vcat=incat1[incat1['FILTER']=='V']
Bcat=incat1[incat1['FILTER']=='B']
Icat=incat1[incat1['FILTER']=='I']

B1cat=Bcat[(Bcat['MJD']>Tmax-10) & (Bcat['MJD']<Tmax+20)]
x_data    = B1cat['MJD']
y_data    = B1cat['MAG']
yerr_data = B1cat['MAGERR']


Rv=1.766

def EarlyFit1(incat_name='hgDcSN2019ein.LCfinal.txt', fitcat_name='SN2019ein-PolyFitRes.dat'):
    """
    multi band, calculate T0
    Minimize the sum of the square of the residual.
    """
    case = 1
    # Polynomial or Spline fitting result
    fitcat = ascii.read(fitcat_name)
    # TBmax & Distance
    TBmax  = fitcat['Tmax'][0] ; TBmaxErr =  fitcat['TmaxErr'][0]
    D_res  = GetDist(fitcat_name,Rv)
    DM     = D_res['DM']
    D      = D_res['D']
    # Light curve data (MW & Host reddening corrected)
    incat    = ascii.read(incat_name)
    earlycat = incat[incat['MJD'] < incat['MJD'][0] + 8.3]
    fitidx   = np.where((earlycat['B'] != -99) & (earlycat['V'] != -99) &
                        (earlycat['R'] != -99) & (earlycat['I'] != -99))[0]
    xdata    = earlycat['MJD'][fitidx]
    ydata1   =   earlycat['B'][fitidx]
    ydata2   =   earlycat['V'][fitidx]
    ydata3   =   earlycat['R'][fitidx]
    ydata4   =   earlycat['I'][fitidx]
    ydataerr1 = earlycat['Berr'][fitidx]
    ydataerr2 = earlycat['Verr'][fitidx]
    ydataerr3 = earlycat['Rerr'][fitidx]
    ydataerr4 = earlycat['Ierr'][fitidx]
    ###
    #https://stackoverflow.com/questions/57486632/optimization-of-variables-to-have-the-best-fit-for-two-curves-at-the-same-time
    def MyFunc(xdata, m0B, m0V, m0R, m0I, aB, aV, aR, aI, t0):
        model1 = m0B - 2.5*aB*np.log10(xdata - t0)
        model2 = m0V - 2.5*aV*np.log10(xdata - t0)
        model3 = m0R - 2.5*aR*np.log10(xdata - t0)
        model4 = m0I - 2.5*aI*np.log10(xdata - t0)
        return np.concatenate((model1, model2, model3, model4))
    import lmfit
    MyModel = lmfit.Model(MyFunc)
    params  = lmfit.Parameters()
    a_init  = 2 # alpha, fireball
    m0_init = 19. # normalize, initial value
    params.add('m0B', m0_init)
    params.add('aB', a_init)
    params.add('m0V', m0_init)
    params.add('aV', a_init)
    params.add('m0R', m0_init)
    params.add('aR', a_init)
    params.add('m0I', m0_init)
    params.add('aI', a_init)
    params.add('t0', 58603) # ealier than first detection time , ex) -1
    MyData  = np.concatenate((ydata1, ydata2, ydata3, ydata4))
    MyError = np.concatenate((ydataerr1, ydataerr2, ydataerr3, ydataerr4))
    o1 = MyModel.fit(MyData, params, xdata=xdata, weights=1./MyError, method='leastsq')
    print(o1.fit_report())
    #residual = ((MyData - MyFunc(xdata, m0B, m0V, m0R, m0I, aB, aV, aR, aI, t0))/MyError)**2
    #ChiSq = residual.sum()
    #Reduced_ChiSq = ChiSq / (len(MyData) - 9)
    #fit_params = o1.params
    m0B, m0V, m0R, m0I = o1.params['m0B'].value, o1.params['m0V'].value, \
                            o1.params['m0R'].value,\
                            o1.params['m0I'].value
    m0BErr, m0VErr, m0RErr, m0IErr = o1.params['m0B'].stderr,\
                                    o1.params['m0V'].stderr, \
                                    o1.params['m0R'].stderr, \
                                    o1.params['m0I'].stderr
    aB, aV, aR, aI = o1.params['aB'].value, \
                    o1.params['aV'].value, \
                    o1.params['aR'].value, \
                    o1.params['aI'].value
    aBErr, aVErr, aRErr, aIErr = o1.params['aB'].stderr, \
                                o1.params['aV'].stderr, \
                                o1.params['aR'].stderr, \
                                o1.params['aI'].stderr
    t0, t0Err = o1.params['t0'].value, \
                o1.params['t0'].stderr
    # Fitting lines (total data)
    #x       = np.linspace(t0, np.max(earlycat['MJD'][fitidx]), num = int(  (np.max(earlycat['MJD'][fitidx]))/(360.*(1./24)*(1./60)))  ) # 360 min int.
    x       = np.linspace(t0, np.max(earlycat['MJD'][fitidx]), num = len(earlycat['MJD'][fitidx])*2 )
    #yB      = m0B - 2.5*aB*np.log10(x - t0)
    #yV      = m0V - 2.5*aV*np.log10(x - t0)
    #yR      = m0R - 2.5*aR*np.log10(x - t0)
    #yI      = m0I - 2.5*aI*np.log10(x - t0)
    bandlist = ['B', 'V', 'R', 'I']
    for band in bandlist :
        m0 = o1.params['m0'+band]
        a  = o1.params['a'+band]
        y  = m0 - 2.5*a*np.log10(x - t0)
        #Save fitting lines
        SPLtbl = Table({'x' : x, 'y' : y}, names=['x', 'y'])
        ascii.write(SPLtbl, output='SN2019ein-SPLXY_case{}_{}.dat'.format(case, band), overwrite=True)
        # Residuals
        allidx = np.where((incat[band] != -99) & (incat['MJD'] < incat['MJD'][0] + 8.3))[0]
        residual = incat[band][allidx] - (m0 - 2.5*a*np.log10(incat['MJD'][allidx] - t0))
        resitbl  = Table({'OBS' : incat['OBS'][allidx],
                        'xr' :incat['MJD'][allidx], 'yr' : residual,
                        'yrerr' : incat[band+'err'][allidx]}, names=['OBS', 'xr', 'yr', 'yrerr'])
        ascii.write(resitbl, output='SN2019ein-SPLXY-Residual_case{}_{}.dat'.format(case, band), overwrite=True)


def GetDist(fitcat_name,Rv):
    """
    fitcat_name='SN2019ein-PolyFitRes.dat'
    DM = m - M
    DMErr = sqrt(mErr**2 + MErr**2)
    """
    fitcat = ascii.read(fitcat_name)
    Bmax    = fitcat['MAGmax'][0]
    BmaxErr = fitcat['MmaxErr'][0]
    AbsBmax    = fitcat['AbsMAGmax'][0]
    AbsBmaxErr = fitcat['AbsMAGmaxErr'][0]
    Host    = GetEBVhost(fitcat_name,Rv)
    Ahost_B = Host['Ahost'][0]
    DM    = Bmax - AbsBmax - Ahost_B
    DMErr = np.sqrt(BmaxErr**2 + AbsBmaxErr**2)
    D     = (10.**((DM+5.)/5))/1.e+06
    DErr  = ((10.**(0.2*(DM+5)))*0.2*np.log(10)*DMErr)/1.e+06
    print('Distance modulus is estimated as {}+/-{} mag'.format(round(DM, 3), round(DMErr, 3)))
    print('Distance is estimated as {}+/-{} Mpc'.format(round(D, 3), round(DErr,3)))
    return {'DM' : DM, 'DMErr' : DMErr, 'D' : D, 'DErr': DErr}

def CalChiSq(fit, x, y, yerr, N, n_free) :
    '''
    fit        : Input array from the fitting. ex. simplePL(t, x, *popt)
    x, y, yerr : Input data set as array.
    N          : Total number of data points. ex. len(y)
    n_free     : The number of parameters that we are fitting.
                 If the model has no free parameter then dof = N.
                 If we are fitting a model with n_free free parameters,
                 we can force the model to exactly agree with n_free data points and dof = N - n_free.
    '''
    #if sum(((fit - y)/yerr)**2) > 1.e-03 :
    if N-n_free <= N :
        if sum(((fit - y)/yerr)**2) == 0 :
            print('fit - y = 0, Rejected')
            RedCSQ = -99
        else :
            RedCSQ = (1.0/(N-n_free))*sum(((fit - y)/yerr)**2)
    else :
        print('N-n_free <= N, this will be rejected.')
        RedCSQ = -99
    return RedCSQ

def GetEBVhost(fitcat_name, Rv) :
    """
    fitcat_name = 'SN2019ein-PolyFitRes.dat'
    Rv=3.1
    Bring the result of polynomial fitting.
    ==== Host galaxy extinction (Philips+99) ==='
    E(B-V)host = (Bmax - Vmax)c - (Bmax - Vmax)0
    (Bmax - Vmax)c = (Bmax -Vmax)obs - E(B-V)gal - K_BmaxVmax
    (Bmax - Vmax)0 = -0.07(+/-0.012) +0.114(+/-0.037)*(dB15-1.1), sigma   = 0.03
    BVRIJHK = 01234567
    """
    import extinction
    fitcat  = ascii.read(fitcat_name)
    MAGmax  = fitcat['MAGmax']; MAGmaxErr = fitcat['MmaxErr']
    Bmax    = MAGmax[0]       ; Vmax      = MAGmax[1]
    BmaxErr = MAGmaxErr[0]    ; VmaxErr   = MAGmaxErr[1]
    dM15    = fitcat['dM15']  ; dM15Err   = fitcat['dM15Err']
    dB15    = dM15[0]         ; dB15Err  = dM15Err[0]
    BVmaxObs      = Bmax - Vmax # Galactic extinction corrected
    # BVmaxObs here is BVmaxobs - (A_B-A_V) in the previous code.
    BVmaxObsErr   = np.sqrt(BmaxErr**2 + VmaxErr**2)
    K_BmaxVmax    = 0. # Assumption for a local object.
    BVmax0        = -0.07 + 0.114*(dB15 - 1.1)
    BVmax0_Sigma  = 0.03
    BVmax0Err     = np.sqrt(BVmax0_Sigma**2 + dB15Err**2)
    BVmaxc        = BVmaxObs - K_BmaxVmax
    EBVhost       = BVmaxc - BVmax0
    EBVhostErr    = np.sqrt(BVmaxObsErr**2 + BVmax0Err**2)
    Rv_gal    = Rv; Rv_HV = 1.55 # +/-0.06 (Wang+09)
    Ahost_V   = EBVhost * Rv_gal
    FilterSet = np.array([4400., 5400, 6500.]#, 8797., 12483.0, 16313.0, 22010.0])
    Ahost     = extinction.fitzpatrick99(FilterSet, Ahost_V, Rv_gal)
    print('*** Host Galaxy Extinction (Philips+99) Results ***')
    print('E(B-V)host = (Bmax - Vmax)c - (Bmax - Vmax)0')
    print('(Bmax - Vmax)c = (Bmax -Vmax)obs - E(B-V)gal - K_BmaxVmax')
    print('(Bmax - Vmax)0 = -0.07(+/-0.012) +0.114(+/-0.037)*(dB15-1.1), sigma = 0.03')
    print('(Bmax - Vmax)0 = {}+/-{} (MW corrected)'.format(round(BVmaxObs, 3), round(BVmaxObsErr,3)))
    print('(Bmax - Vmax)c = {} (K-correction = {})'.format(round(BVmaxc, 3), K_BmaxVmax))
    print('E(B-V)host = {}+/-{}'.format(round(EBVhost,3), round(EBVhostErr, 3)))
    return {'EBVhost': EBVhost, 'EBVhostErr': EBVhostErr, 'Ahost': Ahost}


def SPL(t, t0, a, mag0):
    """
    f = [(t-t0)]**alpha
    t   : time
    t0  : Explosion time
    a   : Power law index.
    mg0 : Magnitude at explosion time t0
    """
    import numpy as np
    mag = mag0 - 2.5*a*np.log10(t - t0)
    return mag

def SPL2(t, t0, mag0):
    """
    f = [(t-t0)]**alpha
    t   : time
    t0  : Explosion time
    a   : Power law index. Fixed to 2.
    mg0 : Magnitude at explosion time t0
    """
    import numpy as np
    mag = mag0 - 2.5*2.0*np.log10(t - t0)
    return mag

#===============================================================================


def EarlyFit2Pre(incat_name='hgDcSN2019ein.LCfinal.txt', fitcat_name='SN2019ein-PolyFitRes.dat', FreeN=True) :
    """
    Perform preliminary Simple power-law fitting for the case 1 (2nd data < t < 8.3d) to decide the reference t0.
    True : power index n
    False : power index 2
    """
    from scipy.optimize import curve_fit
    case   = 2
    incat  = ascii.read(incat_name)
    # Polynomial or Spline fitting result
    fitcat = ascii.read(fitcat_name)
    # TBmax & Distance
    TBmax  = fitcat['Tmax'][0] ; TBmaxErr =  fitcat['TmaxErr'][0]
    D_res  = GetDist(fitcat_name,Rv)
    DM     = D_res['DM']
    D      = D_res['D']
    # Early data fitting with simple power-law (Free power index)
    bandlist     = ['B', 'V', 'R']
    colorlist    = ['royalblue', 'seagreen', 'orange']#, 'orangered']
    if FreeN :
        output_name  = 'SN2019ein-SPLFit2PreResN.dat'
        # alpha is not fixed to 2.
    else :
        output_name  = 'SN2019ein-SPLFit2PreRes2.dat'
        # alpha is fixed to 2.


    f = open(output_name, 'w+')
    f.write('T0 T0Err alpha alphaErr mag0 mag0Err Trise TriseErr\n')
    for cat, color in zip(catlist, colorlist) :
        # Select early data
        earlycat    =  cat[(cat['MJD'] < cat['MJD'][0] + 13) ]
        # Fitting initial paramters
        if FreeN :
            p0          = [58139., 2., 19.] # [T0, alpha, m0]
            # Do Fitting
            popt, pcov  = curve_fit(SPL, earlycat['MJD'], earlycat['MAG'], sigma=earlycat['MAGERR'],
                                    p0=p0, absolute_sigma=True, maxfev=10000, check_finite=True)
            RedChiSq    = CalChiSq(SPL(earlycat['MJD'], *popt),  earlycat['MJD'],
                                    earlycat['MAG'], earlycat['MAGERR'], len(earlycat['MJD']), 2)
            # Fitting results
            t0       = popt[0]
            t0Err    = np.sqrt(np.diag(pcov))[0]
            alpha    = popt[1]
            alphaErr = np.sqrt(np.diag(pcov))[1]
            mag0     = popt[2]
            mag0Err  = np.sqrt(np.diag(pcov))[2]
            dt0      = incat['MJD'][0] - t0
            Trise    = TBmax - t0
            TriseErr = np.sqrt(TBmaxErr**2 + t0Err**2)
            print('!====== Simple Power-law fitting results [{} band] ======!'.format('R'))
            print('(1) Explosion time t0 = {}+/-{}'.format(round(t0,3), round(t0Err,3)))
            print('(1)-1. {} days before the 1st obs. :'.format(round(dt0, 3)))
            print('(2) Power index \u03B1 = {}+/-{}'.format(round(alpha,3), round(alphaErr,3)))
            print('(3) mag0 = {}+/-{}'.format(round(mag0,3), round(mag0Err,3)))
            print('(4) Reduced \u03C7\u00b2 = {}'.format(round(RedChiSq, 3)))
            print('(5) Rise time trise = {}+/-{} days'.format(round(Trise,3), round(TriseErr, 3)))
            f.write('{} {} {} {} {} {} {} {}\n'.format(t0, t0Err, alpha, alphaErr, mag0, mag0Err, Trise, TriseErr))
        else :
            p0          = [58603., 19.]
            # Do Fitting
            popt, pcov  = curve_fit(SPL2, earlycat['MJD'], earlycat[band], sigma=earlycat[band+'err'], p0=p0, absolute_sigma=True, maxfev=10000, check_finite=True)
            RedChiSq    = CalChiSq(SPL2(earlycat['MJD'], *popt),  earlycat['MJD'], earlycat[band], earlycat[band+'err'], len(earlycat['MJD']), 1)
            # Fitting results
            t0       = popt[0]
            t0Err    = np.sqrt(np.diag(pcov))[0]
            mag0     = popt[1]
            mag0Err  = np.sqrt(np.diag(pcov))[1]
            dt0      = incat['MJD'][0] - t0
            Trise    = TBmax - t0
            TriseErr = np.sqrt(TBmaxErr**2 + t0Ersr**2)
            print('!====== Simple Power-law fitting results [{} band] ======!'.format(band))
            print('(1) Explosion time t0 = {}+/-{}'.format(round(t0,3), round(t0Err,3)))
            print('(1)-1. {} days before the 1st obs. :'.format(round(dt0, 3)))
            print('(2) Power index \u03B1 = 2 (fixed)')
            print('(3) mag0 = {}+/-{}'.format(round(mag0,3), round(mag0Err,3)))
            print('(4) Reduced \u03C7\u00b2 = {}'.format(round(RedChiSq, 3)))
            print('(5) Rise time trise = {}+/-{} days'.format(round(Trise,3), round(TriseErr, 3)))
            f.write('{} {} {} {} {} {} {} {}\n'.format(t0, t0Err, 2.0, -99, mag0, mag0Err, Trise, TriseErr))
        if FreeN :
            SPLXYcat_name = 'SN2019ein-SPLXY_case2_{}N.dat'.format(band)
            resicat_name = 'SN2019ein-SPLXY-Residual_case2_{}N.dat'.format(band)
        else :
            SPLXYcat_name = 'SN2019ein-SPLXY_case2_{}2.dat'.format(band)
            resicat_name = 'SN2019ein-SPLXY-Residual_case2_{}2.dat'.format(band)
        # Fitting lines & Residuals
        x = np.linspace(t0, np.max(earlycat['MJD']), num = len(earlycat['MJD']))
        y = SPL(x, *popt)
        xytbl  = Table({'x' : x, 'y' : y}, names=['x', 'y'])
        ascii.write(xytbl, output=SPLXYcat_name, overwrite=True)
        # Residuals earlycat[band]
        residual = earlycat[band] - SPL(earlycat['MJD'], *popt)
        resitbl  = Table({'OBS' : earlycat['OBS'], 'xr' : earlycat['MJD'], 'yr' : residual, 'yrerr' : earlycat[band+'err']}, names=['OBS', 'xr', 'yr', 'yrerr'])
        ascii.write(resitbl, output=resicat_name, overwrite=True)
    f.close()




import matplotlib.pyplot as plt
plt.ion()

plt.errorbar(earlycat['MJD'],earlycat['MAG'],yerr=earlycat['MAGERR']/2,fmt='ro')
SPL(earlycat['MJD'],*popt)
plt.plot(earlycat['MJD'],SPL(earlycat['MJD'],*popt),'--')
plt.show()
plt.ylim(21,15)
from astropy.time import Time
t = Time('2018-01-24T21:57:12', format='isot', scale='utc')
discoverymjd=t.mjd
discoverymjd
plt.vlines(discoverymjd,ymin=15,ymax=21)
plt.title('EARLY LC FIT')


    from lgpy.sn2019ein_fitting import fearly2_kasen
    from lgpy.sn2019ein_fitting import fearly2_rw10
        # Shock heated cooling emission
        td     = np.array(0.00 + 0.01*np.arange(201), dtype=np.float128)
        y_dict = {}
        print('R* = {} Rsol for {} band'.format(round(Rstar[r],3), band))
        # Kasen
        y      = []
        for i in range(len(td)):
            y_dum = fearly2_kasen(td[i], Rstar[r], band)
            y.append(y_dum)
        y_dict.update({Rstar[r] : y})
        # RW10
        y_rw10 = []
        for i in range(len(td)) :
            y_rw10_dum = fearly2_rw10(td[i], Rstar[r], band)
            y_rw10.append(y_rw10_dum)
