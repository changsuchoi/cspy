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
# single filter, for SN 2018kp, used for 'R' band only

def EarlyFit2Pre(name='SN2017ein',incat_name='SN2017ein_B.txt',TBmax=58000, TBmaxErr=0.2, DM=30., D=0.5 FreeN=True) :
    """
    Perform preliminary Simple power-law fitting for the case 1 (2nd data < t < 8.3d) to decide the reference t0.
    True : power index n
    False : power index 2
    """
    incat=ascii.read('SN2017ein_r.dat')

    f=extinction.Fitzpatrick99(2.6)
    Av=0.884
    A_lamda=f(np.array([4400,5400,6500,8100]),Av)
    AR_host = 0.65
    AR_mw = 0.046
    incat['MAG']=incat['MAG'] - AR_host - AR_mw
    incat['UL5']=incat['UL5'] - AR_host - AR_mw

    nodet=incat[incat['DET']=='N']
    incat0=incat[incat['DET']=='Y']

    from scipy.optimize import curve_fit
    case   = 2
    incat  = ascii.read(incat_name)
    # Polynomial or Spline fitting result
    # fitcat = ascii.read(fitcat_name)
    # TBmax & Distance
    #TBmax  = fitcat['Tmax'][0] ; TBmaxErr =  fitcat['TmaxErr'][0]
    # D_res  = GetDist(fitcat_name,Rv)
    #DM     = D_res['DM'] #distance measurement
    #D      = D_res['D']  # DM error
    # Early data fitting with simple power-law (Free power index)
    bandlist     = ['B', 'V', 'R']
    colorlist    = ['royalblue', 'seagreen', 'orange']#, 'orangered']
    if FreeN :
        output_name  = name+'-SPLFit2PreResN.dat'
        # alpha is not fixed to 2.
    else :
        output_name  = name+'-SPLFit2PreRes2.dat'
        # alpha is fixed to 2.


    f = open(output_name, 'w+')
    f.write('T0 T0Err alpha alphaErr mag0 mag0Err Trise TriseErr\n')
    #for cat, color in zip(catlist, colorlist) :
        # Select early data
    earlycat    =  incat0[(incat0[incat0['MJD'] < incat0['MJD'][0] + 10) ]
    # Fitting initial paramters
    if FreeN :
        p0          = [57898., 2., 19.] # [T0, alpha, m0]
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
        print('!====== Simple Power-law fitting results [{} band] ======!'.format('B'))
        print('(1) Explosion time t0 = {}+/-{}'.format(round(t0,3), round(t0Err,3)))
        print('(1)-1. {} days before the 1st obs. :'.format(round(dt0, 3)))
        print('(2) Power index \u03B1 = {}+/-{}'.format(round(alpha,3), round(alphaErr,3)))
        print('(3) mag0 = {}+/-{}'.format(round(mag0,3), round(mag0Err,3)))
        print('(4) Reduced \u03C7\u00b2 = {}'.format(round(RedChiSq, 3)))
        print('(5) Rise time trise = {}+/-{} days'.format(round(Trise,3), round(TriseErr, 3)))
        f.write('{} {} {} {} {} {} {} {}\n'.format(t0, t0Err, alpha, alphaErr, mag0, mag0Err, Trise, TriseErr))
    else :
        p0          = [57898., 19.]
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
        SPLXYcat_name = name+'-SPLXY_case2_{}N.dat'.format(band)
        resicat_name = name+'-SPLXY-Residual_case2_{}N.dat'.format(band)
    else :
        SPLXYcat_name = name+'-SPLXY_case2_{}2.dat'.format(band)
        resicat_name = name+'-SPLXY-Residual_case2_{}2.dat'.format(band)
    # Fitting lines & Residuals
    x = np.linspace(t0, np.max(earlycat['MJD']), num = len(earlycat['MJD']))
    y = SPL(x, *popt)
    xytbl  = Table({'x' : x, 'y' : y}, names=['x', 'y'])
    ascii.write(xytbl, output=SPLXYcat_name, overwrite=True)
    # Residuals earlycat[band]
    residual = earlycat[band] - SPL(earlycat['MJD'], *popt)
    resitbl  = Table({'OBS' : earlycat['OBS'], 'xr' : earlycat['MJD'], 'yr' : residual, 'yrerr' : earlycat[band+'ERR']}, names=['OBS', 'xr', 'yr', 'yrerr'])
    ascii.write(resitbl, output=resicat_name, overwrite=True)
    f.close()




#===============================================================================
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

def single_powerlaw(t, t0, a, mg0):
    """
    t   :
    t0  :
    a   :
    mg0 :
    """
    import numpy as np
    mg = mg0 + 2.5*a*np.log10(t - t0)
    return mg

def fit_single(t, mg, mge, initial, maxfev):
    """
    t       : time relative to specific date.
    mg      : magnitude
    mge     : magnitude error
    initial : list of initial values of t0, a, mg0, enter like [1.5, 2., 16.5]
    """
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(single_powerlaw, t, mg, sigma=mge,  p0=initial, absolute_sigma=True, maxfev=maxfev, check_finite=True)

def planck(wave, temp):
    import numpy as np
    #if len(wave) > 1 :
    #    print('Syntax - bbflux = planck( wave, temp)')
    #    return 0
    #if len(temp) != 1 :
    #    input('Enter a blackbody temperature : ')
    # Gives the blackbody flux (i.e. PI*Intensity) ergs/cm2/s/a
    w = wave/1.e8 # angstroms to cm
    # constants appropriate to cgs units.
    c1 = np.float128(3.7417749e-5)          # =2*!DPI*h*c*c
    c2 = np.float128(1.4387687)             # =h*c/k
    val = c2/w/np.float128(temp)
    bbflux = c1/( (w**5)*(np.exp(val)-1.))
    return bbflux*1.e-8 # convert to ergs cm-2 s-1 A-1

def fearly2_rw10(td, rstar, band) :
    """
    This function produces the shock-heated emssion light curve.
    Written for the SN 2015F work based on the Rabinak & Waxman (2011) model [2015. M. Im].
    Currently, the explosion energy is set to 10^51 erg which is appropriate for Type Ia SN.
    This value may need to be changed for other types of SNe.
    Also, at the end of the program, you will encounter the line
    fearly2_kasen=interpol(mbbflux,xw,6580.)
    Note that the number "6580." is the effective wavelength of the filter in Angstrom.
    In this case, it is set to R-band.
    Please change the number if you want to plot the light curve in different bands. [2018-05-03, added by M. Im]
    Slightly modified at 2018-05-03 to add comments. [2018-05-03, M. Im].
    """
    import numpy as np
    rstar = np.float128(rstar)
    td = np.float64(td)
    # Rsun = 6.955 * 10**10 cm
    r10 = np.float128(rstar*6.955) # rstar in Rsun unit
    r13 = np.float128(rstar*6.955e-3)
    # rstar is the radius of progenitor or companion in solar radius unit.
    # In K10 model, the emission is from an interaction btw the companion and the ejected materials so r13 is that of companion.
    # In RW11 model, the emission is from the progenitor itself. r13 is that of progenitor
    # Progenitor radius in R/10^10 cm
    Mc = np.float128(1.0/1.40)
    # Mass in chandrasekhar mass unit
    # eff  = 1.0
	# efficiency of conversion of mass to light
    Msun = 1.988e33 # g
    c    = 2.9979e10 # cm/s
    # mc2  = np.log10(Msun) + 2.*np.log10(c)
    # Energy in Msun*c**2
    eff = 1.0
    dm = np.log10(Msun) + 2. * np.log10(2.9979e10) - 51. # Difference btw E51 and Msun*c^2
    loge51 = np.log10(eff*Mc) + np.log10(Msun) + 2.*np.log10(c) - 51. - dm # explosion energy in log unit (10**51 erg)
    #ke  = np.float128(0.2) # Opacity in K10, ke = 0.2 cm**2/g -> k02 = 1, appropriate for e- scattering in fully ionized A/Z=2 elements
    k02 = np.float128(1.0) # Opacity in k/0.2cm**2/g
    #fp = 0.05 # Form factor 0.031 - 0.13 (RW11)
    fp = 0.05 # form factor 0.031 - 0.13
    #v9 = np.float128(1.)   # velocity in 10**9 cm/s unit
    logLt = 40. + np.log10(1.2) + np.log10(r10) + 0.85*loge51 - 0.69*np.log10(Mc) -0.85*np.log10(k02) -0.16*np.log10(fp) - 0.35*np.log10(td) # total luminosity : early UV/optical luminosity emitted from ejecta diffusion including radioactive + cooling emission
    logTeff = np.log10(4105.) + 0.25*np.log10(r10) + 0.016*loge51 + 0.03*np.log10(Mc) - 0.27*np.log10(k02) - 0.022*np.log10(fp) -0.47*np.log10(td)  # Effective temperature of the early emission
    c = np.float(2.9979e8) # speed of light in m/s
    # wavelengh in angstrom, 1A = 10e-10 m
    sigb = np.float128(5.6704e-5) # Stephan-Boltzman constant in CGS [egs*cm**-2 * s**-1 * ]
    logFbb = np.log10(sigb) + 4.*logTeff # Blackbody flux of early emission
    d10pc  = np.float128(10. * 3.0856e18) # 10pc in cm to convert Luminosity to Flux
    fluxm    = (logLt - 2.*np.log10(d10pc) - np.log10(4.*np.pi)) - logFbb # Flux of shock heated emission = Total flux - blackbody component of early emission
    #rph24pi = logLt - logFbb
    abzero = np.float128(-48.600) # in erg*s**-1*cm**-2*Hz**-1

    teff = 10.**(logTeff)
    xw   = 1000. + 80.*np.arange(100)
    xw = np.array(xw, dtype='float128')

    #from lgpy.sn2019ein_fitting import planck
    bbflux = planck(xw, teff)
    ff   = -2.5 * (2.*np.log10(xw) -10. -np.log10(c) ) # Angstrom to Hertz conversion factor
    mbbflux = -2.5*np.log10(bbflux) + ff + abzero -2.5*fluxm # convert flux to AB magnitude (Find ref in wiki).
    x_data  = xw
    y_data  = mbbflux
    from scipy.interpolate import UnivariateSpline
    spl     = UnivariateSpline(x_data, y_data, s=0.2, k=5)
    # From 1000A to 8920A, we fit BB spectrum to obtain AB magnitude in specific band from continuous value. ex) What is magnitude in 6580A? 6580 is not generated in the array wx so we fit mbbflux!
    # The last input parameter in the above function (e.g., 4770., and 6580.) are the effective wavelength of the filter in Angstrom. Please modify the value, depending on which filter you use.
    # 4770  : B band, 6580. : R band
    x_array = np.linspace(np.min(x_data), np.max(x_data), num= int((np.max(x_data) -np.min(x_data))/80.))
    if band == 'U':
        fearly2 = np.float128(spl(3656.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'B' :
        #print('Band : '+band+', eff w = 4770.')
        #fearly2 = np.float128(spl(4770)) # LSGT
        fearly2 = np.float128(spl(4353.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'V' :
        fearly2 = np.float128(spl(5477.))
    elif band == 'R' :
        #print('Band : '+band+', eff w = 6580.')
        #fearly2 = np.float128(spl(6580)) # LSGT
        fearly2 = np.float128(spl(6349.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'I' :
        #fearly2 = np.float128(spl(8175.6))
        fearly2 = np.float128(spl(8797.))
    elif band == 'g' :
        fearly2 = np.float128(spl(4770.))
    elif band == 'r' :
        fearly2 = np.float128(spl(6231.))
    return fearly2

def fearly2_kasen(td, rstar, band):
    """
    This function produces the shock-heated emssion light curve. Written for the SN 2015F work based on the Kasen (2010) model [2015. M. Im].
    Currently, the explosion energy is set to 10^51 erg which is appropriate for Type Ia SN.
    This value may need to be changed for other types of SNe.
    Also, at the end of the program, you will encounter the line
    fearly2_kasen=interpol(mbbflux,xw,6580.)
    Note that the number "6580." is the effective wavelength of the filter in Angstrom.
    In this case, it is set to R-band. Please change the number if you want to plot the light curve in different bands. [2018-05-03, added by M. Im]
    Slightly modified at 2018-05-03 to add comments. [2018-05-03, M. Im].

    * Size and Mass relation (Kasen 2010)
    1 - 3 Msun, MS          : R* = 1 - 3 x 10**11cm
    5 - 6 Msun, MS subgiant : R* = 5 x 10**11cm
    1 - 2 Msun, Red giant   : R* ~ 10**13cm
    """
    import numpy as np
    # rstar = 1.0
    # td = 0.5
    # band = 'R'
    rstar = np.float128(rstar)
    td = np.float64(td)
    # Rsun = 6.955 * 10**10 cm
    r10 = np.float128(rstar*6.955) # rstar in Rsun unit
    r13 = np.float128(rstar*6.955e-3)
    # rstar is the radius of progenitor or companion in solar radius unit.
    # In K10 model, the emission is from an interaction btw the companion and the ejected materials so r13 is that of companion.
    # In RW11 model, the emission is from the progenitor itself. r13 is that of progenitor
    # Progenitor radius in R/10^10 cm
    Mc = np.float128(1.0/1.40)
    # Mass in chandrasekhar mass unit
    # eff  = 1.0
	# efficiency of conversion of mass to light
    # Msun = 1.988e33 # g
    # c    = 2.9979e10 # cm/s
    # mc2  = np.log10(Msun) + 2.*np.log10(c)
    # Energy in Msun*c**2
    ke  = np.float128(1.0) # Opacity in K10, ke = 0.2 cm**2/g -> k02 = 1, appropriate for e- scattering in fully ionized A/Z=2 elements
    k02 = np.float128(5.0) # Opacity in k/0.2cm**2/g
    #fp = 0.05 # Form factor 0.031 - 0.13 (RW11)
    v9 = np.float128(1.)   # velocity in 10**9 cm/s unit
    logLt = 43. + np.log10(2*r13) + 0.25*np.log10(Mc) + (7./4.)*np.log10(v9) + (-0.75)*np.log10(ke) + (-0.5)*np.log10(td) # total luminosity : early UV/optical luminosity emitted from ejecta diffusion including radioactive + cooling emission
    logTeff = np.log10(2.5) + 4. + 0.25*np.log10(2.*r13) - (35./36.)*np.log10(ke) - (37./72.)*np.log10(td) # Effective temperature of the early emission
    c = np.float(2.9979e8) # speed of light in m/s
    # wavelengh in angstrom, 1A = 10e-10 m
    sigb = np.float128(5.6704e-5) # Stephan-Boltzman constant in CGS [egs*cm**-2 * s**-1 * ]
    logFbb = np.log10(sigb) + 4.*logTeff # Blackbody flux of early emission
    d10pc  = np.float128(10. * 3.0856e18) # 10pc in cm to convert Luminosity to Flux
    fluxm    = (logLt - 2.*np.log10(d10pc) - np.log10(4.*np.pi)) - logFbb # Flux of shock heated emission = Total flux - blackbody component of early emission
    abzero = np.float128(-48.600) # in erg*s**-1*cm**-2*Hz**-1

    teff = 10.**(logTeff)
    xw   = 1000. + 40.*np.arange(200)
    #xw   = 1000. + 80.*np.arange(100)
    xw = np.array(xw, dtype='float128')

    #from lgpy.sn2019ein_fitting import planck
    bbflux = planck(xw, teff)
    ff   = -2.5 * (2.*np.log10(xw) -10. -np.log10(c) ) # Angstrom to Hertz conversion factor
    mbbflux = -2.5*np.log10(bbflux) + ff + abzero -2.5*fluxm # convert flux to AB magnitude (Find ref in wiki).
    x_data  = xw
    y_data  = mbbflux
    from scipy.interpolate import UnivariateSpline
    spl     = UnivariateSpline(x_data, y_data, s=0.2, k=5)
    # From 1000A to 8920A, we fit BB spectrum to obtain AB magnitude in specific band from continuous value. ex) What is magnitude in 6580A? 6580 is not generated in the array wx so we fit mbbflux!
    # The last input parameter in the above function (e.g., 4770., and 6580.) are the effective wavelength of the filter in Angstrom. Please modify the value, depending on which filter you use.
    # 4770  : B band
    # 5477 : V band https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    # 6580. : R band
    x_array = np.linspace(np.min(x_data), np.max(x_data), num= int((np.max(x_data) -np.min(x_data))/80.))
    if band == 'U':
        fearly2 = np.float128(spl(3656.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    if band == 'B' :
        #print('Band : '+band+', eff w = 4770.')
        #fearly2 = np.float128(spl(4770)) # LSGT
        fearly2 = np.float128(spl(4353.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'V' :
        fearly2 = np.float128(spl(5477.))
    elif band == 'R' :
        #print('Band : '+band+', eff w = 6580.')
        #fearly2 = np.float128(spl(6580)) # LSGT
        fearly2 = np.float128(spl(6349.)) # https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    elif band == 'I' :
        #fearly2 = np.float128(spl(8175.6))
        fearly2 = np.float128(spl(8797.))
    elif band == 'g' :
        fearly2 = np.float128(spl(4770.))
    elif band == 'r' :
        fearly2 = np.float128(spl(6231.))
    # I (CTIO/ANDICAM.I_KPNO): lamb_eff : 8175.6
    # J (UKIRT) : 12483.0
    # H (UKIRT) : 16313.0
    # K (UKIRT) : 22010.0
    return fearly2

def fearly3_kasen(td, rstar, band):
    """
    This function produces the shock-heated emssion light curve. Written for the SN 2015F work based on the Kasen (2010) model [2015. M. Im].
    Currently, the explosion energy is set to 10^51 erg which is appropriate for Type Ia SN.
    This value may need to be changed for other types of SNe.
    Also, at the end of the program, you will encounter the line
    fearly3_kasen=interpol(mbbflux,xw,6580.)
    Note that the number "6580." is the effective wavelength of the filter in Angstrom.
    In this case, it is set to R-band. Please change the number if you want to plot the light curve in different bands. [2018-05-03, added by M. Im]
    Slightly modified at 2018-05-03 to add comments. [2018-05-03, M. Im].

    fearly3_kasen provides the size of the companion when you enter magnitude. [2020-08-24, added by G. Lim]

    * Size and Mass relation (Kasen 2010)
    1 - 3 Msun, MS          : R* = 1 - 3 x 10**11cm
    5 - 6 Msun, MS subgiant : R* = 5 x 10**11cm
    1 - 2 Msun, Red giant   : R* ~ 10**13cm
    """
    import numpy as np
    # rstar = 1.0
    # td = 0.5
    # band = 'R'
    rstar = np.float128(rstar)
    td    = np.float64(td)
    # Rsun = 6.955 * 10**10 cm
    r10 = np.float128(rstar*6.955) # rstar in Rsun unit
    r13 = np.float128(rstar*6.955e-3)
    # rstar is the radius of progenitor or companion in solar radius unit.
    # In K10 model, the emission is from an interaction btw the companion and the ejected materials so r13 is that of companion.
    # In RW11 model, the emission is from the progenitor itself. r13 is that of progenitor
    # Progenitor radius in R/10^10 cm
    Mc = np.float128(1.0/1.40)
    # Mass in chandrasekhar mass unit
    # eff  = 1.0
	# efficiency of conversion of mass to light
    # Msun = 1.988e33 # g
    # c    = 2.9979e10 # cm/s
    # mc2  = np.log10(Msun) + 2.*np.log10(c)
    # Energy in Msun*c**2
    ke  = np.float128(1.0) # Opacity in K10, ke = 0.2 cm**2/g -> k02 = 1, appropriate for e- scattering in fully ionized A/Z=2 elements
    k02 = np.float128(5.0) # Opacity in k/0.2cm**2/g
    #fp = 0.05 # Form factor 0.031 - 0.13 (RW11)
    v9 = np.float128(1.)   # velocity in 10**9 cm/s unit
    c = np.float(2.9979e8) # speed of light in m/s
    # wavelengh in angstrom, 1A = 10e-10 m
    sigb = np.float128(5.6704e-5) # Stephan-Boltzman constant in CGS [egs*cm**-2 * s**-1 * ]

    logLt = 43. + np.log10(2*r13) + 0.25*np.log10(Mc) + (7./4.)*np.log10(v9) + (-0.75)*np.log10(ke) + (-0.5)*np.log10(td)
    # total luminosity : early UV/optical luminosity emitted from ejecta diffusion including radioactive + cooling emission (This is observed data.)
    '''
    # (1) Observed total flux
    flux_obs = 2165.7
    # (2) Cooling emission
    flux_cooling = 839.4
    # (3) Ni (Simple power-law) flux
    flux_Ni = 1326.32

    ### Into AB mag to Flux
    #mbbflux = -2.5*np.log10(bbflux) + ff + abzero -2.5*fluxm
    mag_obs   = 18.913
    mag_Ni    = 19.193
    abzero       = np.float128(-48.600) # in erg*s**-1*cm**-2*Hz**-1
    flux_Ni_Jy   = 10.**(23-(mag_Ni+48.6)/2.5  ) # AB to Jansky
    flux_Ni_Hz   = flux_Ni_Jy*1.e-23 # Jansky to erg s-1 Hz-1 cm-2
    HzToAng      = (c/( (wav*1e-10)**2)) #
    flux_Ni_wav  = flux_Ni_Hz * HzToAng
    logflux_Ni   = np.log10(flux_Ni)
    ### Temp. of Observed flux (Early emission)
    logTeff  = (np.log10(flux_Ni_wav) - np.log10(sigb))/4.
    0.25*np.log10(2.*r13) = logTeff -np.log10(2.5) - 4. +(35./36.)*np.log10(ke) + (37./72.)*np.log10(td)
    ###
    log2r13 =  logLt - 43 -0.25*np.log10(Mc) -(7./4.)*np.log10(v9)+ 0.75*np.log10(ke) + (0.5)*np.log10(td)
    '''
    # ----
    logTeff = np.log10(2.5) + 4. + 0.25*np.log10(2.*r13) - (35./36.)*np.log10(ke) - (37./72.)*np.log10(td) # Effective temperature of the early emission

    logFbb  = np.log10(sigb) + 4.*logTeff # Blackbody flux of early emission (Ni decay part), bolometric
    d10pc   = np.float128(10. * 3.0856e18) # 10pc in cm to convert Luminosity to Flux
    fluxm   = (logLt - 2.*np.log10(d10pc) - np.log10(4.*np.pi)) - logFbb # Flux of shock heated emission = Total flux - blackbody component of early emission (Ni decay part), bolometric
    abzero  = np.float128(-48.600) # in erg*s**-1*cm**-2*Hz**-1

    teff = 10.**(logTeff)
    #xw   = 1000. + 40.*np.arange(200)
    #xw   = np.array(xw, dtype='float128')

    #from lgpy.sn2019ein_fitting import planck
    #bbflux  = planck(xw, teff)
    if  band == 'B' :
        xw = 4353.
        bbflux = planck(xw, teff)
    elif band == 'V' :
        xw = 5477.
        bbflux = planck(xw, teff)
    elif band == 'R' :
        xw = 6349.
        bbflux = planck(xw, teff)
    elif band == 'I' :
        xw = 8797.
        bbflux = planck(xw, teff)
    elif band == 'g' :
        xw = 4770.
        bbflux = planck(xw, teff)
    elif band == 'r' :
        xw = 6231.
        bbflux = planck(xw, teff)
    ff      = -2.5 * (2.*np.log10(xw) -10. -np.log10(c) ) # Angstrom to Hertz conversion factor
    mbbflux = -2.5*np.log10(bbflux) + ff + abzero -2.5*fluxm # convert flux to AB magnitude (Find ref in wiki), Bolometric.
    return mbbflux

'''
def MAG_Kasen(td, rstar, band):
    """
    fearly2_kasen
    #u.M_sun.cgs # Solar mass
    #c.c.cgs # Speed of light
    """
    import numpy as np
    import astropy.units as u
    from astropy import constants as c

    r13 = (rstar*u.R_sun/1.e+13).value

    eff = 1.0
    Mc  = 1.0/1.4
    dm  = np.log10(Msun) + 2. * np.log10(2.9979e10) - 51. # Difference btw E51 and Msun*c^2
'''
###################################################################################


def lcearly(rstar, band, fig=True):
    """
    Draw theoritical model of early light curve caused by Shock-heated emission.
    rstar : 0.1, 1
    band  : 'B', 'R'

    # Marion+16
    RG    = 2*1e13cm
    6M MS = 2*1e12cm
    2M MS = 5*1e11cm

    # Kasen 2010
    1-2M RG = ~1e13cm  (2 times)
    5-6M SG = 5*1e11cm (4 times)
    1-3M MS = 1-3e11cm (~2.5 times)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u
    #from lgpy.sn2019ein_fitting import fearly2_kasen

    #rstar = [0.1, 0.6, 1.0, round((5.*1.e+11*u.cm).to(u.R_sun).value,3),round((2.*1.e+12*u.cm).to(u.R_sun).value,3), round((2.*1.e+13*u.cm).to(u.R_sun).value,3)]
    #rstar = [0.1, 0.6, 1.0, round((5.*1.e+11*u.cm).to(u.R_sun).value,3),round((2.*1.e+12*u.cm).to(u.R_sun).value,3), round((2.*1.e+13*u.cm).to(u.R_sun).value,3)]
    #rstar = [0.1, 0.6, 1.0, round((5.*1.e+11*u.cm).to(u.R_sun).value,3),round((2.*1.e+12*u.cm).to(u.R_sun).value,3), 50.]
    td   = 0.01 + 0.01*np.arange(2000)
    y    = []

    for i in range(len(td)) :
        y_dum = fearly2_kasen(td[i], rstar, band)
        y.append(y_dum)

    #rmg = np.array(y) + 31.16
    #dl=[10.,20.,50.]
    #rmg1=np.array(y) + 5.*np.log10(dl[0])+25.
    #rmg2=np.array(y) + 5.*np.log10(dl[1])+25.
    #rmg3=np.array(y) + 5.*np.log10(dl[2])+25.

    fig, ax1 = plt.subplots(figsize=(6,5))
    ax1.set_ylim(np.max(y), np.min(y))
    ax1.plot(td, y, color='black', linewidth=2, linestyle='--')



#===============================================================================

#===============================================================================
# single filter, for SN 2018kp, used for 'R' band only

def EarlyFit2Pre(name='SN2017ein',incat_name='SN2017ein_B.txt',
                TBmax=58000, TBmaxErr=0.2, DM=30., D=0.5 FreeN=True) :
    """
    Perform preliminary Simple power-law fitting for the case 1 (2nd data < t < 8.3d) to decide the reference t0.
    True : power index n
    False : power index 2
    """
    from scipy.optimize import curve_fit
    case   = 2
    incat  = ascii.read(incat_name)
    f=extinction.Fitzpatrick99(1.7)
    Av=1.06
    A_lamda=f(np.array([4400,5400,6500,8100]),Av)
    AR_host = 0.67
    AR_mw = 0.062
    incat['MAG']=incat['MAG'] - AR_host - AR_mw
    incat['UL5']=incat['UL5'] - AR_host - AR_mw
    band='R'

    nodet=incat[incat['DET']=='N']
    incat0=incat[incat['DET']=='Y']
    # Polynomial or Spline fitting result
    # fitcat = ascii.read(fitcat_name)
    # TBmax & Distance
    #TBmax  = fitcat['Tmax'][0] ; TBmaxErr =  fitcat['TmaxErr'][0]
    # D_res  = GetDist(fitcat_name,Rv)
    #DM     = D_res['DM'] #distance measurement
    #D      = D_res['D']  # DM error
    # Early data fitting with simple power-law (Free power index)
    bandlist     = ['B', 'V', 'R']
    colorlist    = ['royalblue', 'seagreen', 'orange']#, 'orangered']
    if FreeN :
        output_name  = name+'-SPLFit2PreResN.dat'
        # alpha is not fixed to 2.
    else :
        output_name  = name+'-SPLFit2PreRes2.dat'
        # alpha is fixed to 2.


    f = open(output_name, 'w+')
    f.write('T0 T0Err alpha alphaErr mag0 mag0Err Trise TriseErr\n')
    #for cat, color in zip(catlist, colorlist) :
        # Select early data
    earlycat0    =  incat1[(incat1['MJD'] < incat1['MJD'][0] + 12) ]
    earlycat0    =  incat1[(incat1['MJD'] < incat1['MJD'][0] + 12) ]
    earlycat0    =  incat1[(incat1['MJD'] < 58153) ]
    earlycat     =  incat1[(incat1['MJD'] < 58154) ]
    #earlycat    =  earlycat0[:15]
    del earlycat[:2]
    # Fitting initial paramters
        #else : alpha= 2 set
        p0          = [57898., 19.]
        # Do Fitting
        popt2, pcov2  = curve_fit(SPL2, earlycat['MJD'], earlycat['MAG'], sigma=earlycat['MAGERR'],
                                p0=p0, absolute_sigma=True, maxfev=10000, check_finite=True)
        RedChiSq    = CalChiSq(SPL2(earlycat['MJD'], *popt2),  earlycat['MJD'], 
                                earlycat['MAG'], earlycat['MAGERR'], len(earlycat['MJD']), 1)
        # Fitting results
        t0       = popt2[0]
        t0Err    = np.sqrt(np.diag(pcov2))[0]
        mag0     = popt2[1]
        mag0Err  = np.sqrt(np.diag(pcov2))[1]
        dt0      = incat['MJD'][0] - t0
        Trise    = TBmax - t0
        TriseErr = np.sqrt(TBmaxErr**2 + t0Err**2)
        print('!====== Simple Power-law fitting results [{} band] ======!'.format(band))
        print('(1) Explosion time t0 = {}+/-{}'.format(round(t0,3), round(t0Err,3)))
        print('(1)-1. {} days before the 1st obs. :'.format(round(dt0, 3)))
        print('(2) Power index \u03B1 = 2 (fixed)')
        print('(3) mag0 = {}+/-{}'.format(round(mag0,3), round(mag0Err,3)))
        print('(4) Reduced \u03C7\u00b2 = {}'.format(round(RedChiSq, 3)))
        print('(5) Rise time trise = {}+/-{} days'.format(round(Trise,3), round(TriseErr, 3)))
        #f.write('{} {} {} {} {} {} {} {}\n'.format(t0, t0Err, 2.0, -99, mag0, mag0Err, Trise, TriseErr))
    #if FreeN :
        p0          = [57898., 2., 19.] # [T0, alpha, m0]
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
        print('!====== Simple Power-law fitting results [{} band] ======!'.format('B'))
        print('(1) Explosion time t0 = {}+/-{}'.format(round(t0,3), round(t0Err,3)))
        print('(1)-1. {} days before the 1st obs. :'.format(round(dt0, 3)))
        print('(2) Power index \u03B1 = {}+/-{}'.format(round(alpha,3), round(alphaErr,3)))
        print('(3) mag0 = {}+/-{}'.format(round(mag0,3), round(mag0Err,3)))
        print('(4) Reduced \u03C7\u00b2 = {}'.format(round(RedChiSq, 3)))
        print('(5) Rise time trise = {}+/-{} days'.format(round(Trise,3), round(TriseErr, 3)))
        #f.write('{} {} {} {} {} {} {} {}\n'.format(t0, t0Err, alpha, alphaErr, mag0, mag0Err, Trise, TriseErr))

    if FreeN :
        SPLXYcat_name = name+'-SPLXY_case2_{}N.dat'.format(band)
        resicat_name = name+'-SPLXY-Residual_case2_{}N.dat'.format(band)
    else :
        SPLXYcat_name = name+'-SPLXY_case2_{}2.dat'.format(band)
        resicat_name = name+'-SPLXY-Residual_case2_{}2.dat'.format(band)

    # Fitting lines & Residuals
    x = np.linspace(t0, np.max(earlycat['MJD']), num = len(earlycat['MJD']))
    y = SPL(x, *popt)
    xytbl  = Table({'x' : x, 'y' : y}, names=['x', 'y'])
    ascii.write(xytbl, output=SPLXYcat_name, overwrite=True)
    # Residuals earlycat[band]
    residual = earlycat['MAG'] - SPL(earlycat['MJD'], *popt)
    resitbl  = Table({'OBS' : earlycat['OBS'], 'xr' : earlycat['MJD'], 
                        'yr' : residual, 'yrerr' : earlycat['MAGERR']}, 
                        names=['OBS', 'xr', 'yr', 'yrerr'])
    ascii.write(resitbl, output=resicat_name, overwrite=True)
    f.close()


import matplotlib.pyplot as plt
plt.ion()
#nodet=ascii.read('non_det.txt')
plt.errorbar(earlycat0['MJD'],earlycat0['MAG'],yerr=earlycat0['MAGERR'],fmt='ro')
#SPL(earlycat['MJD'],*popt)
#plt.plot(earlycat['MJD'],SPL(earlycat['MJD'],*popt),'--')
x=np.linspace(58138,58154,200 )

plt.plot(x,SPL(x,*popt),linestyle='solid',label='alpha = '+str( round(popt[1],2)))
plt.plot(x,SPL2(x,*popt2),linestyle='solid',label='alpha = 2')

#plt.show()
plt.ylim(21,15)
from astropy.time import Time
#t = Time('2017-05-25T23:16:48', format='isot', scale='utc')
#discoverymjd=t.mjd
#discoverymjd
plt.vlines(discoverymjd,ymin=14,ymax=22, linestyle='dashdot',alpha=0.5)
plt.errorbar(nodet['MJD'],nodet['UL5'],yerr=0.2,lolims=True,fmt='_')
plt.xlim(58140,58154)
plt.ylim(21,15)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

plt.title('SN 2018kp EARLY LC FIT R-band',fontsize=16)


#from lgpy.sn2019ein_fitting import fearly2_kasen
#from lgpy.sn2019ein_fitting import fearly2_rw10
    # Shock heated cooling emission
td     = np.array(0.00 + 0.01*np.arange(401), dtype=np.float128)
#y_dict = {}
#print('R* = {} Rsol for {} band'.format(round(Rstar[r],3), band))
# Kasen

band='R'
Rstar=5
delay=0
y      = []
for i in range(len(td)):
    y_dum = fearly2_kasen(td[i], Rstar, band)
    y.append(y_dum)
#y_dict.update({Rstar : y})
plt.plot(td+t0+delay, np.array(y)+DM , linestyle='dashdot',label='RW '+str(Rstar)+' R ,'+str(delay)+' d')
# RW10
y_rw10 = []
for i in range(len(td)) :
    y_rw10_dum = fearly2_rw10(td[i], Rstar, band)
    y_rw10.append(y_rw10_dum)
plt.plot(td+t0+delay, np.array(y_rw10)+DM ,linestyle='dotted',label='Kasen '+str(Rstar)+' R ,'+str(delay)+' d')

plt.legend(loc='lower right')
