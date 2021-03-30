def GetDist(fitcat_name='SN2019ein-PolyFitRes.dat'):
    """
    DM = m - M
    DMErr = sqrt(mErr**2 + MErr**2)
    """
    fitcat = ascii.read(fitcat_name)    Bmax    = fitcat['MAGmax'][0]
    BmaxErr = fitcat['MmaxErr'][0]    AbsBmax    = fitcat['AbsMAGmax'][0]
    AbsBmaxErr = fitcat['AbsMAGmaxErr'][0]    Host    = GetEBVhost(fitcat_name)
    Ahost_B = Host['Ahost'][0]    DM    = Bmax - AbsBmax - Ahost_B
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
    n_free     : The number of parameters that we are fitting. If the model has no free parameter then dof = N. If we are fitting a model with n_free free parameters, we can force the model to exactly agree with n_free data points and dof = N - n_free.     '''
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

def GetEBVhost(fitcat_name = 'SN2019ein-PolyFitRes.dat') :
    """
    Bring the result of polynomial fitting.
    ==== Host galaxy extinction (Philips+99) ==='
    E(B-V)host = (Bmax - Vmax)c - (Bmax - Vmax)0
    (Bmax - Vmax)c = (Bmax -Vmax)obs - E(B-V)gal - K_BmaxVmax
    (Bmax - Vmax)0 = -0.07(+/-0.012) +0.114(+/-0.037)*(dB15-1.1), sigma   = 0.03
    BVRIJHK = 01234567
    """
    import extinction
    fitcat  = ascii.read('SN2019ein-PolyFitRes.dat')
    MAGmax  = fitcat['MAGmax']; MAGmaxErr = fitcat['MmaxErr']
    Bmax    = MAGmax[0]       ; Vmax      = MAGmax[1]
    BmaxErr = MAGmaxErr[0]    ; VmaxErr   = MAGmaxErr[1]
    dM15    = fitcat['dM15']  ; dM15Err   = fitcat['dM15Err']
    dB15    = dM15[0]         ; dB15Err  = dM15Err[0]    BVmaxObs      = Bmax - Vmax # Galactic extinction corrected
    # BVmaxObs here is BVmaxobs - (A_B-A_V) in the previous code.
    BVmaxObsErr   = np.sqrt(BmaxErr**2 + VmaxErr**2)
    K_BmaxVmax    = 0. # Assumption for a local object.
    BVmax0        = -0.07 + 0.114*(dB15 - 1.1)
    BVmax0_Sigma  = 0.03
    BVmax0Err     = np.sqrt(BVmax0_Sigma**2 + dB15Err**2)
    BVmaxc        = BVmaxObs - K_BmaxVmax
    EBVhost       = BVmaxc - BVmax0
    EBVhostErr    = np.sqrt(BVmaxObsErr**2 + BVmax0Err**2)    Rv_gal    = 3.1; Rv_HV = 1.55 # +/-0.06 (Wang+09)
    Ahost_V   = EBVhost * Rv_gal
    FilterSet = np.array([4353., 5477, 6349., 8797., 12483.0, 16313.0, 22010.0])
    Ahost     = extinction.fitzpatrick99(FilterSet, Ahost_V, Rv_gal)    print('*** Host Galaxy Extinction (Philips+99) Results ***')
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

def EarlyFit2Pre(incat_name='hgDcSN2019ein.LCfinal.txt', fitcat_name='SN2019ein-PolyFitRes.dat', FreeN=True) :
    """
    Perform preliminary Simple power-law fitting for the case 1 (2nd data < t < 8.3d) to decide the reference t0.
    True : power index n
    False : power index 2
    """
    case   = 2
    incat  = ascii.read(incat_name)
    # Polynomial or Spline fitting result
    fitcat = ascii.read(fitcat_name)
    # TBmax & Distance
    TBmax  = fitcat['Tmax'][0] ; TBmaxErr =  fitcat['TmaxErr'][0]
    D_res  = GetDist(fitcat_name=fitcat_name)
    DM     = D_res['DM']
    D      = D_res['D']
    # Early data fitting with simple power-law (Free power index)
    bandlist     = ['B', 'V', 'R', 'I']
    colorlist    = ['royalblue', 'seagreen', 'orange', 'orangered']
    if FreeN :
        output_name  = 'SN2019ein-SPLFit2PreResN.dat'
        # alpha is not fixed to 2.
    else :
        output_name  = 'SN2019ein-SPLFit2PreRes2.dat'
        # alpha is fixed to 2.
    f = open(output_name, 'w+')
    f.write('T0 T0Err alpha alphaErr mag0 mag0Err Trise TriseErr\n')
    for band, color in zip(bandlist, colorlist) :
        # Select early data
        earlycat    =  incat[(incat['MJD'] >= 58603) & (incat['MJD'] < incat['MJD'][0] + 8.3) & (incat[band] != -99)]
        # Fitting initial paramters
        if FreeN :
            p0          = [58603., 2., 19.]
            # Do Fitting
            popt, pcov  = curve_fit(SPL, earlycat['MJD'], earlycat[band], sigma=earlycat[band+'err'], p0=p0, absolute_sigma=True, maxfev=10000, check_finite=True)
            RedChiSq    = CalChiSq(SPL(earlycat['MJD'], *popt),  earlycat['MJD'], earlycat[band], earlycat[band+'err'], len(earlycat['MJD']), 2)
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
            print('!====== Simple Power-law fitting results [{} band] ======!'.format(band))
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
