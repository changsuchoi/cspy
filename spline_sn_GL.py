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
