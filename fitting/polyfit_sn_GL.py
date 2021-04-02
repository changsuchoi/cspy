def GetDist(fitcat_name='SN2019ein-PolyFitRes.dat'):
    """
    DM = m - M
    DMErr = sqrt(mErr**2 + MErr**2)
    """
    fitcat = ascii.read(fitcat_name)
    Bmax    = fitcat['MAGmax'][0]
    BmaxErr = fitcat['MmaxErr'][0]
    AbsBmax    = fitcat['AbsMAGmax'][0]
    AbsBmaxErr = fitcat['AbsMAGmaxErr'][0]
    Host    = GetEBVhost(fitcat_name)
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


def DoPoly(x_data, y_data, yerr_data, degree=7):
    import scipy.stats as stats
    from astropy.stats import sigma_clip
    from astropy.modeling import models, fitting
    model    = models.Polynomial1D(degree=degree)
    fitter   = fitting.LevMarLSQFitter()
    best_fit = fitter(model, x_data, y_data)
    ReChisq  = CalChiSq(best_fit(x_data), x_data, y_data, yerr_data, len(x_data), degree+1)
    #print(best_fit)
    print('Reduced Chi Sq. = {}'.format(ReChisq))
    return [best_fit, ReChisq]

def PolyBTS(incat_name='gDcSN2019ein.LCfinal.txt',
            nircat_name= 'gSN2019ein.NIR.LCfinal.txt',
            repeats=1000, verbose=False) :
    """
    Perform 1d polynomial fitting for input catalog.
    idx   = np.where(  (nircat[band] != -99) & (nircat['MJD'] > 58619.35 - 10) & (nircat['MJD'] < 58619.35 + 16))[0]
    """
    from tqdm import tqdm
    import scipy.stats as stats
    from astropy.table import Table
    from astropy.stats import sigma_clip
    from astropy.modeling import models, fitting
    import astropy.io.ascii as ascii
    import numpy as np


    incat  = ascii.read(incat_name)  # Optical bands
    # nircat = ascii.read(nircat_name) # NIR bands

    # Spline interpolation of Phillips Curve
    PhCurve =  ascii.read('philips_curve.txt', guess=True)
    PhX     = PhCurve['X'] ; PhY = PhCurve['Y']
    from scipy.interpolate import UnivariateSpline
    f_Ph    = UnivariateSpline(PhCurve['X'], PhCurve['Y'], s=None, k=5)
    xnew = np.linspace(np.min(PhX), np.max(PhX), num=51, endpoint=True)
    xp_array, yp_array = [], []
    bandlist  = ['B', 'V', 'R', 'I', 'J', 'H', 'K']
    bandlist  = ['B', 'V', 'R', 'I','g','r','i']#, 'J', 'H', 'K']
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
    V1cat=Vcat[(Vcat['MJD']>Tmax-10) & (Vcat['MJD']<Tmax+20)]
    x_data    = V1cat['MJD']
    y_data    = V1cat['MAG']
    yerr_data = V1cat['MAGERR']
    R1cat=Rcat[(Rcat['MJD']>Tmax-10) & (Rcat['MJD']<Tmax+20)]
    x_data    = R1cat['MJD']
    y_data    = R1cat['MAG']
    yerr_data = R1cat['MAGERR']
    I1cat=Icat[(Icat['MJD']>Tmax-10) & (Icat['MJD']<Tmax+20)]
    x_data    = I1cat['MJD']
    y_data    = I1cat['MAG']
    yerr_data = I1cat['MAGERR']
    f = open('SN2018kp-PolyFitRes.dat', 'w+')
    f.write('Tmax TmaxErr MAGmax MmaxErr dM15 dM15Err AbsMAGmax AbsMAGmaxErr\n')
    for band in tqdm(bandlist) :
        Tmax = 58159.5 # 58160 (sn 2018kp, around peak)
        if band == 'B' :
            degree = 7
            x_data    = B1cat['MJD']
            y_data    = B1cat['MAG']
            yerr_data = B1cat['MAGERR']
            #idx   = np.where( (incat['MJD'] > Tmax - 10) & (incat['MJD'] < Tmax + 20))[0]
        elif band == 'V' :
            degree = 7
            x_data    = V1cat['MJD']
            y_data    = V1cat['MAG']
            yerr_data = V1cat['MAGERR']
            #idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax - 10) & (incat['MJD'] < Tmax + 20))[0]
        elif band == 'R' :
            degree = 9
            x_data    = R1cat['MJD']
            y_data    = R1cat['MAG']
            yerr_data = R1cat['MAGERR']
            #idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax - 10) & (incat['MJD'] < Tmax + 20))[0]
        elif band == 'I' :
            degree = 9
            x_data    = I1cat['MJD']
            y_data    = I1cat['MAG']
            yerr_data = I1cat['MAGERR']

            #idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax - 10) & (incat['MJD'] < Tmax + 20))[0]
        #elif band in ['J', 'H', 'K'] :
        #    degree = 3
        #    incat = nircat
        #    idx   = np.where(  (incat[band] != -99) & (incat['MJD'] > Tmax - 13) & (incat['MJD'] < Tmax + 7))[0]
        #x_data    = incat['MJD'][idx]
        #y_data    = incat[band][idx]
        #yerr_data = incat[band+'err'][idx]

        tbl          = Table({'x' : x_data, 'y' : y_data, 'yerr': yerr_data}, names = ['x', 'y', 'yerr'])
        N_data       = len(x_data)
        Sampled_data = np.random.choice(tbl, size=(repeats, N_data), replace=True) # ?? ???? repeats? ?? ??? ???.
        print('!====== Bootstap resampling ======!')
        print('Band for bootstrap : {}'.format(band))
        print('Number of samples  : {}'.format(N_data))
        print('Repeats            : {}'.format(repeats))
        print('!=================================!')
        Tmax_array, MAGmax_array, dM15_array = [], [], []
        AbsMAGmax_array, DM_array = [], []
        for i in range(repeats) :
            sample = Sampled_data[i] # i ?? ???? ??? ???
            sample = np.unique(sample[np.argsort(sample['x'])]) # ???? ???? ??? ?? i?? ???
            if verbose :
                print('N of data without overlapping : {}'.format(len(sample)))
            PolyRes  = DoPoly(sample['x'], sample['y'], sample['yerr'], degree=degree)
            best_fit = PolyRes[0]
            ReChisq  = PolyRes[1]
            if ReChisq == -99 :
                pass
            else :
                xp     =  np.linspace(np.min(sample['x']), np.max(sample['x']),
                                    num= int((np.max(sample['x'])-np.min(sample['x']))/(5.*(1./24)*(1./60)))) # 1min interval
                yp     =  best_fit(xp)
                if i == 0 :
                    xp_array = list(xp)
                    yp_array = list(yp)
                Tmax   = xp[yp == np.min(yp)]
                MAGmax = yp[yp == np.min(yp)]
                if Tmax.size > 1 :
                    Tmax   = np.mean(Tmax)
                    MAGmax = np.mean(MAGmax)
                else:
                    Tmax   = Tmax[0]
                    MAGmax = MAGmax[0]
                    pass
                Tmax_array.append(Tmax)
                MAGmax_array.append(MAGmax)
                MAG15  = best_fit(Tmax + 15)
                dM15   = MAG15 - MAGmax      ; dM15_array.append(dM15)
                if band == 'B' :
                    AbsMAGmax = f_Ph(dM15); AbsMAGmax_array.append(AbsMAGmax.base[0])
                    DM = MAGmax - AbsMAGmax.base[0]; DM_array.append(DM)
                else :
                    pass
                Tmax_array.sort(); MAGmax_array.sort(); dM15_array.sort()
        MAGmax = np.percentile(MAGmax_array, [50.0])[0]
        Tmax   = np.percentile(Tmax_array, [50.0])[0]
        dM15   = np.percentile(dM15_array, [50.0])[0]
        TmaxErr   = np.std(Tmax_array)
        MAGmaxErr = np.std(MAGmax_array)
        dM15Err   = np.std(dM15_array)
        Polytbl = Table({'xp' : xp_array, 'yp' : yp_array}, names=['xp', 'yp'])
        ascii.write(Polytbl, output='SN2018kp-PolyXY_{}.dat'.format(band), overwrite=True)
        if band == 'B' :
            AbsMAGmax    = np.percentile(AbsMAGmax_array, [50.0])[0]
            AbsMAGmaxErr = np.std(AbsMAGmax_array)
            DM    = np.percentile(DM_array, [50.0])[0]
            DMErr = np.std(DM_array)
        else :
            AbsMAGmax    = -99
            AbsMAGmaxErr = -99
            DM    = -99
            DMErr = -99
        f.write('{} {} {} {} {} {} {} {}\n'.format(round(Tmax, 3), round(TmaxErr, 3), round(MAGmax, 3), round(MAGmaxErr, 3), round(dM15, 3), round(dM15Err, 3), round(AbsMAGmax,3), round(AbsMAGmaxErr,3)))
        #f.write('{} {} {} {} {} {} {} {} {} {}\n'.format(round(Tmax, 3), round(TmaxErr, 3), round(MAGmax, 3), round(MAGmaxErr, 3), round(dM15, 3), round(dM15Err, 3), round(AbsMAGmax,3), round(AbsMAGmaxErr,3), round(DM, 3), round(DMErr, 3) ))
    f.close()
    print('SN2018kp-PolyFitRes.dat is created.')
