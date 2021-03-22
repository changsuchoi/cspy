def PolyBTS(incat_name='gDcSN2019ein.LCfinal.txt', nircat_name= 'gSN2019ein.NIR.LCfinal.txt', repeats=1000, verbose=False) :
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
    f = open('SN2018kp-PolyFitRes.dat', 'w+')
    f.write('Tmax TmaxErr MAGmax MmaxErr dM15 dM15Err AbsMAGmax AbsMAGmaxErr\n')
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

    for band in tqdm(bandlist) :
        Tmax_Kawabata = 58618.24 # 58160 (sn 2018kp, around peak)
        if band in ['B', 'V'] :
            degree = 7
            idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 5) & (incat['MJD'] < Tmax_Kawabata + 18))[0]
        elif band in ['R'] :
            degree = 9
            idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 5) & (incat['MJD'] < Tmax_Kawabata + 18))[0]
        elif band in ['I'] :
            degree = 9
            idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 7) & (incat['MJD'] < Tmax_Kawabata + 18))[0]
        #elif band in ['J', 'H', 'K'] :
        #    degree = 3
        #    incat = nircat
        #    idx   = np.where(  (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 13) & (incat['MJD'] < Tmax_Kawabata + 7))[0]
        x_data    = incat['MJD'][idx]
        y_data    = incat[band][idx]
        yerr_data = incat[band+'err'][idx]
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
        ascii.write(Polytbl, output='SN2019ein-PolyXY_{}.dat'.format(band), overwrite=True)
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
    print('SN2019ein-PolyFitRes.dat is created.')
