def DoSpline(x_data, y_data, yerr_data, s=None, k=5):
    from scipy.interpolate import UnivariateSpline
    Spline = UnivariateSpline(x_data, y_data, s=s, k=k)
    ys     = Spline(x_data)
    #ReChisq  = CalChiSq(Spline(x_data), x_data, y_data, yerr_data, len(x_data),  len(x_data)-1)
   # print('Reduced Chi Sq. = {}'.format(ReChisq))
    return [Spline]
def SplineBTS(incat_name='gDcSN2019ein.LCfinal.txt',
            nircat_name= 'gSN2019ein.NIR.LCfinal.txt', repeats=1000, verbose=False):
    from tqdm import tqdm
    import scipy.stats as stats
    from astropy.table import Table
    from astropy.stats import sigma_clip
    from astropy.modeling import models, fitting    incat  = ascii.read(incat_name)  # Optical bands
    nircat = ascii.read(nircat_name) # NIR bands
    f = open('SN2019ein-SplineFitRes.dat', 'w+')
    #f.write('Tmax TmaxErr MAGmax MmaxErr dM15 dM15Err AbsMAGmax AbsMAGmaxErr DM DMErr\n')
    f.write('Tmax TmaxErr MAGmax MmaxErr dM15 dM15Err AbsMAGmax AbsMAGmaxErr\n')    # Spline interpolation of Phillips Curve
    PhCurve =  ascii.read('/data1/SN2019ein/work/philips/philips_curve.txt', guess=True)
    PhX     = PhCurve['X'] ; PhY = PhCurve['Y']
    from scipy.interpolate import UnivariateSpline
    f_Ph    = UnivariateSpline(PhCurve['X'], PhCurve['Y'], s=None, k=5)
    xnew = np.linspace(np.min(PhX), np.max(PhX), num=51, endpoint=True)
    for band in tqdm(['B', 'V', 'R', 'I', 'J', 'H', 'K']) :
        Tmax_Kawabata = 58618.24
        if band in ['B', 'V'] :
            k=5
            idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 5) & (incat['MJD'] < Tmax_Kawabata + 18))[0]
        elif band in ['R'] :
            k=5
            idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 5) & (incat['MJD'] < Tmax_Kawabata + 18))[0]
        elif band in ['I'] :
            k=5
            idx   = np.where( (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 7) & (incat['MJD'] < Tmax_Kawabata + 18))[0]
        elif band in ['J', 'H', 'K'] :
            k=5
            incat = nircat
            idx   = np.where(  (incat[band] != -99) & (incat['MJD'] > Tmax_Kawabata - 13) & (incat['MJD'] < Tmax_Kawabata + 8))[0]
        x_data    = incat['MJD'][idx]
        y_data    = incat[band][idx]
        yerr_data = incat[band+'err'][idx]
        tbl          = Table({'x' : x_data, 'y' : y_data, 'yerr': yerr_data}, names = ['x', 'y', 'yerr'])
        N_data       = len(x_data)        Sampled_data = np.random.choice(tbl, size=(repeats, N_data), replace=True)
                        # 중복 허용하여 repeats번 만큼 임의로 뽑는다.        print('!====== Bootstap resampling ======!')
        print('Band for bootstrap : {}'.format(band))
        print('Number of samples  : {}'.format(N_data))
        print('Repeats            : {}'.format(repeats))
        print('!=================================!')        Tmax_array, MAGmax_array, dM15_array = [], [], []
        AbsMAGmax_array, DM_array = [], []
        for i in range(repeats) :
            sample = Sampled_data[i] # i 번째 반복에서 피팅할 데이터
            sample = np.unique(sample[np.argsort(sample['x'])]) # 시간으로 정렬하고 중복을 없앤 i번째 데이터
            if verbose :
                print('N of data without overlapping : {}'.format(len(sample)))
            Spline = UnivariateSpline(sample['x'], sample['y'], s=None, k=5)
            xs     =  np.linspace(np.min(sample['x']), np.max(sample['x']),
                                    num= int((np.max(sample['x'])-np.min(sample['x']))/(5.*(1./24)*(1./60))))
                                    # 5 min interval
            ys     =  Spline(xs)
            Tmax   = xs[ys == np.min(ys)]
            MAGmax = ys[ys == np.min(ys)]
            if Tmax.size > 1 :
                Tmax   = np.mean(Tmax)
                MAGmax = np.mean(MAGmax)
            else:
                Tmax   = Tmax[0]
                MAGmax = MAGmax[0]
                pass
            Tmax_array.append(Tmax)
            MAGmax_array.append(MAGmax)
            MAG15  = Spline(Tmax + 15)
            dM15   = MAG15.base[0] - MAGmax      ; dM15_array.append(dM15)
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
    f.close()
    print('SN2019ein-SplineFitRes.dat is created.')
    return
