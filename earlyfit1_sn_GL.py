def EarlyFit1(incat_name='hgDcSN2019ein.LCfinal.txt', fitcat_name='SN2019ein-PolyFitRes.dat'):
    """
    Minimize the sum of the square of the residual.
    """
    case = 1
    # Polynomial or Spline fitting result
    fitcat = ascii.read(fitcat_name)
    # TBmax & Distance
    TBmax  = fitcat['Tmax'][0] ; TBmaxErr =  fitcat['TmaxErr'][0]
    D_res  = GetDist(fitcat_name=fitcat_name)
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
