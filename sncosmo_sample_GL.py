# SNCOSMO LC fitting

def sncosmoTable(incat_name='gDcSN2019ein.LCfinal.txt'):
    incat = ascii.read(incat_name)
    FLT   = ['B','V','R','I']
    zp0 = 25.
    incat['Scale_Factor_B'] = 10.**(0.4*(zp0 - incat['ZP_B']))
    incat['Scale_Factor_V'] = 10.**(0.4*(zp0 - incat['ZP_V']))
    incat['Scale_Factor_R'] = 10.**(0.4*(zp0 - incat['ZP_R']))
    incat['Scale_Factor_I'] = 10.**(0.4*(zp0 - incat['ZP_I']))
    incat['FLUXCAL_B'] = 10.**(-0.4*(incat['B']-incat['ZP_B'])) * incat['Scale_Factor_B']
    incat['FLUXCAL_V'] = 10.**(-0.4*(incat['V']-incat['ZP_V'])) * incat['Scale_Factor_V']
    incat['FLUXCAL_R'] = 10.**(-0.4*(incat['R']-incat['ZP_R'])) * incat['Scale_Factor_R']
    incat['FLUXCAL_I'] = 10.**(-0.4*(incat['I']-incat['ZP_I'])) * incat['Scale_Factor_I']
    incat['FLUXCALERR_B'] = np.abs(incat['FLUXCAL_B'])*np.abs(-0.4*np.log(10)*np.sqrt(incat['Berr']**2+incat['ZPERR_B']**2))
    incat['FLUXCALERR_V'] = np.abs(incat['FLUXCAL_V'])*np.abs(-0.4*np.log(10)*np.sqrt(incat['Verr']**2+incat['ZPERR_V']**2))
    incat['FLUXCALERR_R'] = np.abs(incat['FLUXCAL_R'])*np.abs(-0.4*np.log(10)*np.sqrt(incat['Rerr']**2+incat['ZPERR_R']**2))
    incat['FLUXCALERR_I'] = np.abs(incat['FLUXCAL_I'])*np.abs(-0.4*np.log(10)*np.sqrt(incat['Ierr']**2+incat['ZPERR_I']**2))
    f = open('sncosmo_SN2019ein.txt', 'w+')
    f.write('mjd band flux fluxerr zp zpsys\n')
    for i in range(len(incat)):
        for band in FLT :
            #idx  = np.where(incat['ZP_'+band] != -99)[0]
            if band == 'B' :
                #bandname = 'bessellb'
                bandname = 'standard::b'
            elif band == 'V':
                #bandname = 'bessellv'
                bandname = 'standard::v'
            elif band == 'R' :
                #bandname = 'bessellr'
                bandname = 'standard::r'
            elif band == 'I' :
                #bandname = 'besselli'
                bandname = 'standard::i'
            if incat['ZP_'+band][i] != -99 :
                #print('B : {} x {} = {}'.format(round(incat['FLUXCAL_'+band][i],3), round(incat['Scale_Factor_'+band][i],3), round(incat['FLUXCAL_'+band][i]*incat['Scale_Factor_'+band][i],3)))
                #incat['FLUXCAL_'+band][i] = incat['FLUXCAL_'+band][i] * incat['Scale_Factor_'+band][i]
                f.write('{} {} {} {} {} ab\n'.format(incat['MJD'][i], bandname, incat['FLUXCAL_'+band][i], incat['FLUXCALERR_'+band][i], zp0))
                #f.write('{} {} {} {} {} ab\n'.format(incat['MJD'][idx][i], band, incat['FLUXCAL_'+band][idx][i], incat['FLUXCALERR_'+band][idx][i], incat['ZP_'+band][idx][i]))
            else :
                pass
    f.close()


def sncosmoFit(datacat_name='sncosmo_SN2019ein.exact.txt'):
    """
    z  =  0.007166 ; NED
    """
    import sncosmo
    data   = ascii.read(datacat_name)
    magsys = sncosmo.CompositeMagSystem( bands={'standard::b' : ('ab', 9.851778333549941), 'standard::v' : ('ab', 10.165691850734973), 'standard::r' : ('ab', 9.780891653643735), 'standard::i' : ('ab', 10.00617773098994)    })
    dust   = sncosmo.F99Dust(r_v=3.1)
    model  = sncosmo.Model(source='salt2', effects=[dust], effect_names=['host'], effect_frames=['rest'])
    model.set(z = 0.007166)
    #model.set_source_peakabsmag(-19,'standard::b', 'ab')
    model.param_names
    model.parameters
    print('Set the model like this;', model)
    #source = sncosmo.SALT2Source(modeldir='/home/lim9/.astropy/cache/sncosmo/models/salt2/salt2-4')
    ### Applying cuts
    minsnr = 3
    tl     = data['mjd'][0]
    tu     = data['mjd'][0] +50.
    #mask = (data['mjd'] >= tl) & (data['mjd'] < tu)
    #mask = (data['mjd'] >= tl) & (data['mjd'] < tu) & (data['flux'] / data['fluxerr'] > minsnr)
    #data   = data[mask]
    params_to_fit = ['t0', 'x0', 'x1', 'c', 'hostebv']
    bounds = {'hostebv':(0, 0.2)}
    #bounds = { 't0':(58618, 58620), 'x0':(0.01,0.05), 'x1':(-2,-0.5), 'c':(0.0, 0.15), 'hostebv':(0, 0.1)}
    result, fitted_model = sncosmo.fit_lc(data, model, params_to_fit,  guess_amplitude=False, guess_t0=True, guess_z=False, bounds=bounds, method='minuit', modelcov=True, phase_range=(-15, 60), verbose=True)
    print("Number of chi^2 function calls :", result.ncall)
    print("Number of degrees of freedom in fit :", result.ndof)
    print("Chi^2 value at minimum :", result.chisq)
    print("Reduced Chi^2          :", result.chisq/result.ndof)
    print("Model parameters :", result.param_names)
    print("Best-fit values  :", result.parameters)
    print("The result contains the following attributes:\n", result.keys())
    sncosmo.plot_lc(data, model=fitted_model, errors=result.errors, pulls=True, figtext='SN 2019ein; SALT2', show_model_params=True)
