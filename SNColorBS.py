def MakeColorBS(incat_name='hgSN2021hpr.LCfinal.txt', band1='B', band2='V', sameobs=True):
    """
    Make color catalog using bisect.
    Put host, mw reddening-corrected catalog

    Sameobs == True : 같은 천문대끼리 묶어줌.
    """
    import bisect
    incat         = ascii.read(incat_name)
    if sameobs :
        obslist = ['SAO_KL4040', 'DOAO', 'LOAO', 'SOAO', 'CBNUO', 'CMO']
        cat     = Table(names=('mjd', 'obs', band1+'-'+band2, 'sigma'+band1+band2), dtype=['f', 'S', 'f', 'f'])
        for obs in obslist :
            print('For '+obs)
            idx1 = np.where( (incat['filter'] == band1) & (incat['mag'] != -99) & (incat['det'] == 1) & (incat['obs'] == obs))[0]
            idx2 = np.where( (incat['filter'] == band2) & (incat['mag'] != -99) & (incat['det'] == 1) & (incat['obs'] == obs))[0]
            for i in range(len(incat[idx2])):
                b = bisect.bisect_right(incat[idx1]['mjd'], incat[idx2]['mjd'][i])
                if b == 0 :
                    pass
                else :
                    cat.add_row((np.mean([incat[idx1]['mjd'][b-1], incat[idx2]['mjd'][i]]), incat[idx2]['obs'][i], incat[idx1]['mag'][b-1] - incat[idx2]['mag'][i], np.sqrt(incat[idx1]['magerr'][b-1]**2 + incat[idx2]['magerr'][i]**2)))
    elif sameobs == False :
        cat     = Table(names=('mjd', band1+'-'+band2, 'sigma'+band1+band2), dtype=['f', 'f', 'f'])
        idx1 = np.where( (incat['filter'] == band1) & (incat['mag'] != -99) & (incat['det'] == 1) )[0]
        idx2 = np.where( (incat['filter'] == band2) & (incat['mag'] != -99) & (incat['det'] == 1) )[0]
        for i in range(len(incat[idx2])):
            b = bisect.bisect_right(incat[idx1]['mjd'], incat[idx2]['mjd'][i])
            if b == 0 :
                pass
            else :
                cat.add_row((np.mean([incat[idx1]['mjd'][b-1], incat[idx2]['mjd'][i]]), incat[idx1]['mag'][b-1] - incat[idx2]['mag'][i], np.sqrt(incat[idx1]['magerr'][b-1]**2 + incat[idx2]['magerr'][i]**2)))        
    print('Merging cat and sorting...')
    cat.sort('mjd')  
    if sameobs : 
        ascii.write(cat, 'SN2021hpr-'+band1+band2+'.cat', names=['mjd', 'obs',band1+'-'+band2, 'sigma'+band1+band2], overwrite=True)
    elif sameobs == False :
        ascii.write(cat, 'SN2021hpr-'+band1+band2+'.cat', names=['mjd', band1+'-'+band2, 'sigma'+band1+band2], overwrite=True)
    print( 'SN2021hpr-'+band1+band2+'.cat is created.')
