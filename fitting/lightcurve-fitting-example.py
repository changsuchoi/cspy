
from lightcurve_fitting.lightcurve import LC
from pkg_resources import resource_filename
from lightcurve_fitting import lightcurve, models, fitting, bolometric
from lightcurve_fitting.models import ShockCooling2
from lightcurve_fitting.models import ShockCooling
from lightcurve_fitting.models import CompanionShocking
from lightcurve_fitting.fitting import lightcurve_mcmc, lightcurve_corner
from lightcurve_fitting import fitting
from lightcurve_fitting.bolometric import calculate_bolometric, plot_bolometric_results, plot_color_curves
import matplotlib.pyplot as plt


redshift = 0.002
redshift = 0.01042
filename='SN2018kp_lightcurve_fitting_input.dat'
lc = lightcurve.LC.read(filename,format='ascii')

lc.meta['dm']=33.44
lc.meta['extinction']={
'B':0.104,
'V':0.079,
'R':0.062,
'I':0.043,
'g':0.095,
'r':0.066,
'i':0.049}
# RV=1.766, E(B-V)=0.68, Fitzpatrick99
lc.meta['hostext']={
'B':1.631,
'V':1.066,
'R':0.670,
'I':0.457,
'g':0.939,
'r':0.479,
'i':0.328}

lc.calcAbsMag()
lc.calcLum()
lc.plot(xcol='MJD')

# BOLOMETRIC
#===========================================================================
outpath = '/data7/cschoi/sngal/NGC3367/phot/cut/fitting/lcfitting/bolometric'
t = calculate_bolometric(lc, redshift, outpath, colors=['B-V', 'V-R', 'B-R','R-I', 'g-r', 'r-i'])
print(t)
plot_bolometric_results(t)
plot_color_curves(t)


# Fit only the early light curve
lc_early = lc.where(MJD_min=58142., MJD_max=58154.)

# Define the priors and initial guesses
# ShockCooling2
#p_min = [0., 0., 0., 57468.]
#p_max = [100., 100., 100., 57468.7]
#p_lo = [20., 2., 20., 57468.5]
#p_up = [50., 5., 50., 57468.7]

p_min = [0., 0., 0., 58142.]
p_max = [100., 100., 100., 58142.2]
p_lo = [20., 2., 20., 58142.0]
p_up = [50., 5., 50., 58142.2]

models=['ShockCooling','ShockCooling2','CompanionShocking',]
sampler = fitting.lightcurve_mcmc(lc_early, ShockCooling2, model_kwargs={'z': redshift,'RW':True},
                              p_min=p_min, p_max=p_max, p_lo=p_lo, p_up=p_up,
                              nwalkers=100, nsteps=1000, nsteps_burnin=1000, show=True)
lightcurve_corner(lc_early, ShockCooling2, sampler.flatchain, model_kwargs={'z': redshift,'RW':True})

p_min = [0., 0., 0.,0, 58142.]
p_max = [100., 100., 100.,100, 58142.2]
p_lo = [20., 2.,2, 20., 58142.0]
p_up = [50., 5.,5, 50., 58142.2]
models=['ShockCooling','ShockCooling2','CompanionShocking',]
sampler = fitting.lightcurve_mcmc(lc_early, ShockCooling, model_kwargs=
                {'z': redshift,'RW':True},
                p_min=p_min, p_max=p_max, p_lo=p_lo, p_up=p_up,
                nwalkers=100, nsteps=1000, nsteps_burnin=1000, show=True)
lightcurve_corner(lc_early, ShockCooling, sampler.flatchain, model_kwargs={'z': redshift,'RW':True})


#SiFTo UBVgri only
lc0=lc[lc['filt']!='R']
lc0=lc0[lc0['filt']!='I']

lc_early = lc0.where(MJD_min=58142., MJD_max=58150.)
p_min = [0.,   0.,   0.,   0,  0,    0,    0, 58142.]
p_max = [100., 100., 100.,100, 100., 100.,100, 58142.2]
p_lo =  [20.,  2.,    2,   20., 1,    1,   1, 58142.0]
p_up =  [50.,  5.,    5,   50., 100,  100, 100, 58142.2]
models=['ShockCooling','ShockCooling2','CompanionShocking',]
sampler = fitting.lightcurve_mcmc(lc_early, CompanionShocking, model_kwargs=
                {'z': redshift},
                p_min=p_min, p_max=p_max, p_lo=p_lo, p_up=p_up,
                nwalkers=100, nsteps=1000, nsteps_burnin=1000, show=True)
lightcurve_corner(lc_early, CompanionShocking, sampler.flatchain, model_kwargs={'z': redshift})
