
from lightcurve_fitting.lightcurve import LC
from pkg_resources import resource_filename


from lightcurve_fitting.bolometric import calculate_bolometric, plot_bolometric_results, plot_color_curves

redshift = 0.002
outpath = '/Users/griffin/Desktop/SN2016bkv_bolometric'
t = calculate_bolometric(lc, redshift, outpath, colors=['B-V', 'g-r', 'r-i'])
print(t)
plot_bolometric_results(t)
plot_color_curves(t)


from lightcurve_fitting.models import ShockCooling2
from lightcurve_fitting.models import ShockCooling
from lightcurve_fitting.fitting import lightcurve_mcmc, lightcurve_corner
from lightcurve_fitting import fitting

# Fit only the early light curve
lc_early = lc.where(MJD_min=57468., MJD_max=57485.)

# Define the priors and initial guesses
p_min = [0., 0., 0., 57468.]
p_max = [100., 100., 100., 57468.7]
p_lo = [20., 2., 20., 57468.5]
p_up = [50., 5., 50., 57468.7]

redshift = 0.002

sampler = fitting.lightcurve_mcmc(lc_early, ShockCooling, model_kwargs={'z': redshift},
                              p_min=p_min, p_max=p_max, p_lo=p_lo, p_up=p_up,
                              nwalkers=10, nsteps=100, nsteps_burnin=100, show=True)
lightcurve_corner(lc_early, ShockCooling, sampler.flatchain, model_kwargs={'z': redshift})
