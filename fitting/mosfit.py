mosfit -m ia -e ../SN2018kp_lightcurve_fitting_input.dat -s ia -c


%matplotlib inline
%config InlineBackend.figure_format='retina'

import corner
import json
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tqdm import tqdm_notebook
from collections import OrderedDict
from mosfit.plotting import bandcolorf

sns.reset_orig()

plt.rcParams["font.family"] = "serif"
plt.rcParams.update({'font.size': 14})

with open('../products/walkers.json', 'r', encoding = 'utf-8') as f:
    data = json.loads(f.read())
    if 'name' not in data:
        data = data[list(data.keys())[0]]

photo = data['photometry']
model = data['models'][0]

real_data = len([x for x in photo if 'band' in x and 'magnitude' in x and (
    'realization' not in x or 'simulated' in x)]) > 0

band_attr = ['band', 'instrument', 'telescope', 'system', 'bandset']
band_list = list(set([tuple(x.get(y, '')
                            for y in band_attr) for x in photo
                            if 'band' in x and 'magnitude' in x]))
real_band_list = list(set([tuple(x.get(y, '')
                                 for y in band_attr) for x in photo
                                 if 'band' in x and 'magnitude' in x and (
                                     'realization' not in x or 'simulated' in x)]))
xray_instrument_attr = ['instrument', 'telescope']
xray_instrument_list = list(set([tuple(x.get(y, '')
                            for y in xray_instrument_attr) for x in photo
                            if 'instrument' in x and 'countrate' in x]))
real_band_list = list(set([tuple(x.get(y, '')
                                 for y in band_attr) for x in photo
                                 if 'band' in x and 'magnitude' in x and (
                                     'realization' not in x or 'simulated' in x)]))
real_xray_instrument_list = list(set([tuple(x.get(y, '')
                                 for y in xray_instrument_attr) for x in photo
                                 if 'instrument' in x and 'countrate' in x and (
                                     'realization' not in x or 'simulated' in x)]))





#==============================================================================
# Uncomment line below to only plot from the specified instruments.
# inst_exclusive_list = ['UVOT']

fig = plt.figure(figsize=(12,8))
plt.gca().invert_yaxis()
plt.gca().set_xlim(55700,55820)
plt.gca().set_ylim(bottom=25, top=19)
plt.gca().set_xlabel('MJD')
plt.gca().set_ylabel('Apparent Magnitude')
used_bands = []
for full_band in tqdm_notebook(band_list, desc='Photo', leave=False):
    (band, inst, tele, syst, bset) = full_band
    try:
        inst_exclusive_list
    except:
        pass
    else:
        if inst not in inst_exclusive_list:
            continue
    extra_nice = ', '.join(list(filter(None, OrderedDict.fromkeys((inst, syst, bset)).keys())))
    nice_name = band + ((' [' + extra_nice + ']') if extra_nice else '')

    realizations = [[] for x in range(len(model['realizations']))]
    for ph in photo:
        rn = ph.get('realization', None)
        si = ph.get('simulated', False)
        if rn and not si:
            if tuple(ph.get(y, '') for y in band_attr) == full_band:
                realizations[int(rn) - 1].append((
                    float(ph['time']), float(ph['magnitude']), [
                        float(ph.get('e_lower_magnitude', ph.get('e_magnitude', 0.0))),
                        float(ph.get('e_upper_magnitude', ph.get('e_magnitude', 0.0)))],
                ph.get('upperlimit')))
    numrz = np.sum([1 for x in realizations if len(x)])
    for rz in realizations:
        if not len(rz):
            continue
        xs, ys, vs, us = zip(*rz)
        label = '' if full_band in used_bands or full_band in real_band_list else nice_name
        if max(vs) == 0.0:
            plt.plot(xs, ys, color=bandcolorf(band),
                             label=label, linewidth=0.5)
        else:
            xs = np.array(xs)
            ymi = np.array(ys) - np.array([np.inf if u else v[0] for v, u in zip(vs, us)])
            yma = np.array(ys) + np.array([v[1] for v in vs])
            plt.fill_between(xs, ymi, yma, color=bandcolorf(band), edgecolor=None,
                             label=label, alpha=1.0/numrz, linewidth=0.0)
            plt.plot(xs, ys, color=bandcolorf(band),
                             label=label, alpha=1.0, linewidth=0.5)
        if label:
            used_bands = list(set(used_bands + [full_band]))
    if real_data:
        for s in range(2):
            if s == 0:
                cond = False
                symb = 'o'
            else:
                cond = True
                symb = 'v'
            vec = [(float(x['time']), float(x['magnitude']),
                    0.0 if 'upperlimit' in x else float(x.get('e_lower_magnitude', x.get('e_magnitude', 0.0))),
                    float(x.get('e_upper_magnitude', x.get('e_magnitude', 0.0)))) for x in photo
                   if 'magnitude' in x and ('realization' not in x or 'simulated' in x) and
                   'host' not in x and 'includeshost' not in x and
                   x.get('upperlimit', False) == cond and
                   tuple(x.get(y, '') for y in band_attr) == full_band]
            if not len(vec):
                continue
            xs, ys, yls, yus = zip(*vec)
            label = nice_name if full_band not in used_bands else ''
            plt.errorbar(xs, ys, yerr=(yus, yls), color=bandcolorf(band), fmt=symb,
                         label=label,
                         markeredgecolor='black', markeredgewidth=1, capsize=1,
                         elinewidth=1.5, capthick=2, zorder=10)
            plt.errorbar(xs, ys, yerr=(yus, yls), color='k', fmt=symb, capsize=2,
                         elinewidth=2.5, capthick=3, zorder=5)
            if label:
                used_bands = list(set(used_bands + [full_band]))
plt.margins(0.02, 0.1)
plt.legend()
plt.show()
fig.savefig('../products/lc.pdf')



#==============================================================================

fig = plt.figure(figsize=(12,8))
ax = plt.gca()
used_instruments = []
color = {}
for i in xray_instrument_list:
    color[i[0]] =  next(ax._get_lines.prop_cycler)['color'] # using this to match MOSFiT color scheme

for all_instrument_attr in tqdm_notebook(xray_instrument_list, desc='Photo', leave=False):
    (inst, telscp) = all_instrument_attr
    try:
        inst_exclusive_list
    except:
        pass
    else:
        if inst not in inst_exclusive_list:
            continue
    extra_nice = ', '.join(list(filter(None, OrderedDict.fromkeys((inst, syst)).keys())))
    nice_name = telscp + (('[' + extra_nice + ']') if extra_nice else '')
    realizations = [[] for x in range(len(model['realizations']))]
    alias = ['' for x in range(len(model['realizations']))]
    for ph in photo:
        rn = ph.get('realization', None)
        if rn:
            if tuple(ph.get(y, '') for y in xray_instrument_attr) == all_instrument_attr:
                    realizations[int(rn) - 1].append((float(ph['time']),
                                                      float(ph['countrate']),
                                                      float(ph.get('e_lower_countrate', 0.0)),
                                                      float(ph.get('e_upper_countrate', 0.0))))
                    alias[int(rn) - 1] = ph['realization']
    for i,rz in enumerate(realizations):
        if not len(rz):
            continue
        xs, ys, vls, vus = zip(*rz)
        label = ('' if all_instrument_attr in used_instruments or
                 all_instrument_attr in real_xray_instrument_list else nice_name)
        if max(vs) == 0.0:
            plt.plot(xs, ys, color=color[inst],
                             label=label, linewidth=0.5)
        else:
            xs = np.array(xs)
            ymi = np.array(ys) - np.array(vls)
            yma = np.array(ys) + np.array(vus)
            plt.fill_between(xs, ymi, yma, color=color[inst],
                             label=label, alpha=2.0/len(realizations))

        if label:
            used_instruments = list(set(used_instruments + [all_instrument_attr]))
    if real_data:
        for s in range(2):
            if s == 0:
                cond = False
                symb = 'o'
            else:
                cond = True
                symb = 'v'
            vec = [(float(x['time'][0]), float(x['countrate']), float(x.get('e_countrate', 0.0))) for x in photo
                   if 'countrate' in x and 'realization' not in x and
                   'host' not in x and 'includeshost' not in x and
                   x.get('upperlimit', False) == cond and
                   tuple(x.get(y, '') for y in xray_instrument_attr) == all_instrument_attr]

            if not len(vec):
                continue
            xs, ys, yes = zip(*vec)
            label = nice_name if all_instrument_attr not in used_instruments else ''

            plt.errorbar(xs, ys, yerr=yes, color = color[inst], fmt=symb,
                         label=label,
                         markeredgecolor='black', markeredgewidth=1, capsize=5,
                         elinewidth=2, capthick=2, zorder=10)
            plt.errorbar(xs, ys, yerr=yes, color='k', fmt=symb, capsize=6,
                         elinewidth=3, capthick=3, zorder=5)
            if label:
                used_instruments = list(set(used_instruments + [all_instrument_attr]))
plt.gca().set_xlabel('MJD')
plt.gca().set_ylabel('countrates')
#plt.xlim(55200, 56000)
plt.yscale('log')
#plt.ylim(1e-320,1e20)
plt.margins(0.1,0.1)
plt.legend(loc = 'best')
plt.show()
fig.savefig('../products/lc_xrays.pdf')


#==============================================================================

with open('../products/chain.json', 'r', encoding = 'utf-8') as f:
    all_chain = np.array(json.load(f))

nparam = len(all_chain[0,0,0,:]);
fig = plt.figure(figsize=(4. * np.ceil(nparam / 4.), 8));
for pi, param in enumerate(range(nparam)):
    my_chain = all_chain[0, :, :, param]
    ax = fig.add_subplot(np.ceil(nparam / 4.), 4, pi + 1);
    ax.plot(my_chain.T);
    ax.plot(np.mean(my_chain, axis=0), color='k');
plt.tight_layout()


#==============================================================================
import logging
logging.disable(logging.WARNING)

# Construct walker arrays for corner
corner_input = []
pars = [x for x in model['setup'] if model['setup'][x].get('kind') == 'parameter' and
        'min_value' in model['setup'][x] and 'max_value' in model['setup'][x]]
weights = []
for realization in model['realizations']:
    par_vals = realization['parameters']
    if 'weight' in realization:
        weights.append(float(realization['weight']))
    var_names = ['$' + ('\\log\\, ' if par_vals[x].get('log') else '') +
                 par_vals[x]['latex'] + '$' for x in par_vals if x in pars and 'fraction' in par_vals[x]]
    corner_input.append([np.log10(par_vals[x]['value']) if
                         par_vals[x].get('log') else par_vals[x]['value'] for x in par_vals
                         if x in pars and 'fraction' in par_vals[x]])
weights = weights if len(weights) else None
ranges = [0.999 for x in range(len(corner_input[0]))]
cfig = corner.corner(corner_input, labels=var_names, quantiles=[0.16, 0.5, 0.84],
                     show_titles=True, weights=weights, range=ranges)
cfig.savefig('../products/corner.pdf')

#===============================================================================
