#!/usr/bin/env python
# coding: utf-8

# In[48]:


import os
os.chdir('/data7/cschoi/sngal/NGC3367/phot/cut/fitting/sncosmo')
os.getcwd()
os.listdir()


# In[49]:


import astropy.io.ascii as ascii
import numpy as np
import sys
import sncosmo as sc
import sncosmo


# In[3]:


gax=ascii.read('SN2017gax_griz_vstack.dat')


# In[4]:


gax1=gax[gax['MAG']!=-99]


# In[5]:


ab = sncosmo.get_magsystem('ab')


# In[6]:


ab.zpbandflux('sdssg')
ab.zpbandflux('sdssr')
ab.zpbandflux('sdssi')
ab.zpbandflux('sdssz')


# In[10]:


#import what_the_flux as wtf


# In[11]:


gmag=gax1['MAG'][gax1['FILTER']=='g']
gmagerr=gax1['MAGERR'][gax1['FILTER']=='g']


# In[12]:


gmagflux=ab.band_mag_to_flux(gmag, 'sdssg')


# In[13]:


gmagfluxerr=gmagerr* np.log(10)/2.5*gmagflux


# In[14]:


gmagflux[0],gmagflux[0]/5.466008340859812e-05


# In[15]:


gmagfluxerr[0],gmagfluxerr[0]/5.466008340859812e-05


# flux 1.0= 10 **(-25/2.5) * ab_zp_flux(at badpass)

# In[16]:


scale_g=10**(-10)*ab.zpbandflux('sdssg')


# In[17]:


#gmagflux=gmagflux/scale_g
#gmagfluxerr=gmagfluxerr/scale_g


# In[18]:


scale_r=10**(-10)*ab.zpbandflux('sdssr')


# In[19]:


scale_i=10**(-10)*ab.zpbandflux('sdssi')


# In[20]:


scale_z=10**(-10)*ab.zpbandflux('sdssz')


# In[21]:


rmag=gax1[gax1['FILTER']=='r']['MAG']
rmagerr=gax1[gax1['FILTER']=='r']['MAGERR']
imag=gax1[gax1['FILTER']=='i']['MAG']
imagerr=gax1[gax1['FILTER']=='i']['MAGERR']
zmag=gax1[gax1['FILTER']=='z']['MAG']
zmagerr=gax1[gax1['FILTER']=='z']['MAGERR']


# In[22]:


rmagflux=ab.band_mag_to_flux(rmag, 'sdssr')
rmagfluxerr=rmagerr* np.log(10)/2.5*rmagflux
#rmagflux=rmagflux/scale_r
#rmagfluxerr=rmagfluxerr/scale_r


# In[23]:


imagflux=ab.band_mag_to_flux(imag, 'sdssi')
imagfluxerr=imagerr* np.log(10)/2.5*imagflux
#imagflux=imagflux/scale_i
#imagfluxerr=imagfluxerr/scale_i


# In[24]:


zmagflux=ab.band_mag_to_flux(zmag, 'sdssz')
zmagfluxerr=zmagerr* np.log(10)/2.5*zmagflux
#zmagflux=zmagflux/scale_z
#zmagfluxerr=zmagfluxerr/scale_z


# In[25]:


from astropy.table import Table, vstack


# In[26]:


t=Table()
dd=gax1[gax1['FILTER']=='g']
band=['sdssg']*len(dd)
zp=[25]*len(dd)
zpsys=['ab']*len(dd)
t.add_column(dd['MJD'],name='time')
t.add_column(band,name='band')
t.add_column(gmagflux,name='flux')
t.add_column(gmagfluxerr,name='fluxerr')
t.add_column(zp,name='zp')
t.add_column(zpsys,name='zpsys')
tg=t


# In[27]:


t=Table()
dd=gax1[gax1['FILTER']=='r']
band=['sdssr']*len(dd)
zp=[25]*len(dd)
zpsys=['ab']*len(dd)
t.add_column(dd['MJD'],name='time')
t.add_column(band,name='band')
t.add_column(rmagflux,name='flux')
t.add_column(rmagfluxerr,name='fluxerr')
t.add_column(zp,name='zp')
t.add_column(zpsys,name='zpsys')
tr=t


# In[28]:


t=Table()
dd=gax1[gax1['FILTER']=='i']
band=['sdssi']*len(dd)
zp=[25]*len(dd)
zpsys=['ab']*len(dd)
t.add_column(dd['MJD'],name='time')
t.add_column(band,name='band')
t.add_column(imagflux,name='flux')
t.add_column(imagfluxerr,name='fluxerr')
t.add_column(zp,name='zp')
t.add_column(zpsys,name='zpsys')
ti=t


# In[29]:


t=Table()
dd=gax1[gax1['FILTER']=='z']
band=['sdssz']*len(dd)
zp=[25]*len(dd)
zpsys=['ab']*len(dd)
t.add_column(dd['MJD'],name='time')
t.add_column(band,name='band')
t.add_column(zmagflux,name='flux')
t.add_column(zmagfluxerr,name='fluxerr')
t.add_column(zp,name='zp')
t.add_column(zpsys,name='zpsys')
tz=t


# In[30]:


#t.write('SN2017ein-FLI-B.dat',format='ascii.commented_header',overwrite=True)


# In[31]:


t=vstack([tg,tr,ti])
t.sort('time')#,'OBS')

t.write('SN2017gax-sncosmo.dat',format='ascii.commented_header',overwrite=True)


# In[32]:


ls


# In[33]:


data=t
data


# In[34]:


dust = sncosmo.CCM89Dust()
model = sncosmo.Model(source='v19-2007gr-corr',
                         effects=[dust,dust],
                         effect_names=['host','mw'],
                         effect_frames=['rest','obs'])
#‘nugent-sn1bc’
#‘v19-2007gr-corr’
#‘v19-2007gr’
model.param_names

#model.parameters


# In[35]:


model.set(z=0.004464)
model.set(mwr_v=3.1)
model.set(mwebv=0.021)
model.set(hostr_v=3.1)


# In[36]:


result, fitted_model = sncosmo.fit_lc(
    data, model,
    ['t0', 'amplitude', 'hostebv'], 
    bounds={'z':(0., 0.005)}
)
print('fitting is done')
print("Number of chi^2 function calls:", result.ncall)
print("Number of degrees of freedom in fit:", result.ndof)
print("chi^2 value at minimum:", result.chisq)
print("model parameters:", result.param_names)
print("best-fit values:", result.parameters)
print("The result contains the following attributes:\n", result.keys())


# In[37]:


result


# In[38]:


sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
sncosmo.plot_lc(data, model=fitted_model, errors=result.errors,fname='SN2017gax_sncosmo.png')


# In[ ]:





# In[ ]:




