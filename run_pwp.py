import PWP_helper as phf
import PWP

def run_niw(suffix):
    '''
    This is a wrapper to run PWP on EM-APEX float data in the Western Pacific
    '''
    forcing_fname = 'niw_7788a_met.nc'
    prof_fname = 'niw_7788a_profile.nc'
    print("Running Test Case with data from EM-APEX floats...")

    p = {}
    p['dt'] = 1.
    p['dz'] = 2.
    p['mld_thresh'] =1e-4
    p['rkz'] = 0
    p['emp_ON'] = False

    forcing, pwp_out = PWP.run(met_data=forcing_fname,
                               prof_data=prof_fname,
                               suffix=suffix,
                               save_plots=True,
                               diagnostics=False,
                               param_kwds=p)

run_niw('test00')

# %%
import xarray as xr

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from cmocean import cm

# set up figure params
sns.set(style='ticks',context='paper')
mpl.rc('figure', dpi=100, figsize=[8.5,11])
mpl.rc('savefig',dpi=500,bbox='tight')
mpl.rc('legend',frameon=False)

# %%
results = xr.open_dataset('./output/pwp_output_test00.nc')
results = results.assign_coords(z=(-1)*results.z)

forcing = xr.open_dataset('./input_data/niw_7788a_met.nc')
forcing['qnet'] = -forcing.qlat-forcing.qsens+forcing.sw+forcing.lw

results['hke'] = 0.5*(results.uvel**2+results.vvel**2)
# %%
def compute_mld(data):
	'''
	mld criterion: depth of z|_(rho_10m +0.03kg/m3)

	'''
	import numpy as np
	import xarray as xr
	from scipy.interpolate import interp1d

	mld = np.zeros(data.time.size)
	for t in range(data.time.size):
		test = data.isel(time=t)
		# test = test.dropna(dim='z')
		if test.dens.size>0:
			f = interp1d(test.dens,test.z)
			a = test.dens.interp(z=-10)+0.03
			a
			if (a>test.dens.min()) & (a<test.dens.max()):
				mld[t] = f(a)
			else:
				mld[t]=np.nan
		else:
			mld[t]=np.nan
	data['mld'] = ('time',mld)
	return data


results=compute_mld(results)
results['dens']=results.dens-1000

# %%
from mpl_toolkits.axes_grid1 import make_axes_locatable
zmin=-80
tmin=1
tmax=10
f,ax=plt.subplots(6,1,sharex=True)
forcing.tx.sel(time=slice(tmin,tmax)).plot(ax=ax[0],label=r'$\tau_x$')
forcing.ty.sel(time=slice(tmin,tmax)).plot(ax=ax[0],label=r'$\tau_y$')
ax[0].legend()
ax[0].set_xlabel(None)
ax[0].set_ylabel('z[m] ')
ax[0].set_ylabel(r'$\tau$ [Nm$^{-2}$]')

divider = make_axes_locatable(ax[0])
cax= divider.append_axes("right", size="18%", pad=.05)
cax.remove()

forcing.qnet.sel(time=slice(tmin,tmax)).plot(ax=ax[1],label=r'$\tau_y$')
ax[1].axhline(0)
divider = make_axes_locatable(ax[1])
cax= divider.append_axes("right", size="18%", pad=.05)
cax.remove()
ax[1].set_xlabel(None)
ax[1].set_ylabel(r'Q$_{net}$ [Wm$^{-2}$]')

results.dens.sel(time=slice(tmin,tmax)).plot(y='z',ax=ax[2],ylim=(zmin,0),
                                        vmin=22.8,vmax=25,
                                        cbar_kwargs={'pad':0.01,'label':r'sigma$_0$ [kgm$^{-3}$]'})
results.mld.plot(ax=ax[2],color='black')
ax[2].set_xlabel(None)
ax[2].set_ylabel('z [m]')

results.uvel.sel(time=slice(tmin,tmax)).plot(y='z',ax=ax[3],ylim=(zmin,0),
                                             cbar_kwargs={'pad':0.01,'label':r'u [ms$^{-1}$]'})
results.mld.plot(ax=ax[3],color='black',label='MLD')
ax[3].legend()
ax[3].set_xlabel(None)
ax[3].set_ylabel('z [m]')

results.vvel.sel(time=slice(tmin,tmax)).plot(y='z',ax=ax[4],ylim=(zmin,0),
                                             cbar_kwargs={'pad':0.01,'label':r'v [ms$^{-1}$]'})
results.mld.plot(ax=ax[4],color='black',label='MLD')
ax[4].legend()
ax[4].set_xlabel(None)
ax[4].set_ylabel('z [m]')

results.hke.sel(time=slice(tmin,tmax)).plot(y='z',ax=ax[5],ylim=(zmin,0),
                                            cbar_kwargs={'pad':0.01,'label':r'HKE [m$^2$s$^{-2}$]'})
results.mld.plot(ax=ax[5],color='white',label='MLD')
leg = ax[5].legend()
plt.setp(leg.get_texts(), color='w')
ax[5].set_xlabel('Time [days]')
ax[5].set_ylabel('z [m]')
# plt.tight_layout()
plt.savefig('../../figures/pwp_7788a.pdf')
plt.show()
