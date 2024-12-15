from functools import partial
import xarray as xr
import numpy as np
import scipy as sp
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from lanczos_filt import *
import seapy
from numba import float64, guvectorize
import dask
from dask.diagnostics import ProgressBar

#import dask.multiprocessing
pid=11
#@guvectorize(
#        "(float64[:], float64[:], float64[:], float64[:], float64[:])",
#        "(m), (n),(n),(n) -> (m)",
#        forceobj=True,        
#        )
#def lfilter_np_gufunc(data, cf,Nf,M, out):
#    out[:] = lanczos_filter(data,cf,Nf,M)
#
#
#def xr_lfilter(data,Cf,Nf,M):
#    return xr.apply_ufunc(lfilter_np_gufunc, data,Cf,Nf,M,
#                           input_core_dims=[["ocean_time"],[],[],[]],
#                           output_core_dims=[["ocean_time"]],
#                           #vectorize=True,
#                           dask="parallelized",
#                           output_dtypes=[Var_da.dtype]
#                          )
#
#

#dask.config.set(schedular='multiprocessing')
#dask.config.set(schedular='threads', num_workers=4)
#from dask.distributed import Client, LocalCluster
# client = cluster.get_client()
#client = Client(threads_pre_worker=2, n_worker=8)
Inp_files = ['/Volumes/WD_3/outputs_SCORRECTION/ocean_ecs_avg_'\
              +str(i).rjust(4,'0')+'.nc' for i in range(1,30)]
#Inp_files = ['E:\\outputs_SCORRECTION\\ocean_ecs_avg_'\
#             +str(i).rjust(4,'0')+'.nc' for i in range(180,447)]
print(len(Inp_files))
fields = set(seapy.roms.fields)

grid = seapy.model.asgrid(Inp_files[0])    
nc = seapy.netcdf(Inp_files[0])

fields = set(nc.variables).intersection(fields)

epoch, time_var = seapy.roms.get_reftime(nc)

def _preprocess(x, lon_inds, lat_inds):
    return x.isel(xi_rho=slice(*lon_inds),eta_rho=slice(*lat_inds),
                   xi_u=slice(*lon_inds),eta_u=slice(*lat_inds),
                 xi_v=slice(*lon_inds),eta_v=slice(*lat_inds),
                 xi_psi=slice(*lon_inds),eta_psi=slice(*lat_inds),)


lon_inds, lat_inds = (0,49), (0,61)
partial_func = partial(_preprocess, lon_inds=lon_inds, lat_inds=lat_inds)

Inp_ds = xr.open_mfdataset(Inp_files,combine='nested',\
                            coords='minimal',compat='override',concat_dim='ocean_time',\
                            preprocess=partial_func,
                            drop_variables = {'h','f','pm','pn','omega','w','ssflux','shflux','sustr','svstr'},
                            data_vars=['temp','salt','u','v','zeta'],parallel=True,\

                            chunks={'eta_rho':31,'xi_rho':25,'s_rho':20,\
                                  'eta_u':31,'xi_u':25,\
                                  'eta_v':31,'xi_v':25,\
                                  'eta_psi':31,'xi_psi':25,}).chunk(dict(ocean_time=-1))

#with ProgressBar():
#    test1d = Temp_ds.isel(s_rho=-1,xi_rho=172,eta_rho=66).data.compute()

# filter data using matlab lanczos filter
#dt = 1
#fs = 1/(dt)
#Nf = fs/2
#
#Cff = fs/18
#M = 50
#N = len(test1d)
#coeff = lanczos_filter_coef(Cff,M)
#coeff = coeff[:,0]
#window, Ff = spectral_window(coeff,N)
## print(window.shape)
#y,Cx = spectral_filtering(test1d,window);

#filter data using scipy lanczos filter
dt = 1
fs = 1/(dt)
Nf = fs/2

Cf = fs/36
M = 100

#with ProgressBar():
#    filter_temp1d = filt_var.isel(s_rho=-1,xi_rho=172,eta_rho=66).data.compute()

#plt.plot(test1d[:],label='raw')
## plt.plot(filter_test1d[99:99+1000])
#plt.plot(filter_temp1d[:],label='xr_lanczos_filt')
#plt.plot(y[:],label='mat_lanczos_filt')
#plt.legend()
#plt.show()

#Var_da = Inp_ds.temp
#print(Var_da)
#
#filt_var_ds = xr.apply_ufunc(lanczos_filter,Var_da,Cf,Nf,M,
#               input_core_dims=[["ocean_time"],[],[],[]],
#               output_core_dims=[["ocean_time"]],
#               vectorize=True,
#               dask="parallelized",
#               output_dtypes=[Var_da.dtype]
#              )
#
#
#
#
#filt_var_mongrp = filt_var_ds.groupby("ocean_time.month")
#print(filt_var_mongrp.groups.keys())
#time = np.array([float(i)*1e-9 for i in filt_var_mongrp[1].ocean_time.data])
#
#filt_var_monstd = filt_var_mongrp.std("ocean_time")
#print(filt_var_monstd)
filt_var = {}
with ProgressBar():
    for month in [1]:
        print('generating std file in month %s ..'%(str(month)))
#        ncout = seapy.roms.ncgen.create_da_ini_std('F:\\roms4dvar_ecs\\i4dvar_outputs\\bg_std\\ecs_std_i_v4sc_%s.nc'%(str(month)),
            # If there are any fields that are not in the standard output file,
            # add them to the output file
#        for f in fields.difference(ncout.variables):
#            ncout.createVariable(f, np.float32,
#                                 ('ocean_time',"s_rho", "eta_rho", "xi_rho"))



        for v in ['zeta']:
            print(v)
            Var_da = Inp_ds[v]
            print(Var_da)

#            filt_var_ds = xr_lfilter(Var_da,Cf,Nf,M)

            filt_var_ds = xr.apply_ufunc(lanczos_filter,Var_da,Cf,Nf,M,
                           input_core_dims=[["ocean_time"],[],[],[]],
                           output_core_dims=[["ocean_time"]],
                           vectorize=True,
                           dask="parallelized",
                           output_dtypes=[Var_da.dtype]
                          )

            filt_var_mongrp = filt_var_ds.groupby("ocean_time.month")
            #filt_var_yeargrp = filt_var_ds.groupby("ocean_time.year")
            #filt_var
    #        print(filt_var_mongrp.groups.keys())
            time = np.array([float(i)*1e-9 for i in filt_var_mongrp[month].ocean_time.data])

            filt_var_monstd = filt_var_mongrp.std("ocean_time")
#            filt_dat = np.zeros(filt_var_monstd.sel(month=month).shape)
#            print(filt_var_monstd)
            if v in ['u','v','temp','salt']: 
#                for iz in range(0,20):
#                    filt_dat[i,:,:] = filt_var_monstd.sel(month=month).isel(s_rho=iz).data.compute()
                filt_dat = filt_var_monstd.sel(month=month).data.compute()
                filt_var[v] = filt_dat.copy()

            else:
                filt_dat = filt_var_monstd.sel(month=month).data.compute()
                filt_var[v] = filt_dat.copy()
    sp.io.savemat('ecs_std_part%s.mat'%(pid),filt_var)


#    for month in filt_var_mongrp.groups.keys():
#        filt_var['mon_'+str(1)] = filt_var_monstd[month].data.compute()
#        filt_var_monstd.sel(s_rho=-0.025,month=month).plot()

#import scipy.io as sio
#sio.savemat('filt_var_std.mat', filt_var)
