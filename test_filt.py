import xarray as xr
import numpy as np
import scipy as sp
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import os
from dask.distributed import Client
from dask.diagnostics import ProgressBar
import seapy
import dask.multiprocessing
dask.config.set(schedular='multiprocessing')

dt = 1
fs = 1/(dt)
Nf = fs/2

Cf = fs/36
M = 100

def lanczos_filter(data,Cf,Nf,M):
    coef = sp.signal.firwin(M+1, Cf/Nf,width = 2/len(data),window='lanczos', pass_zero='lowpass')
    return sp.signal.filtfilt(coef,1.0,data)


Inp_files = ['/Volumes/WD_3/outputs_SCORRECTION/ocean_ecs_avg_'\
             +str(i).rjust(4,'0')+'.nc' for i in range(1,446)]
# Inp_files = ['E:\\roms4dvar_ecs\\freerun_detide_v4.1\\filter\\ocean_ecs_avg_'\
             # +str(i).rjust(4,'0')+'.nc' for i in range(1,30)]
Inp_files

Inp_ds = xr.open_mfdataset(Inp_files,combine='nested',\
                            coords='minimal',compat='override',concat_dim='ocean_time',\
                            data_vars=['temp','salt','u','v','zeta'],parallel=True,\
                           chunks={'eta_rho':181,'xi_rho':121,'s_rho':1,\
                                  'eta_u':181,'xi_u':121,\
                                  'eta_v':181,'xi_v':121,\
                                  'eta_psi':181,'xi_psi':121,})


Temp_da = Inp_ds.temp
print(Temp_da)


filt_temp = xr.apply_ufunc(lanczos_filter,Temp_da,Cf,Nf,M,
               input_core_dims=[["ocean_time"],[],[],[]],
               output_core_dims=[["ocean_time"]],
               vectorize=True,
               dask="allowed",
               output_dtypes=[Temp_da.dtype]
              )
print(filt_temp)

