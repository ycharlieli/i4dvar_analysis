import xarray as xr
import numpy as np
import scipy as sp
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from lanczos_filt import *
import seapy
import dask
from dask.diagnostics import ProgressBar
#import dask.multiprocessing
#dask.config.set(schedular='multiprocessing')
#dask.config.set(schedular='processes', num_workers=8)
dask.config.set(schedular='threads', num_workers=8)
#from dask.distributed import Client, LocalCluster
#client = Client()
#from dask_mpi import initialize
#initialize()
#
#from distributed import Client
#client = Client()

Inp_files = ['/Volumes/TO_1/outputs_SCORRECTION/ocean_ecs_avg_'\
              +str(i).rjust(4,'0')+'_t.nc' for i in range(1,3)]
template_file = '/Volumes/WD_3/outputs_SCORRECTION/ocean_ecs_avg_0001.nc'
#Inp_files = ['E:\\outputs_SCORRECTION\\ocean_ecs_avg_'\
#             +str(i).rjust(4,'0')+'.nc' for i in range(180,447)]
print(len(Inp_files))
fields = set(seapy.roms.fields)

grid = seapy.model.asgrid(template_file)
nc = seapy.netcdf(template_file)

fields = set(nc.variables).intersection(fields)

epoch, time_var = seapy.roms.get_reftime(nc)



Inp_ds = xr.open_mfdataset(Inp_files,
                            coords='minimal',compat='override',
                            data_vars=['temp','salt','u','v','zeta'],parallel=True,
                           chunks={'eta_rho':181,'xi_rho':121,'s_rho':1,
                           'eta_u':181,'xi_u':121,\
                           'eta_v':181,'xi_v':121,\
                           'eta_psi':181,'xi_psi':121,},
                           engine = 'netcdf4' )

from dask.diagnostics import ProgressBar

dt = 1
fs = 1/(dt)
Nf = fs/2

Cf = fs/36
M = 100


with ProgressBar():
    for v in fields:
        print(v)
        Var_da = Inp_ds[v]
        print(Var_da.shape)
        (ly,lx,lz, lt) = Var_da.shape 
        Var_std = np.zeros((ly,lx,lz))
        for iy in range(10):
            for ix in range(10):
                # extract data column's time series at this 
                # point(eta_rho(iy), xi_rho(ix))
                if v in ['temp', 'salt']:
                    var_t_col = Var_da.isel(eta_rho=ix,xi_rho=iy).data.compute()
                    print(var_t_col.shape)
                    filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,Cf,Nf,M)
                elif v == 'u':
                    var_t_col = Var_da.isel(eta_u=ix,xi_u=iy).data.compute()
                    print(var_t_col.shape)
                    filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,Cf,Nf,M)
 
                elif v == 'v':
                    var_t_col = Var_da.isel(eta_v=ix,xi_v=iy).data.compute()
                    print(var_t_col.shape)
                    filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,Cf,Nf,M)
                elif v == 'zeta':
                    var_t_col = Var_da.isel(eta_rho=ix,xi_rho=iy).data.compute()
                    print(var_t_col.shape)
                    filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,Cf,Nf,M)
 
