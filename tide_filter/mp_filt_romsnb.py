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
from concurrent import futures

#import dask.multiprocessing
#dask.config.set(schedular='multiprocessing')
#dask.config.set(schedular='processes', num_workers=6)
dask.config.set(schedular='threads', num_workers=64)
#from dask.distributed import Client, LocalCluster
#client = Client()
#from dask_mpi import initialize
#initialize()
#
#from distributed import Client
#client = Client()
def mp_filt_var(Var_da,ix,iy,Cf,Nf,M):
#    print(iy,ix)
    global Var_filt
    if  ix ==0:
        print("processing %f%%..."%(iy/362*100))
    # extract data column's time series at this 
    # point(eta_rho(iy), xi_rho(ix))
    if v in ['temp', 'salt']:
        var_t_col = Var_da.isel(eta_rho=iy,xi_rho=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,Cf,Nf,M)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif v == 'u':
        var_t_col = Var_da.isel(eta_u=ix,xi_u=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,Cf,Nf,M)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif v == 'v':
        var_t_col = Var_da.isel(eta_v=iy,xi_v=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,Cf,Nf,M)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif v == 'zeta':
        var_t_col = Var_da.isel(eta_rho=iy,xi_rho=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,0,var_t_col,Cf,Nf,M)

#        print(np.nanmean(filt_var_t_col))
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:] = filt_var_t_col.copy()




#for v in fields:
Inp_files = ['/Volumes/TO_1/outputs_SCORRECTION/ocean_ecs_avg_'\
              +str(i).rjust(4,'0')+'_t.nc' for i in range(1,30)]
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

v = 'zeta'
print(v)
Var_da = Inp_ds[v]
#        global Var_da
#    print(Var_da.shape)
if v == 'zeta':
    (ly,lx,lz) = Var_da.shape
    Var_filt = np.zeros((ly,lx,lz))
else:
    (ly,lx,lz, lt) = Var_da.shape 
    Var_filt = np.zeros((ly,lx,lz,lt))


with futures.ThreadPoolExecutor(max_workers=64) as executor:
#    with futures.ThreadPoolExecutor(max_workers=8) as executor:
    for iy in range(ly):
        for ix in range(lx):
            executor.submit(mp_filt_var, Var_da,ix,iy,Cf,Nf,M)
import matplotlib.pyplot as plt
sp.io.savemat('test_filt_%s.mat'%(v),{'%s_std'%(v): Var_filt})    
np.save('test_filt_%s'%(v),Var_filt)
