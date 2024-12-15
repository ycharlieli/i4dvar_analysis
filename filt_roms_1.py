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
dask.config.set(schedular='processes', num_workers=4)
#dask.config.set(schedular='threads', num_workers=48)
#from dask.distributed import Client, LocalCluster
#client = Client()
#from dask_mpi import initialize
#initialize()
#
#from distributed import Client
#client = Client()

Inp_files = ['/Volumes/WD_3/outputs_SCORRECTION/ocean_ecs_avg_'\
              +str(i).rjust(4,'0')+'.nc' for i in range(1,30)]
#Inp_files = ['E:\\outputs_SCORRECTION\\ocean_ecs_avg_'\
#             +str(i).rjust(4,'0')+'.nc' for i in range(180,447)]
print(len(Inp_files))
fields = set(seapy.roms.fields)

grid = seapy.model.asgrid(Inp_files)
nc = seapy.netcdf(Inp_files)

fields = set(nc.variables).intersection(fields)

epoch, time_var = seapy.roms.get_reftime(nc)



Inp_ds = xr.open_mfdataset(Inp_files,combine='nested',\
                            coords='minimal',compat='override',concat_dim='ocean_time',\
                            data_vars=['temp','salt','u','v','zeta'],parallel=True,\
                            chunks={'eta_rho':181,'xi_rho':61,'s_rho':10,\
                                  'eta_u':181,'xi_u':61,\
                                  'eta_v':181,'xi_v':61,\
                                  'eta_psi':181,'xi_psi':61,}).chunk(dict(ocean_time=-1))

from dask.diagnostics import ProgressBar
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
#filt_var = {}
with ProgressBar():
    for month in range(1,13):
        print('generating std file in month %s ..'%(str(month)))
#        ncout = seapy.roms.ncgen.create_da_ini_std('F:\\roms4dvar_ecs\\i4dvar_outputs\\bg_std\\ecs_std_i_v4sc_%s.nc'%(str(month)),
        ncout = seapy.roms.ncgen.create_da_ini_std('./bg_std/ecs_std_i_v4sc_%s.nc'%(str(month)),
                                                  eta_rho=grid.ln, xi_rho=grid.lm,
                                                  s_rho=grid.n,
                                                  reftime=epoch,
                                                  title='std from roms_ecsv4.1_scorrection')
        print(ncout.variables)
        grid.to_netcdf(ncout)
            # If there are any fields that are not in the standard output file,
            # add them to the output file
#        for f in fields.difference(ncout.variables):
#            ncout.createVariable(f, np.float32,
#                                 ('ocean_time',"s_rho", "eta_rho", "xi_rho"))
        depth = ncout.createVariable('s_rho',np.float32,('s_rho',))
        depth.long_name =  "S-coordinate at RHO-points"

        ncout['s_rho'][:] = grid.s_rho.copy()


        for v in fields:
            print(v)
            Var_da = Inp_ds[v]
#            print(Var_da)

            filt_var_ds = xr.apply_ufunc(lanczos_filter,Var_da,Cf,Nf,M,
                           input_core_dims=[["ocean_time"],[],[],[]],
                           output_core_dims=[["ocean_time"]],
                           vectorize=True,
                           dask="parallelized",
                          dask_gufunc_kwargs={'allow_rechunk',True},
                           output_dtypes=[Var_da.dtype]
                          )


            filt_var_ds = filt_var_ds.persist()

            filt_var_mongrp = filt_var_ds.groupby("ocean_time.month")
            #filt_var_yeargrp = filt_var_ds.groupby("ocean_time.year")
            #filt_var
    #        print(filt_var_mongrp.groups.keys())
            time = np.array([float(i)*1e-9 for i in filt_var_mongrp[month].ocean_time.data])

            filt_var_monstd = filt_var_mongrp.std("ocean_time")
#            print(filt_var_monstd)
            if v in ['u','v','temp','salt']:
                (lz,ly,lx) = filt_var_monstd.sel(month=month).shape

                for iz in range(lz):
                    for iy in range(ly):
                        for ix in range(lx):
                            if v in ['temp','salt']:

                                filt_dat = filt_var_monstd.sel(month=month).isel(s_rho=iz,
                                                                                 eta_rho = iy,
                                                                                 xi_rho = ix)
                                                                                 .data.compute()
                            elif v == 'u':
                                filt_dat = filt_var_monstd.sel(month=month).isel(s_rho=iz,
                                                                                 eta_u = iy,
                                                                                 xi_u = ix)
                                                                                 .data.compute()
 
                            elif v == 'v':
                                filt_dat = filt_var_monstd.sel(month=month).isel(s_rho=iz,
                                                                                 eta_v = iy,
                                                                                 xi_v = ix)
                                                                                 .data.compute()
                            nc.variables[v][0,iz,iy,ix]= filt_dat
 
#                filt_dat = filt_var_monstd.sel(month=month).data.compute()
#                ncout.variables[v][0,:] = filt_dat

            else:
                (ly,lx) = filt_var_monstd.sel(month=month).shape
                for iy in range(ly):
                    for ix in range(lx):
                        filt_dat = filt_var_monstd.sel(month=month).isel(
                                                                                 eta_rho = iy,
                                                                                 xi_rho = ix)
                                                                                 .data.compute()
                        nc.variables[v][0,iy,ix]=filt_dat
 
#                filt_dat = filt_var_monstd.sel(month=month).data.compute()
#                ncout.variables[v][0,:] = filt_dat
            ncout.variables[time_var] = np.mean(time)
            ncout.sync()

    ncout.close()
nc.close()

#    for month in filt_var_mongrp.groups.keys():
#        filt_var['mon_'+str(1)] = filt_var_monstd[month].data.compute()
#        filt_var_monstd.sel(s_rho=-0.025,month=month).plot()

#import scipy.io as sio
#sio.savemat('filt_var_std.mat', filt_var)
