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
import sys
from concurrent import futures
import multiprocessing
#from multiprocessing import Lock
from multiprocessing.dummy import Pool 
#from multiprocessing import Pool
import ctypes
#import dask.multiprocessing
#dask.config.set(schedular='multiprocessing')
#dask.config.set(schedular='processes', num_workers=6)
dask.config.set(schedular='threads', num_workers=5)
#from dask.distributed import Client, LocalCluster
#client = Client()
#global Var_filt
#def init(shared_array_base,shape,v,Var_da,Cf,Nf,M):
#    global Var_filt
#    global sv,sVar_da, sCf,sNf,sM
#    sv = v
#    sVar_da = Var_da.copy()
#    sCf = Cf
#    sNf = Nf
#    sM = M
#    Var_filt = np.ctypeslib.as_array(shared_array_base.get_obj())
#    Var_filt = Var_filt.reshape(shape)
#

def mp_filt_var(sVar_da,ix,iy,sCf,sNf,sM):
    
#    print(iy,ix)
    if  ix ==0:
        print("processing %f%%..."%(iy/362*100))
    # extract data column's time series at this 
    # point(eta_rho(iy), xi_rho(ix))
    if v in ['temp', 'salt']:
        var_t_col = sVar_da.isel(eta_rho=iy,xi_rho=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,sCf,sNf,sM)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif v == 'u':
        var_t_col = sVar_da.isel(eta_u=ix,xi_u=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,sCf,sNf,sM)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif v == 'v':
        var_t_col = sVar_da.isel(eta_v=iy,xi_v=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,sCf,sNf,sM)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif v == 'zeta':
        var_t_col = sVar_da.isel(eta_rho=iy,xi_rho=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,0,var_t_col,sCf,sNf,sM)

        Var_filt[iy,ix,:] = filt_var_t_col.copy()

#        print(np.nanmean(filt_var_t_col))
#        print(filt_var_t_col.shape)



if __name__ == '__main__':
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

    global Cf, Nf, M, Var_da, Var_filt, v
    dt = 1
    fs = 1/(dt)
    Nf = fs/2

    Cf = fs/36
    M = 100

    v = 'zeta'
    print(v)
    Var_da = Inp_ds[v]
    #        global Var_da
    print(Var_da.shape)
    if v == 'zeta':
        (ly,lx,lt) = Var_da.shape
        shape = Var_da.shape
        lz = 1
    else:
        (ly,lx,lz, lt) = Var_da.shape 
        shape = Var_da.shape
    Var_filt = np.zeros(shape)    
#    iter_x = []
#    [[iter_x.append(ix) for ix in range(lx)] for iy in range(ly)]
#    iter_y = []
#    [[iter_y.append(iy) for ix in range(lx)] for iy in range(ly)]
#    shared_array_base = multiprocessing.Array(ctypes.c_double,ly*lx*lz*lt)
   #
    #for v in fields:
    #with futures.ThreadPoolExecutor(max_workers=8) as executor:
    #    for iy in range(ly):
    #        for ix in range(lx):
    #            executor.submit(mp_filt_var,Var_da,ix,iy,Cf,Nf,M)
    p = Pool(252,)
#    p.map(mp_filt_var, zip(iter_x,iter_y))
#    threads = []
#    count = 0
    for iy in range(ly):
        for ix in range(lx):
            p.apply_async(mp_filt_var, args=(Var_da,ix,iy,Cf,Nf,M))
#            t = p.apply_async(mp_filt_var,args =(Var_da,ix,iy,Cf,Nf,M))
    #        threads.append(t)
    #        count +=1
    p.close()
    p.join()
#    Var_filt = np.ctypeslib.as_array(shared_array_base.get_obj())
#    Var_filt = Var_filt.reshape(shape)
    
    print(Var_filt)#
    #    p.close()
    import matplotlib.pyplot as plt
    sp.io.savemat('testmt_filt_%s.mat'%(v),{'%s_std'%(v): Var_filt})    
    np.save('testmt_filt_%s'%(v),Var_filt)
