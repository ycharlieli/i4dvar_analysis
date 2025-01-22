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
#from multiprocessing.dummy import Pool 
from multiprocessing import Pool
import ctypes
#import dask.multiprocessing
#dask.config.set(schedular='multiprocessing',num_workers=8)
#dask.config.set(schedular='processes', num_workers=6)

#from dask.distributed import Client, LocalCluster
#client = Client()
#global Var_filt
def init(shared_array_base,shape,v,Var_da,Cf,Nf,M):
    global Var_filt
    global sv,sVar_da, sCf,sNf,sM
    sv = v
    sVar_da = Var_da.copy()
    sCf = Cf
    sNf = Nf
    sM = M
    Var_filt = np.ctypeslib.as_array(shared_array_base.get_obj())
    Var_filt = Var_filt.reshape(shape)


def mp_filt_var(iterer):
    (ix,iy) = iterer
    
#    print(iy,ix)
#    if  ix ==0:
#        print("processing %f%%..."%(iy/362*100))
    # extract data column's time series at this 
    # point(eta_rho(iy), xi_rho(ix))
    if sv in ['temp', 'salt',]:
        var_t_col = sVar_da.isel(eta_rho=iy,xi_rho=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,sCf,sNf,sM)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif sv == 'u':
        var_t_col = sVar_da.isel(eta_u=ix,xi_u=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,sCf,sNf,sM)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif sv == 'v':
        var_t_col = sVar_da.isel(eta_v=iy,xi_v=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,1,var_t_col,sCf,sNf,sM)
#        print(filt_var_t_col.shape)
        Var_filt[iy,ix,:,:] = filt_var_t_col.copy()
    elif sv == 'zeta':
        var_t_col = sVar_da.isel(eta_rho=iy,xi_rho=ix).data.compute()
#        print(var_t_col.shape)
        filt_var_t_col = np.apply_along_axis(lanczos_filter,0,var_t_col,sCf,sNf,sM)

        Var_filt[iy,ix,:] = filt_var_t_col.copy()

#        print(np.nanmean(filt_var_t_col))
#        print(filt_var_t_col.shape)



if __name__ == '__main__':
    Tcount = {} 
    Var_std = {}
    Var_mean = {}
   
    var_list = ['zeta','u','v','temp','salt']
    for i in range(1,13):
        Tcount[i] = {}
        Var_std[i] = {}
        Var_mean[i] = {}
        for iv in var_list:
            Tcount[i][iv] = 0
 
    Find = 172
    Fstep = 5 # file number depends on the machine IO ability. Here 5 is optimal , beyond this  IO ability will contraint CPU efficiency
    for i in range(3):
        Inp_files = ['/Volumes/TO_1/outputs_SCORRECTION/ocean_ecs_avg_'\
                       +str(i).rjust(4,'0')+'_t.nc' for i in range(Find,Find+Fstep)]
        template_file = '/Volumes/TO_1/outputs_SCORRECTION/ocean_ecs_avg_0001.nc'
#        Inp_files = ['F:\\outputs_SCORRECTION\\ocean_ecs_avg_'\
#                      +str(i).rjust(4,'0')+'_t.nc' for i in range(Find,Find+Fstep)]
#        template_file = 'E:\\outputs_SCORRECTION\\ocean_ecs_avg_0001.nc'
        #Inp_files = ['E:\\outputs_SCORRECTION\\ocean_ecs_avg_'\
        #             +str(i).rjust(4,'0')+'.nc' for i in range(180,447)]
        #print(len(Inp_files))
        print('loading ')
        [print(ifile) for ifile in Inp_files]
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

        global Cf, Nf, M, Var_da, v
        dt = 1
        fs = 1/(dt)
        Nf = fs/2

        Cf = fs/36
        M = 100
        for v in var_list:
            print(v)
            Var_da = Inp_ds[v]
            #        global Var_da
            # print(Var_da.shape)
            if v == 'zeta':
                (ly,lx,lt) = Var_da.shape
                shape = Var_da.shape
                lz = 1
            else:
                (ly,lx,lz, lt) = Var_da.shape 
                shape = Var_da.shape
            
            iter_x = []
            [[iter_x.append(ix) for ix in range(lx)] for iy in range(ly)]
            iter_y = []
            [[iter_y.append(iy) for ix in range(lx)] for iy in range(ly)]
            shared_array_base = multiprocessing.Array(ctypes.c_double,ly*lx*lz*lt)
           #
            #for v in fields:
            #with futures.ThreadPoolExecutor(max_workers=8) as executor:
            #    for iy in range(ly):
            #        for ix in range(lx):
            #            executor.submit(mp_filt_var,Var_da,ix,iy,Cf,Nf,M)
            p = Pool(6,initializer=init, initargs=(shared_array_base,shape,v,Var_da,Cf,Nf,M))
            p.map(mp_filt_var, zip(iter_x,iter_y))
        #    threads = []
        #    count = 0
        #    for iy in range(ly):
        #        for ix in range(lx):
        #            p.apply_async(mp_filt_var, args=(Var_da,Var_filt,ix,iy,Cf,Nf,M))
        #            t = p.apply_async(mp_filt_var,args =(Var_da,ix,iy,Cf,Nf,M))
            #        threads.append(t)
            #        count +=1
            p.close()
            p.join()
            Var_filt = np.ctypeslib.as_array(shared_array_base.get_obj())
            Var_filt = Var_filt.reshape(shape)
             
            if v == 'zeta':
                Var_filt = Var_filt[:,:,48:-47]
                Filt_da = xr.DataArray(Var_filt,
                                       coords = {'eta_rho':Inp_ds.eta_rho.data,
                                                 'xi_rho': Inp_ds.xi_rho.data,
                                                 'ocean_time': Inp_ds.ocean_time.data[48:-47]}
                                       ).chunk({'eta_rho':181,'xi_rho':121})
            elif v in ['temp', 'salt']:
                Var_filt = Var_filt[:,:,:,48:-47]
                Filt_da = xr.DataArray(Var_filt,
                                       coords = {'eta_rho':Inp_ds.eta_rho.data,
                                                 'xi_rho': Inp_ds.xi_rho.data,
                                                 's_rho': Inp_ds.s_rho.data,
                                                 'ocean_time': Inp_ds.ocean_time.data[48:-47]}
                                       ).chunk({'eta_rho':181,'xi_rho':121,'s_rho':1})
            elif v == 'u':
                Var_filt = Var_filt[:,:,:,48:-47]
                Filt_da = xr.DataArray(Var_filt,
                                       coords = {'eta_u':Inp_ds.eta_u.data,
                                                 'xi_u': Inp_ds.xi_u.data,
                                                 's_rho': Inp_ds.s_rho.data,
                                                 'ocean_time': Inp_ds.ocean_time.data[48:-47]}
                                       ).chunk({'eta_u':181,'eta_u':121,'s_rho':1})
            elif v == 'v':
                Var_filt = Var_filt[:,:,:,48:-47]
                Filt_da = xr.DataArray(Var_filt,
                                       coords = {'eta_v':Inp_ds.eta_v.data,
                                                 'xi_v': Inp_ds.xi_v.data,
                                                 's_rho': Inp_ds.s_rho.data,
                                                 'ocean_time': Inp_ds.ocean_time.data[48:-47]}
                                       ).chunk({'eta_v':181,'eta_v':121,'s_rho':1})
            
            Filt_mongrp = Filt_da.groupby('ocean_time.month')
            for imonth in Filt_mongrp.groups.keys():
                print('month: ', imonth)
                print('accumulating data...')

                if Tcount[imonth][v] == 0:
#                    print('calculating std...')
#                    print(Filt_mongrp.groups[imonth])
                    Var_std[imonth][v]  = (Filt_mongrp[imonth]**2).sum(axis=-1)
#                        print('calculating mean...')
                    Var_mean[imonth][v]  = (Filt_mongrp[imonth]).sum(axis=-1)
                else:
#                    print('calculating std...')
                    Var_std[imonth][v] += (Filt_mongrp[imonth]**2).sum(axis=-1)
#                   print('calculating mean...')
                    Var_mean[imonth][v] += (Filt_mongrp[imonth]).sum(axis=-1)

                print('\tdata length of this month in this set: ', Filt_mongrp[imonth].ocean_time.shape[0])
                Tcount[imonth][v] +=  Filt_mongrp[imonth].ocean_time.shape[0]
                print(Filt_mongrp[imonth].ocean_time.shape[0])
        Find += Fstep
    for imonth in Var_std.keys():
        print('month: ', imonth)
        for v in var_list:
            print(Tcount[imonth][v])
            Var_mean[imonth][v] = Var_mean[imonth][v]/Tcount[imonth][v]
            Var_std[imonth][v] = np.sqrt(Var_std[imonth][v]/(Tcount[imonth][v]-1) - 
                                        Tcount[imonth][v]*(Var_mean[imonth][v]**2)/(Tcount[imonth][v]-1))
        sp.io.savemat('./bg_std/roms_ecs_v4sc_std_i_%d.mat'%(imonth),Var_std[imonth])



#    print(Var_filt)#
    
    #    p.close()
#    import matplotlib.pyplot as plt
#    sp.io.savemat('testmp_filt_%s.mat'%(v),{'%s_std'%(v): Var_filt})    
#    np.save('testmp_filt_%s'%(v),Var_filt)
