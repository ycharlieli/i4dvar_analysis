#%%
import sys
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import scipy.io as sio

#%%
# assign the target folder name
suffix='nonud'
OUTPUT_dir = '/Volumes/Elements SE/roms4dvar_ecs/i4dvar_outputs/workspace/STORAGE'
Inpb_dir = OUTPUT_dir + '/prior'
Inpa_dir = OUTPUT_dir + '/posterior'
#Inpn_files =  ['/Volumes/Elements SE/roms4dvar_ecs/freerun_detide/ocean_ecsdtd_avg_'\
#             +str(i).rjust(4,'0')+'.nc' for i in range(289,293)]
Inpn_files = ['/Volumes/Elements SE/roms4dvar_ecs/run17/outputs/ocean_ecs_avg_'\
              +str(i).rjust(4,'0')+'.nc' for i in range(289,291)] 
# assign the target area
area = sys.argv[1]
FREE_RUN = True
Inpb_prefix = 'Inpb_'+area+'_'+suffix
Inpa_prefix = 'Inpa_'+area+'_'+suffix
Inpn_prefix = 'Inpn_'+area

time_range =slice('2017-01-03','2017-02-02')
if area == 'cre':
    eta_range = slice(90,180)
    xi_range  = slice(15,100)
elif area == 'ys':
    eta_range = slice(220,305)
    xi_range  = slice(105,200)
elif area == 'bh':
    eta_range = slice(280,360)
    xi_range  = slice(0,80)
elif area == 'ks':
    eta_range = slice(10,120)
    xi_range  = slice(120,220)
else:
    print('''This area has not been included. ONLY
            Changjiang River Estuary: 'cre'
            Yellow Sea              : 'ys'
            BoHai Sea               : 'bh'
            Kuroshio Current        : 'ks'
            is VALID.''')
    sys.exit()

#%%
# build the pointer of multi-ncfiles
Inpb_ds = xr.open_mfdataset(Inpb_dir+'/ocean_ecs_fwd*.nc',combine='nested',\
                            coords='minimal',compat='override',concat_dim='ocean_time',\
                            data_vars=['temp','salt','u','v','zeta'],parallel=True)

Inpa_ds = xr.open_mfdataset(Inpa_dir+'/ocean_ecs_fwd*.nc',combine='nested',\
                            coords='minimal',compat='override',concat_dim='ocean_time',\
                            data_vars=['temp','salt','u','v','zeta'],parallel=True)

Inpn_ds = xr.open_mfdataset(Inpn_files,combine='nested',\
                            coords='minimal',compat='override',concat_dim='ocean_time',\
                            data_vars=['temp','salt','u','v','zeta'],parallel=True)

#%%
# extract from background
Inpb = Inpb_ds.sel(eta_rho=eta_range,xi_rho=xi_range,\
                   eta_u=eta_range,xi_u=xi_range,\
                   eta_v=eta_range,xi_v=xi_range)
Inpb = Inpb.sel(ocean_time=time_range)
rlon = Inpb.lon_rho.data.compute()
rlat = Inpb.lat_rho.data.compute()
ulon = Inpb.lon_u.data.compute()
ulat = Inpb.lat_u.data.compute()
vlon = Inpb.lon_v.data.compute()
vlat = Inpb.lat_v.data.compute()


# save mat file
#sio.savemat(Inpb_mat,{'rlon':Inpb.lon_rho.data.compute(),\
#                                  'rlat':Inpb.lat_rho.data.compute(),\
#                                  'zeta':Inpb.zeta.data.compute(),\
#                                  'temp':Inpb.temp.data.compute(), 'salt':Inpb.salt.data.compute(),\
#                                  'u':Inpb.u.data.compute(),'v':Inpb.v.data.compute()})

#%%
#extract from analysis
Inpa = Inpa_ds.sel(eta_rho=eta_range,xi_rho=xi_range,\
                   eta_u=eta_range,xi_u=xi_range,\
                   eta_v=eta_range,xi_v=xi_range)
Inpa = Inpa.sel(ocean_time=time_range)
timeab = Inpa.ocean_time.data.astype('float64')/86400*1e-9

#save mat file
#sio.savemat(Inpa_mat,{'rlon':Inpa.lon_rho.data.compute(),\
#                                  'rlat':Inpa.lat_rho.data.compute(),\
#                                  'zeta':Inpa.zeta.data.compute(),\
#                                  'temp':Inpa.temp.data.compute(), 'salt':Inpa.salt.data.compute(),\
#                                  'u':Inpa.u.data.compute(),'v':Inpa.v.data.compute()})
#%%
#extract from free fun
#if FREE_RUN:
Inpn = Inpn_ds.sel(eta_rho=eta_range,xi_rho=xi_range,\
                   eta_u=eta_range,xi_u=xi_range,\
                   eta_v=eta_range,xi_v=xi_range)
Inpn = Inpn.sel(ocean_time=time_range)
timen = Inpn.ocean_time.data.astype('float64')/86400*1e-9
#    
#    sio.savemat(Inpn_mat,{'rlon':Inpn.lon_rho.data.compute(),\
#                                  'rlat':Inpn.lat_rho.data.compute(),\
#                                  'zeta':Inpn.zeta.data.compute(),\
#                                 'temp':Inpn.temp.data.compute(), 'salt':Inpn.salt.data.compute(),\
#                                  'u':Inpn.u.data.compute(),'v':Inpn.v.data.compute()})
for i,varname in enumerate(['zeta','temp','salt','u','v']):
    var_mat = 'Inp_'+area+'_'+suffix+'_'+varname+'_rdrgd03'+'.mat'
    print('Saving '+var_mat+' ...')

    varb = eval("Inpb."+varname+".data.compute()")
    vara = eval("Inpa."+varname+".data.compute()")
    varn = eval("Inpn."+varname+".data.compute()")
    #print("Inpb."+varname+".data.compute()")
    #print("Inpa."+varname+".data.compute()")
    #print("Inpn."+varname+".data.compute()")

    eval("sio.savemat(var_mat,{'rlon':rlon,"\
                             +"'rlat':rlat,"\
                             +"'ulon':ulon,"\
                             +"'ulat':ulat,"\
                             +"'vlon':vlon,"\
                             +"'vlat':vlat,"\
                             +"'timeab':timeab,"\
                             +"'timen':timen,"\
                             +"'"+varname+"b'"+":varb,"\
                             +"'"+varname+"a'"+":vara,"\
                             +"'"+varname+"n'"+":varn})")


