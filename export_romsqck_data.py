import roms_workspace as rw
import scipy as sci
import os
##========================================================================##
#  Stick to the guideline of beginning with the simplest approach,
#  only sea surface temperature (SST) was assimilated, and only
#  initial condition was adjusted.
#  Here we first validated the performance of ROMS I4D-VAR configured
#  in this strategy.
##========================================================================##
my_rootnl = '/Volumes/Elements SE/roms4dvar_ecs'
my_rootda = '/Volumes/Elements SE/roms4dvar_ecs/i4dvar_outputs'
nl_workspace_info = {
                     'run04/outputs':'FORWARD',
                     }
da_workspace_info = {
                     'workspace02':'Adj INI',
                     }

var = 'temp'
areas = [
         'bh',
#         'ys',
#         'cre',
#         'ks'
        ]

nl_Vars = {}
nl_Workspace = rw.Nl_workspace(          MY_ROOT=my_rootnl,
                                  workspace_name=list(nl_workspace_info.keys())[0],
                                workspace_detail=list(nl_workspace_info.values())[0]
                               )
nl_Workspace.load_mtfroms(isAna=False,Surf=True)

daa_Vars = {}
daa_Workspace = rw.Da_workspace(         MY_ROOT=my_rootda, 
                                  workspace_name= list(da_workspace_info.keys())[0],
                                workspace_detail= list(da_workspace_info.values())[0]
                               )
daa_Workspace.load_mtfroms(isAna=True, Surf=True)

obs_ds = xr.open_dataset(os.path.join(my_rootda,
                                      list(da_workspace_info.keys())[0],
                                      'OBS',
                                      'ecs_oisst_super_201701_full.nc'
                                      ))
obs_time = obs_ds.obs_time.data.astype('float')/86400/1e9
obs_value = obs_ds.obs_value.data
obs_lon = obs_ds.obs_lon.data
obs_lat = obs_ds.obs_lat.data
obs_Vars = {} 
# area iteration 
for area in areas:
    print('loading %s in %s  area...'%(var,area))
    # loading data from FORWARD model outputs
    nl_Vars[area]  =  nl_Workspace.pull_var2d(         varname=var, area=area,time_range=('2017-01-03','2017-02-02'))
    # loading data from posterior outputs
    daa_Vars[area] = daa_Workspace.pull_var2d(varname='%s_sur'%var, area=area,time_range=('2017-01-03','2017-02-02'))
    index_insea = np.where(~np.isnan(daa_Vars[area].value[0,:,:]))
    # extract data from observation file
    obs_Vars[area] = np.zeros(
                            [len(obs_ds.survey_time.data), 
                             len(index_insea)]
                            ) 
    for itime, itime_obs in enumerate(obs_ds.survey_time.data.astype('float')/86400/1e9):
        index_thistime = np.where(obs_time == itime_obs)[0]
        obs_Vars[area][itime,:] = sci.interpolate.griddata(
                                                           np.c_[obs_lon[index_thistime],
                                                                obs_lat[index_thistime]],
                                                           obs_value[index_thistime],
                                                           np.c_[daa_Vars[area].lon[index_insea],
                                                                 daa_Vars[area].lat[index_insea]],
                                                           method ='nearest')


