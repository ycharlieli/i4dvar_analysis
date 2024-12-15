import numpy as np
import xarray as xr
import scipy.io as sio
import os 

class Nl_workspace():
    def __init__(self, MY_ROOT, workspace_name, workspace_detail=None):
        self.my_root = MY_ROOT
        self.workspace_name = workspace_name
        self.workspace_detail = workspace_detail
        self.workspace_subdir = os.path.join(self.my_root, workspace_name)
        print('Workspace: ' + self.workspace_subdir)
        # given the simulation time?
    #build sub class for variable. 
    # var.name , var.lon, var.lat, var.value.....
    class Var():
        def __init__(self, name=None, dim=None):
            self.name = name
            self.dim =  dim + 1 #time dim

    def start_pull2d(self,var, time_range = None):
        coords2d = {'rho':('eta_rho', 'xi_rho'),
                      'u':(  'eta_u',   'xi_u'),
                      'v':(  'eta_v',   'xi_v')}
        var.time_range = slice(time_range[0], time_range[1])
        if ('temp' in var.name) or ( 'salt' in var.name) or ('zeta' in var.name ):
            using_coord = 'rho'
        elif var.name == 'u' or var.name == 'u_sur':
            using_coord = 'u'
        elif var.name == 'v' or var.name == 'v_sur':
            using_coord = 'v'

        if self.isAna:
            ds = 'ana'
        else:
            ds = 'bg'
        if (var.name == 'u' or var.name == 'u_sur') or (var.name == 'v' or var.name == 'v_sur'):
            expr_value = "self.%s_ds.%s.sel(ocean_time=var.time_range,%s=var.eta_range,%s=var.xi_range).drop_duplicates(dim='ocean_time',keep='last').data.compute()"\
                  %(ds,var.name,
                    coords2d[using_coord][0],
                    coords2d[using_coord][1])
            expr_lon = 'self.%s_ds.lon_%s.sel(%s=var.eta_range,%s=var.xi_range).data.compute()'\
                    %(ds,
                      using_coord,
                      coords2d[using_coord][0],
                      coords2d[using_coord][1]
                      )
            expr_lat = 'self.%s_ds.lat_%s.sel(%s=var.eta_range,%s=var.xi_range).data.compute()'\
                    %(ds,
                      using_coord,
                      coords2d[using_coord][0],
                      coords2d[using_coord][1]
                      )
            expr_time = "self.%s_ds.ocean_time.sel(ocean_time=var.time_range).drop_duplicates(dim='ocean_time',keep='last')"%(ds)
            #print(' ' ,expr_lon,'\n',expr_lat,'\n',expr_time,'\n',expr_value)
            var.lon = eval(expr_lon); var.lat= eval(expr_lat); var.ocean_time = eval(expr_time); var.value=eval(expr_value)
            expr_lonrho = 'self.%s_ds.lon_%s.sel(%s=var.eta_range,%s=var.xi_range).data.compute()'\
                    %(ds,
                      'rho',
                      'eta_rho',
                      'xi_rho' 
                      )
            expr_latrho = 'self.%s_ds.lat_%s.sel(%s=var.eta_range,%s=var.xi_range).data.compute()'\
                    %(ds,
                      'rho',
                      'eta_rho',
                      'xi_rho' 
                      )
            expr_angle = 'self.%s_ds.angle.sel(%s=var.eta_range,%s=var.xi_range).data.compute()'\
                    %(ds,
                      'eta_rho',
                      'xi_rho' 
                      )
            var.lonrho = eval(expr_lonrho); var.latrho = eval(expr_latrho); var.angle = eval(expr_angle)
            var_rho = np.zeros((var.value.shape[0],var.value.shape[1], var.value.shape[2]))
            for t in range(len(var.ocean_time)):
                for i in range(0,var.lon.shape[1]-1):
                    for j in range(0,var.lon.shape[0]-1):
                        if var.name == 'u' or var.name == 'u_sur':
                            var_rho[t,j,i] = 0.5*(var.value[t,j,i+1] + var.value[t,j,i])
                        elif var.name == 'v' or var.name == 'v_sur':
                            var_rho[t,j,i] = 0.5*(var.value[t,j+1,i] + var.value[t,j,i])
            var.lon = var.lonrho[:-1,:-1]; var.lat = var.latrho[:-1,:-1]; var.angle = var.angle[:-1,:-1]; var.value = var_rho[:,:-1,:-1]



        else:
            expr_value = "self.%s_ds.%s.sel(ocean_time=var.time_range,%s=var.eta_range,%s=var.xi_range).drop_duplicates(dim='ocean_time',keep='last').data.compute()"\
                  %(ds,var.name,
                    coords2d[using_coord][0],
                    coords2d[using_coord][1])
            expr_lon = 'self.%s_ds.lon_%s.sel(%s=var.eta_range,%s=var.xi_range).data.compute()'\
                    %(ds,
                      using_coord,
                      coords2d[using_coord][0],
                      coords2d[using_coord][1]
                      )
            expr_lat = 'self.%s_ds.lat_%s.sel(%s=var.eta_range,%s=var.xi_range).data.compute()'\
                    %(ds,
                      using_coord,
                      coords2d[using_coord][0],
                      coords2d[using_coord][1]
                      )
            expr_time = "self.%s_ds.ocean_time.sel(ocean_time=var.time_range).drop_duplicates(dim='ocean_time',keep='last')"%(ds)
            #print(' ' ,expr_lon,'\n',expr_lat,'\n',expr_time,'\n',expr_value)
            var.lon = eval(expr_lon); var.lat= eval(expr_lat); var.ocean_time = eval(expr_time); var.value=eval(expr_value)
        return var
    # load single file?
    # def load_sfroms
    def load_mtfroms(self,isAna=False,Surf=True):
        self.isAna = isAna
        if not self.isAna:
            
            if Surf:
                try:
                    self.bg_ds = xr.open_mfdataset(
                                                     os.path.join(self.workspace_subdir,'ocean_ecs_avg*.nc'),
                                                     combine='nested',
                                                     coords='minimal',compat='override',concat_dim='ocean_time',
                                                     data_vars=['temp','salt','u','v','zeta'],
                                                     parallel=True
                                                     )
                except:
                    self.bg_ds = xr.open_mfdataset(
                                                     os.path.join(self.workspace_subdir,'ocean_ecs_his*.nc'),
                                                     combine='nested',
                                                     coords='minimal',compat='override',concat_dim='ocean_time',
                                                     data_vars=['temp','salt','u','v','zeta'],
                                                     parallel=True
                                                     )
                self.bg_ds = self.bg_ds.isel(s_rho=-1)
            else:
                self.bg_ds = xr.open_mfdataset(
                                               os.path.join(self.workspace_subdir,'ocean_ecs_avg*.nc'),
                                               combine='nested',
                                               coords='minimal',compat='override',concat_dim='ocean_time',
                                               data_vars=['temp','salt','u','v','zeta'],
                                               parallel=True
                                               )

    # pull var profile
    # def pull_varprof

    # pull var transection
    # def pull_vartransc
    
    def pull_var2d(self,varname=None, time_range = None,s_rho = -1, eta_range=None, xi_range=None, area=None):
        var = self.Var(name=varname,dim=2)

        if not eta_range:
            if not area:
                print('Warning: if eta_range, xi_range is not specified, area need to be specified. vice versa.')
            else:
                var.area = area
                areaPull = True
                rangePull = False
        else:
            rangePull = True 
            areaPull = False

        if not varname:
             print('Warning: variable name need to be specified.')
        else:
             varPull = True

        if varPull:
            if varname not in ['temp', 'salt', 'u', 'v', 'zeta']:
                print('Warning: ' + varname + ' may not in the dataset.')
            

            if areaPull:
                if var.area == 'cre':
#                    var.eta_range = slice(90,180)
#                    var.xi_range = slice(15,100)
                     var.eta_range = slice(70,220)
                     var.xi_range = slice(15,150)
                elif var.area == 'ys':
                    var.eta_range = slice(220,305)
                    var.xi_range = slice(105,200)
                elif var.area == 'bh':
                    var.eta_range = slice(280,360)
                    var.xi_range = slice(0,80)
                elif var.area == 'ks':
                    var.eta_range = slice(10,120)
                    var.xi_range = slice(151,230)
                elif var.area == 'all':
                    var.eta_range = slice(0,361)
                    var.xi_range = slice(0,241)
                else:
                    print('no such area, using default cre area.')
                    var.eta_range = slice(90,180)
                    var.xi_range = slice(15,100)

            if rangePull:
                var.eta_range = eta_range
                var.xi_range = xi_range

            var = self.start_pull2d(var,time_range)

        return var

class Da_workspace(Nl_workspace):
    def __init__(self, MY_ROOT, workspace_name, workspace_detail=None):
        super(Da_workspace, self).__init__(MY_ROOT, workspace_name, workspace_detail)

        self.prior_path = os.path.join(self.workspace_subdir, 'STORAGE/prior')
        #print('\t    Prior: ' + self.prior_path)
        
        self.posterior_path = os.path.join(self.workspace_subdir, 'STORAGE/posterior')
        #print('\tPosterior: ' + self.posterior_path)
        


               
    def load_mtfroms(self,isAna=True,Surf=False):
        self.isAna = isAna
        if not self.isAna:
            
            if Surf:
                self.bg_ds = xr.open_mfdataset(
                                                 os.path.join(self.prior_path,'ocean_ecs_qck*.nc'),
                                                 combine='nested',
                                                 coords='minimal',compat='override',concat_dim='ocean_time',
                                                 data_vars=['temp_sur','salt_sur','u_sur','v_sur','zeta'],
                                                 parallel=True
                                                 )
            else:
                self.bg_ds = xr.open_mfdataset(
                                               os.path.join(self.prior_path,'ocean_ecs_fwd*.nc'),
                                               combine='nested',
                                               coords='minimal',compat='override',concat_dim='ocean_time',
                                               data_vars=['temp','salt','u','v','zeta'],
                                               parallel=True
                                               )
        if self.isAna:    
            
            if Surf:
                self.ana_ds = xr.open_mfdataset(
                                                 os.path.join(self.posterior_path,'ocean_ecs_qck*.nc'),
                                                 combine='nested',
                                                 coords='minimal',compat='override',concat_dim='ocean_time',
                                                 data_vars=['temp_sur','salt_sur','u_sur','v_sur','zeta'],
                                                 parallel=True
                                                 )
            else:
                self.ana_ds = xr.open_mfdataset(
                                           os.path.join(self.posterior_path,'ocean_ecs_fwd*.nc'),
                                           combine='nested',
                                           coords='minimal',compat='override',concat_dim='ocean_time',
                                           data_vars=['temp','salt','u','v','zeta'],
                                           parallel=True
                                           )



    def pull_var2d(self,varname=None,time_range = None, eta_range=None, xi_range=None, area=None):
        var = self.Var(name=varname,dim=2)

        if not eta_range:
            if not area:
                print('Warning: if eta_range, xi_range is not specified, area need to be specified. vice versa.')
            else:
                var.area = area
                areaPull = True
                rangePull = False
        else:
            rangePull = True 
            areaPull = False

        if not varname:
             print('Warning: variable name need to be specified.')
        else:
             varPull = True

        if varPull:
            if varname not in ['temp_sur', 'salt_sur', 'u_sur', 'v_sur', 'zeta']:
                print('Warning: ' + varname + ' may not in the dataset.')
            

            if areaPull:
                if var.area == 'cre':
#                    var.eta_range = slice(90,180)
#                    var.xi_range = slice(15,100)
                     var.eta_range = slice(70,220)
                     var.xi_range = slice(15,150)
                elif var.area == 'ys':
                    var.eta_range = slice(220,305)
                    var.xi_range = slice(105,200)
                elif var.area == 'bh':
                    var.eta_range = slice(280,360)
                    var.xi_range = slice(0,80)
                elif var.area == 'ks':
                    var.eta_range = slice(10,120)
                    var.xi_range = slice(151,230)
                elif var.area == 'all':
                    var.eta_range = slice(0,361)
                    var.xi_range = slice(0,241)
                else:
                    print('no such area, using default cre area.')
                    var.eta_range = slice(90,180)
                    var.xi_range = slice(15,100)

            if rangePull:
                var.eta_range = eta_range
                var.xi_range = xi_range

            var = self.start_pull2d(var,time_range)

        return var
            

    #def pull_profile(...)
# class Obs_workspace(Nl_workspace):
# Since only SST has been assimilated. data extraction is quite simple.
# this part needs put more thought on the design of framework
