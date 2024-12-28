import datetime
import os
import re
import numpy as np
from seabird import fCNV
import seapy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import scipy.io as sio
import cmaps
from matplotlib.colors import LinearSegmentedColormap
from cotede import datasets, qctests

class insitu_workspace():
    def __init__(self,MY_ROOT,worksapce_name, workspace_detail=None):
        self.my_root = MY_ROOT
        self.workspace_name = worksapce_name
        self.workspace_detail = workspace_detail
        self.workspace_subdir = os.path.join(self.my_root, self.workspace_name)
        print('Workspace: ' + self.workspace_subdir)
        self.ctd_files = seapy.list_files(os.path.join(self.workspace_subdir+'/CTD/*acfld.cnv'))
        self.ctd_stationinfo = os.path.join(self.workspace_subdir,self.workspace_name+'-station.txt')
        self.isempty = True
        
        
    class Data():
        def __init__(self):
            self.filetype = 'ctd data structure'
            




    def pull_downcast(self,ctd_file,VAR):
        def lanczos_filter(data,Cf,Nf,M):
            coef = sp.signal.firwin(M+1, Cf/Nf,width = 2/len(data),window='lanczos', pass_zero='lowpass')
            return sp.signal.filtfilt(coef,1.0,data)

        dt = 0.04 #sec
        fs = 1/(dt)
        Nf = fs/2
        Cf = fs/1500
        M = 100
        # print(ctd_file)
        profile = fCNV(ctd_file)
        time    = profile['timeS'][:]
        if VAR == 'PSAL2':
            var = profile.data[-6][:]
        else:
            var     = profile[VAR][:]
        pres    = profile['PRES'][:]
        depth   = profile['DEPTH'][:]
        # smooth the noise further
        # print(np.nanmax(time))
        # smooth data based on the data length 
        if np.nanmax(time)>3000:
            pres_conv = np.convolve(pres,np.ones(24*6*5)/24/6/5,mode='full')[72*5:-72*5+1]
            pres_conv = seapy.filt.bandpass(pres_conv,dt,low_cutoff=dt*2000,order=4)
        elif (self.workspace_name == 'IOCAS2013SPRING') and (('3400-05' in ctd_file) or ('3300-03' in ctd_file)):
            # specified data file, not follow the general rules.
            pres_conv = np.convolve(pres,np.ones(24*3)/24/3,mode='full')[36:-35]
            pres_conv = seapy.filt.bandpass(pres_conv,dt,low_cutoff=dt*50,order=12)
        elif (np.nanmax(time) > 370)  or ((self.workspace_name == 'NSFC2012SPRING-ES') and ('ZB-13' in ctd_file)):
             # specified data file, not follow the general rules.
            pres_conv = np.convolve(pres,np.ones(24*6)/24/6,mode='full')[72:-71]
            pres_conv = seapy.filt.bandpass(pres_conv,dt,low_cutoff=dt*200,order=8)
        else:
            pres_conv = np.convolve(pres,np.ones(24*3)/24/3,mode='full')[36:-35]
            pres_conv = seapy.filt.bandpass(pres_conv,dt,low_cutoff=dt*50,order=12)
        
        diff_pres = np.diff(pres_conv)
    #    diff_pres_conv = np.convolve(diff_pres,np.ones(24*6)/24/6,mode='full')[72:-71]
        diff_pres_conv = diff_pres.copy()
        count_down = 0
        count_up = 0
        for i in range(1,len(diff_pres_conv)):
            # print(diff_pres[i])
            # print(diff_pres[i]*diff_pres[i-1])
            if (diff_pres_conv[i]*diff_pres_conv[i-1]) < 0 :
                if diff_pres_conv[i-1]>0:
                    if  not count_down:
                        start_down = i
                    else:
                        start_down = np.r_[start_down, i]
                    count_down+=1
                else:
                    if not count_up :
                        start_up = i
                    else:
                        start_up = np.r_[start_up, i]
                    count_up+=1


    #        print(start_up)
    #        print(start_down)


        # print(pres[start_down[1]])
        larger = 1
        smaller = -1
        equal = 0
    #    print(len(start_down))
    #    print(len(start_up))
        if len(start_up) > len(start_down):
            end_indx = len(start_down)
            status = larger
        elif len(start_up) < len(start_down):
            end_indx = len(start_up)
            status = smaller
        else:
            status = equal


        if start_down[0] < start_up[0]:
            if status == equal:
                mono_duration = start_up[:-1] - start_down[1:]
            elif status == larger:
                mono_duration = start_up[:end_indx-1] - start_down[1:]
            elif status == smaller:
                mono_duration = start_up[:] - start_down[1:]
    #        print(mono_duration)
            downcast_index = np.argmax(np.abs(mono_duration))
    #        print(downcast_index)
            if ((pres[start_down[downcast_index+1]] -
                      pres[start_up[downcast_index]]) < 5):
                mono_duration[downcast_index] = 0
                downcast_index = np.argmax(np.abs(mono_duration))
            mono_duration[downcast_index] = 0

            downcast2_index = np.argmax(np.abs(mono_duration))

            if (downcast2_index > downcast_index)\
                and ((pres[start_down[downcast2_index+1]] - 
                      pres[start_up[downcast2_index]]) > 5)\
                and (pres[start_up[downcast2_index]]>pres[start_down[downcast_index+1]]):
                time_downcast = np.r_[
                                      time[start_up[downcast_index]:start_down[downcast_index+1]],
                                      time[start_up[downcast2_index]:start_down[downcast2_index+1]]

                                      ]
                pres_downcast = np.r_[
                                      pres[start_up[downcast_index]:start_down[downcast_index+1]],
                                      pres[start_up[downcast2_index]:start_down[downcast2_index+1]]
                                      ]
                depth_downcast = np.r_[
                                      depth[start_up[downcast_index]:start_down[downcast_index+1]],
                                      depth[start_up[downcast2_index]:start_down[downcast2_index+1]]
                                      ]
                var_downcast = np.r_[
                                       var[start_up[downcast_index]:start_down[downcast_index+1]],
                                       var[start_up[downcast2_index]:start_down[downcast2_index+1]]
                                    ]
            elif (downcast2_index < downcast_index)\
                   and((pres[start_down[downcast2_index+1]] - 
                  pres[start_up[downcast2_index]]) < np.max(pres) *0.7)\
                and ((((pres[start_down[downcast2_index+1]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres) *0.4) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index+1]]) <np.max(time)*0.1))\
                or (((pres[start_down[downcast2_index+1]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres) *0.2) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index+1]]) <np.max(time)*0.05))\
                or((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index+1]]) <np.max(time)*0.005)): # judge that the interval between two downcast is small enough
                time_downcast = np.r_[
                                      time[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                      time[start_up[downcast_index]:start_down[downcast_index+1]]

                                      ]
                pres_downcast = np.r_[
                                      pres[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                      pres[start_up[downcast_index]:start_down[downcast_index+1]]
                                      ]
                depth_downcast = np.r_[
                                      depth[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                      depth[start_up[downcast_index]:start_down[downcast_index+1]]
                                      ]
                var_downcast = np.r_[
                                       var[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                       var[start_up[downcast_index]:start_down[downcast_index+1]]
                                    ]
            else:

                time_downcast = time[start_up[downcast_index]:start_down[downcast_index+1]]
                pres_downcast = pres[start_up[downcast_index]:start_down[downcast_index+1]]
                depth_downcast = depth[start_up[downcast_index]:start_down[downcast_index+1]]
                var_downcast = var[start_up[downcast_index]:start_down[downcast_index+1]]
                
            



        else:
            if status == equal:
                mono_duration = start_up - start_down
            elif status == larger:
                mono_duration = start_up[:end_indx] - start_down
            elif status == smaller:
                mono_duration = start_up[:] - start_down[:end_indx]

    #        print(mono_duration)
            downcast_index = np.argmax(np.abs(mono_duration))
    #        print(downcast_index)
            if ((pres[start_down[downcast_index]] -
                      pres[start_up[downcast_index]]) < 5):
                mono_duration[downcast_index] = 0
                downcast_index = np.argmax(np.abs(mono_duration))
            mono_duration[downcast_index] = 0
            downcast2_index = np.argmax(np.abs(mono_duration))
            if (downcast2_index > downcast_index)\
                and ((pres[start_down[downcast2_index]] - 
                      pres[start_up[downcast2_index]]) > 5)\
                and (pres[start_up[downcast2_index]]>pres[start_down[downcast_index]]):
                time_downcast = np.r_[
                                      time[start_up[downcast_index]:start_down[downcast_index]],
                                      time[start_up[downcast2_index]:start_down[downcast2_index]]

                                      ]
                pres_downcast = np.r_[
                                      pres[start_up[downcast_index]:start_down[downcast_index]],
                                      pres[start_up[downcast2_index]:start_down[downcast2_index]]
                                      ]
                depth_downcast = np.r_[
                                      depth[start_up[downcast_index]:start_down[downcast_index]],
                                      depth[start_up[downcast2_index]:start_down[downcast2_index]]
                                      ]
                
                var_downcast = np.r_[
                                       var[start_up[downcast_index]:start_down[downcast_index]],
                                       var[start_up[downcast2_index]:start_down[downcast2_index]]
                                    ]
            elif (downcast2_index < downcast_index)\
                   and ((pres[start_down[downcast2_index]] - 
                  pres[start_up[downcast2_index]]) < np.max(pres)*0.7)\
                 and ((((pres[start_down[downcast2_index]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres)*0.4) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index]]) <np.max(time)*0.1))\
              or (((pres[start_down[downcast2_index]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres)*0.2) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index]]) <np.max(time)*0.05))\
                or((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index]]) <np.max(time)*0.005)): # judge that the interval between two downcast is small enough
                time_downcast = np.r_[
                                      time[start_up[downcast2_index]:start_down[downcast2_index]],
                                      time[start_up[downcast_index]:start_down[downcast_index]]

                                      ]
                pres_downcast = np.r_[
                                      pres[start_up[downcast2_index]:start_down[downcast2_index]],
                                      pres[start_up[downcast_index]:start_down[downcast_index]]
                                      ]
                depth_downcast = np.r_[
                                      depth[start_up[downcast2_index]:start_down[downcast2_index]],
                                      depth[start_up[downcast_index]:start_down[downcast_index]]
                                      ]
                var_downcast = np.r_[
                                       var[start_up[downcast2_index]:start_down[downcast2_index]],
                                       var[start_up[downcast_index]:start_down[downcast_index]]
                                    ]
            else:

                time_downcast = time[start_up[downcast_index]:start_down[downcast_index]]
                pres_downcast = pres[start_up[downcast_index]:start_down[downcast_index]]
                depth_downcast = depth[start_up[downcast_index]:start_down[downcast_index]]
                var_downcast = var[start_up[downcast_index]:start_down[downcast_index]]
            
            
        return time_downcast.data, pres_downcast.data, depth_downcast.data, var_downcast.data

    def pull_ctdinfos(self):
        self.ctd_infos = {}
        ctdid = open(self.ctd_stationinfo)
        info = ctdid.readline().split()
        cid = 0
        while info:
            self.ctd_infos[info[1]]={}
            self.ctd_infos[info[1]]['TIME'] = (datetime.datetime.strptime(info[2],'%Y%m%d%H%M') - 
                                              datetime.datetime(1970,1,1)).total_seconds()/3600/24
            self.ctd_infos[info[1]]['LONGITUDE']=info[3]
            self.ctd_infos[info[1]]['LATITUDE']=info[4]
            self.ctd_infos[info[1]]['DEPTH'] = info[5]
            self.ctd_infos[info[1]]['CTDID'] = cid
            cid +=1
            info = ctdid.readline().split()
        ctdid.close()
        
    def pull_alldata(self,VARS, refine_withqc_id=None):
        if self.isempty:
            self.data = {}
            self.isempty = False
        for iv, var in enumerate(VARS):
            print(var)
            self.data[var] = self.Data()
            if len(refine_withqc_id[var]) > 0:
                print('qc specified: ')
                print(refine_withqc_id[var])
                needqc=True
            else:
                needqc=False
            for ic, ctd_file in enumerate(self.ctd_files):
                # print(ctd_file)
                
                ctdid = re.findall(r"%s/CTD/(.*)acfld.cnv"%(self.workspace_subdir),ctd_file)[0]
                # print(ctdid)
                # print(needqc)
                if needqc:
                    if ctdid in refine_withqc_id[var]:
                        # print('qc is activated')
                        idneedqc=True
                    else:
                        idneedqc=False
                else:
                    idneedqc=False
                
                if ic == 0:
                    _, self.data[var].pres, self.data[var].depth, self.data[var].value = self.pull_downcast(ctd_file,var)
                    #remove negative depth
                    valid_depth = np.where(self.data[var].depth>0)
                    self.data[var].pres = self.data[var].pres[valid_depth]
                    self.data[var].depth = self.data[var].depth[valid_depth]
                    self.data[var].value = self.data[var].value[valid_depth]
                    self.data[var].station = np.tile(ctdid,[len(self.data[var].value),]).ravel()
                    self.data[var].time = np.tile(float(self.ctd_infos[ctdid]['TIME']),[len(self.data[var].value),]).ravel()
                    self.data[var].lon = np.tile(float(self.ctd_infos[ctdid]['LONGITUDE']),[len(self.data[var].value),]).ravel()
                    self.data[var].lat = np.tile(float(self.ctd_infos[ctdid]['LATITUDE']),[len(self.data[var].value),]).ravel()
                    
                    # interpolating
                    
                    # print('interp is activate')

                    data_len = int(np.floor(0.95*len(self.data[var].value)))
                    pres_cut = self.data[var].pres[:data_len] 
                    # splitted based on pressure
                    # cut the tail of data
                    pres_levels = np.arange(np.ceil(np.min(pres_cut)),
                                               np.floor(np.max(pres_cut)),
                                               )

                    self.data[var].station_i = np.tile(ctdid,[len(pres_levels),]).ravel()
                    self.data[var].time_i = np.tile(float(self.ctd_infos[ctdid]['TIME']),[len(pres_levels),]).ravel()
                    self.data[var].lon_i = np.tile(float(self.ctd_infos[ctdid]['LONGITUDE']),[len(pres_levels),]).ravel()
                    self.data[var].lat_i = np.tile(float(self.ctd_infos[ctdid]['LATITUDE']),[len(pres_levels),]).ravel()
                    self.data[var].depth_i = seapy.seawater.depth(pres_levels, lat=float(self.ctd_infos[ctdid]['LATITUDE']))
                    self.data[var].pres_i = pres_levels.copy()
                    self.data[var].value_i = np.interp(pres_levels,pres_cut,self.data[var].value[:data_len])
                    
                    # refining
                    grad_value_r = np.diff(self.data[var].value_i)
                    grad_value_r = np.abs(np.r_[grad_value_r[0],grad_value_r])
                    grad_value_l = np.diff(self.data[var].value_i)
                    grad_value_l = np.abs(np.r_[grad_value_l,grad_value_l[-1]])
                    grad_value = (grad_value_r+grad_value_l)/2
                    if var =='PSAL2':
                        grad_value = 1/(1+np.exp(-grad_value**2))
                    elif var =='TEMP':
                        grad_value = 1/(1+np.exp(-grad_value))
                    
                    refine_levels = np.linspace(np.min(pres_levels),np.max(pres_levels),1000)
                        
    
                    if idneedqc:
                        
                        print("quality control is activate: %s"%(ctdid))
                        value_qc = self.data[var].value_i[grad_value<np.mean(grad_value)]
                        pres_qc = self.data[var].pres_i[grad_value<np.mean(grad_value)]
                            
                    else:
                        value_qc = self.data[var].value_i.copy()
                        pres_qc = self.data[var].pres_i.copy()
                    

                    self.data[var].value_r = np.interp(refine_levels,pres_qc,value_qc)
                    self.data[var].pres_r = refine_levels.copy()
                    self.data[var].depth_r = seapy.seawater.depth(refine_levels, lat=float(self.ctd_infos[ctdid]['LATITUDE']))
                    self.data[var].time_r =  np.tile(float(self.ctd_infos[ctdid]['TIME']),[len(refine_levels),]).ravel()
                    self.data[var].lon_r =  np.tile(float(self.ctd_infos[ctdid]['LONGITUDE']),[len(refine_levels),]).ravel()
                    self.data[var].lat_r = np.tile(float(self.ctd_infos[ctdid]['LATITUDE']),[len(refine_levels),]).ravel()
                    self.data[var].station_r = np.tile(ctdid,[len(refine_levels),]).ravel()




                        
                else:
                    
                    _, thispres, thisdepth, thisvalue = self.pull_downcast(ctd_file,var)
                    #remove negative depth
                    valid_depth = np.where(thisdepth>0)
                    thispres = thispres[valid_depth]
                    thisvalue = thisvalue[valid_depth]
                    thisdepth = thisdepth[valid_depth]
                    thisstation = np.tile(ctdid,[len(thisvalue),]).ravel()
                    thistime = np.tile(float(self.ctd_infos[ctdid]['TIME']),[len(thisvalue),]).ravel()
                    thislon = np.tile(float(self.ctd_infos[ctdid]['LONGITUDE']),[len(thisvalue),]).ravel()
                    thislat = np.tile(float(self.ctd_infos[ctdid]['LATITUDE']),[len(thisvalue),]).ravel()
                    self.data[var].station  = np.r_[self.data[var].station ,thisstation]
                    self.data[var].time  = np.r_[self.data[var].time ,thistime]
                    self.data[var].pres  = np.r_[self.data[var].pres ,thispres]
                    self.data[var].depth = np.r_[self.data[var].depth,thisdepth]
                    self.data[var].value = np.r_[self.data[var].value,thisvalue]
                    self.data[var].lon   = np.r_[self.data[var].lon,thislon]
                    self.data[var].lat   = np.r_[self.data[var].lat,thislat]
                    
                    # interpolating
                    # print('interp is activate')

                    data_len = int(np.floor(0.95*len(thisvalue)))
                    pres_cut = thispres[:data_len] 
                    # splitted based on pressure
                    # cut the tail of data
                    pres_levels = np.arange(np.ceil(np.min(pres_cut)),
                                               np.floor(np.max(pres_cut)),
                                               )

                    thisstation_i = np.tile(ctdid,[len(pres_levels),]).ravel()
                    thistime_i = np.tile(float(self.ctd_infos[ctdid]['TIME']),[len(pres_levels),]).ravel()
                    thislon_i = np.tile(float(self.ctd_infos[ctdid]['LONGITUDE']),[len(pres_levels),]).ravel()
                    thislat_i = np.tile(float(self.ctd_infos[ctdid]['LATITUDE']),[len(pres_levels),]).ravel()
                    thisdepth_i = seapy.seawater.depth(pres_levels, lat=float(self.ctd_infos[ctdid]['LATITUDE']))
                    thispres_i = pres_levels.copy()
                    thisvalue_i = np.interp(pres_levels,pres_cut,thisvalue[:data_len])
                    self.data[var].station_i = np.r_[self.data[var].station_i, thisstation_i]
                    self.data[var].time_i = np.r_[self.data[var].time_i, thistime_i]
                    self.data[var].lon_i = np.r_[self.data[var].lon_i, thislon_i]
                    self.data[var].lat_i = np.r_[self.data[var].lat_i, thislat_i]
                    self.data[var].depth_i = np.r_[self.data[var].depth_i, thisdepth_i]
                    self.data[var].pres_i = np.r_[self.data[var].pres_i, thispres_i]
                    self.data[var].value_i = np.r_[self.data[var].value_i,thisvalue_i]
                        
                    
                    # refining  
                    
                    grad_value_r = np.diff(thisvalue_i)
                    grad_value_r = np.abs(np.r_[grad_value_r[0],grad_value_r])
                    grad_value_l = np.diff(thisvalue_i)
                    grad_value_l = np.abs(np.r_[grad_value_l,grad_value_l[-1]])
                    grad_value = (grad_value_r+grad_value_l)/2
                    if var =='PSAL2':
                        grad_value = 1/(1+np.exp(-grad_value**2))
                    elif var =='TEMP':
                        grad_value = 1/(1+np.exp(-grad_value))
                    
                    
                    refine_levels = np.linspace(np.min(pres_levels),np.max(pres_levels),1000)


                            
                          
                    if idneedqc:
                        print("quality control is activate: %s"%(ctdid))
                        value_qc = thisvalue_i[grad_value<np.mean(grad_value)]
                        pres_qc = thispres_i[grad_value<np.mean(grad_value)]
                    else:
                        value_qc = thisvalue_i.copy()
                        pres_qc = thispres_i.copy()

                    thisvalue_r = np.interp(refine_levels,pres_qc,value_qc)
                    thispres_r = refine_levels.copy()
                    thisdepth_r = seapy.seawater.depth(refine_levels, lat=float(self.ctd_infos[ctdid]['LATITUDE']))
                    thistime_r =  np.tile(float(self.ctd_infos[ctdid]['TIME']),[len(refine_levels),]).ravel()
                    thislon_r =  np.tile(float(self.ctd_infos[ctdid]['LONGITUDE']),[len(refine_levels),]).ravel()
                    thislat_r = np.tile(float(self.ctd_infos[ctdid]['LATITUDE']),[len(refine_levels),]).ravel()
                    thisstation_r = np.tile(ctdid,[len(refine_levels),]).ravel()
                    self.data[var].station_r = np.r_[self.data[var].station_r, thisstation_r]
                    self.data[var].time_r = np.r_[self.data[var].time_r, thistime_r]
                    self.data[var].lon_r = np.r_[self.data[var].lon_r, thislon_r]
                    self.data[var].lat_r = np.r_[self.data[var].lat_r, thislat_r]
                    self.data[var].depth_r = np.r_[self.data[var].depth_r, thisdepth_r]
                    self.data[var].pres_r = np.r_[self.data[var].pres_r, thispres_r]
                    self.data[var].value_r = np.r_[self.data[var].value_r,thisvalue_r]
                    
    
        
            
                    
                    
                   
                
                
            
        