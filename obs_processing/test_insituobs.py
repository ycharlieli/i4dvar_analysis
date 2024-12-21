import seabird
import datetime
import numpy as np
import scipy as sp
import seapy
import gsw
from seabird.cnv import fCNV
import  matplotlib.pyplot as plt
import os

CTDDIR='IOCAS2012SPRING'
os.chdir('/Volumes/TO_1/roms4dvar_ecs/i4dvar_outputs/INSITU_OBS/%s/CTD/'%(CTDDIR))

ctd_files = seapy.list_files('./*acfld.cnv')

def lanczos_filter(data,Cf,Nf,M):
    coef = sp.signal.firwin(M+1, Cf/Nf,width = 2/len(data),window='lanczos', pass_zero='lowpass')
    return sp.signal.filtfilt(coef,1.0,data)

dt = 0.04 #sec
fs = 1/(dt)
Nf = fs/2
Cf = fs/1500
M = 100
VAR = 'PSAL2'
for ctd_file in ctd_files:
    print(ctd_file)
    profile = fCNV(ctd_file)
    invalid =np.ma.getmaskarray(profile['flag'])
    time = profile['timeS'][:]
    #temp = profile[VAR][:]
    temp = profile.data[-6].data
    pres = profile['PRES'][:]
    # smooth the noise further
    pres_conv = np.convolve(pres,np.ones(24*6)/24/6,mode='full')[72:-71]
    pres_conv = seapy.filt.bandpass(pres_conv,dt,low_cutoff=dt*200,order=8)
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
                  pres[start_up[downcast2_index]]) > 5):
            time_downcast = np.r_[
                                  time[start_up[downcast_index]:start_down[downcast_index+1]],
                                  time[start_up[downcast2_index]:start_down[downcast2_index+1]]

                                  ]
            pres_downcast = np.r_[
                                  pres[start_up[downcast_index]:start_down[downcast_index+1]],
                                  pres[start_up[downcast2_index]:start_down[downcast2_index+1]]
                                  ]
            temp_downcast = np.r_[
                                   temp[start_up[downcast_index]:start_down[downcast_index+1]],
                                   temp[start_up[downcast2_index]:start_down[downcast2_index+1]]
                                ]
            invalid_downcast = np.r_[
                                   invalid[start_up[downcast_index]:start_down[downcast_index+1]],
                                   invalid[start_up[downcast2_index]:start_down[downcast2_index+1]]
                                ]
        elif (downcast2_index < downcast_index)\
            and (((pres[start_down[downcast2_index+1]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres) *0.4) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index+1]]) <np.max(time)*0.1))\
                or (((pres[start_down[downcast2_index+1]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres) *0.2) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index+1]]) <np.max(time)*0.05)):
            time_downcast = np.r_[
                                  time[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                  time[start_up[downcast_index]:start_down[downcast_index+1]]

                                  ]
            pres_downcast = np.r_[
                                  pres[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                  pres[start_up[downcast_index]:start_down[downcast_index+1]]
                                  ]
            temp_downcast = np.r_[
                                   temp[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                   temp[start_up[downcast_index]:start_down[downcast_index+1]]
                                ]
            invalid_downcast = np.r_[
                                   invalid[start_up[downcast2_index]:start_down[downcast2_index+1]],
                                   invalid[start_up[downcast_index]:start_down[downcast_index+1]]
                                ]
        else:
                
            time_downcast = time[start_up[downcast_index]:start_down[downcast_index+1]]
            pres_downcast = pres[start_up[downcast_index]:start_down[downcast_index+1]]
            temp_downcast = temp[start_up[downcast_index]:start_down[downcast_index+1]]
            invalid_downcast = invalid[start_up[downcast_index]:start_down[downcast_index+1]]



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
                  pres[start_up[downcast2_index]]) > 5):
            time_downcast = np.r_[
                                  time[start_up[downcast_index]:start_down[downcast_index]],
                                  time[start_up[downcast2_index]:start_down[downcast2_index]]

                                  ]
            pres_downcast = np.r_[
                                  pres[start_up[downcast_index]:start_down[downcast_index]],
                                  pres[start_up[downcast2_index]:start_down[downcast2_index]]
                                  ]
            temp_downcast = np.r_[
                                   temp[start_up[downcast_index]:start_down[downcast_index]],
                                   temp[start_up[downcast2_index]:start_down[downcast2_index]]
                                ]
            invalid_downcast = np.r_[
                                   invalid[start_up[downcast_index]:start_down[downcast_index]],
                                   invalid[start_up[downcast2_index]:start_down[downcast2_index]]
                                ]
        elif (downcast2_index < downcast_index)\
            and (((pres[start_down[downcast2_index]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres)*0.4) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index]]) <np.max(time)*0.1))\
              or (((pres[start_down[downcast2_index]] - 
                  pres[start_up[downcast2_index]]) > np.max(pres)*0.2) \
                and ((time[start_up[downcast_index]] -
                  time[start_down[downcast2_index]]) <np.max(time)*0.05)) :
            time_downcast = np.r_[
                                  time[start_up[downcast2_index]:start_down[downcast2_index]],
                                  time[start_up[downcast_index]:start_down[downcast_index]]

                                  ]
            pres_downcast = np.r_[
                                  pres[start_up[downcast2_index]:start_down[downcast2_index]],
                                  pres[start_up[downcast_index]:start_down[downcast_index]]
                                  ]
            temp_downcast = np.r_[
                                   temp[start_up[downcast2_index]:start_down[downcast2_index]],
                                   temp[start_up[downcast_index]:start_down[downcast_index]]
                                ]
            invalid_downcast = np.r_[
                                   invalid[start_up[downcast2_index]:start_down[downcast2_index]],
                                   invalid[start_up[downcast_index]:start_down[downcast_index]]
                                ]
        else:
                
            time_downcast = time[start_up[downcast_index]:start_down[downcast_index]]
            pres_downcast = pres[start_up[downcast_index]:start_down[downcast_index]]
            temp_downcast = temp[start_up[downcast_index]:start_down[downcast_index]]
            invalid_downcast = invalid[start_up[downcast_index]:start_down[downcast_index]]


#    print(start_up[downcast_index])
    # print(np.argmax(np.abs(start_down-start_up)))
    
    plt.plot(time, pres,'b')
    plt.plot(time,pres_conv,'k--')
    plt.plot(time_downcast,pres_downcast,'r')
    plt.plot(time[start_up],pres[start_up],'r.')
    plt.plot(time[start_down],pres[start_down],'g.')
    plt.xlabel('TIME')
    plt.ylabel('PRES')
    plt.title('%s'%(ctd_file[2:-9]))
    #plt.show()
    plt.savefig('../fig_downcast_detection/%s_%sdowncast_detection.jpeg'%(ctd_file[2:-4],VAR), dpi=400)
    plt.close()
    pres_intp = np.arange(np.ceil(np.min(pres_downcast)), np.floor(np.max(pres_downcast)))
    temp_intp = np.interp(pres_intp,pres_downcast,temp_downcast)
    plt.plot(temp_downcast,pres_downcast,'k', label = 'Raw')
    plt.plot(temp_intp, pres_intp,'r',label='Bin Averaged')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.xlabel(VAR)
    plt.ylabel('PRES')
    plt.title('%s'%(ctd_file[2:-9]))
    #plt.show()
    plt.savefig('../fig_downcast_avg/%s_%sdowncast_avg.jpeg'%(ctd_file[2:-4],VAR), dpi=400)
    plt.close()



