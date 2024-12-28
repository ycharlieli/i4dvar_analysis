import matplotlib.pyplot as plt
import Obs_workspace as ow

datasets = {
           'NSFC2012SPRING-YS':
               {
                'PSAL2':['H03','H16','H25','H34'],
                'TEMP':[]
                'stations':[
                            ['E01','E02','E03','E04','E05'],
                            ['H01','H02','H03','H04','H05','H06','H07','H08'],
                            ['H25','H26','H27','H28','H29'],
                            ['H33','H34','H35','H36'],
                            ['H38','H39','H40','H41','H42'],
                            ['HF1','HF2','HF3','H44','H43',]
                           ],
                },
           'NSFC2012SPRING-ES':{'PSAL2':[],'TEMP':[]},
           'IOCAS2012SPRING':{'PSAL2':['3600-06','3800-01','3800-02','B-03','DH6-2'],'TEMP':['DH6-2',]},
           'NSFC2012AUTUMN-ES':{'PSAL2':[],'TEMP':[]},
           'IOCAS2012AUTUMN':{'PSAL2':['3500-09',],'TEMP':['3500-09','DH4-1']},
           'NSFC2013SPRING-CJ':{'PSAL2':['A01-03','A01-09','A02-04','A02-08','A03-06','A03-09','A06-05','A08-06','A13-01'],'TEMP':[]},
           'IOCAS2013SPRING':{'PSAL2':['3600-05','3800-04','B-04'],'TEMP':[]}
           }


for dataset in list(datasets.keys()):

