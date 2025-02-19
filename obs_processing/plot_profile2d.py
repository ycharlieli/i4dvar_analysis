import matplotlib.pyplot as plt
import Obs_workspace as ow

datasets = {
           'NSFC2012SPRING-YS':
           {
                'PSAL2':['H03','H16','H25','H34'],
                'TEMP':[],
                'Profiles':[
                            ['E01','E02','E03','E04','E05'],
                            ['H01','H02','H03','H04','H05','H06','H07','H08'],
                            ['H25','H26','H27','H28','H29'],
                            ['H33','H34','H35','H36'],
                            ['H38','H39','H40','H41','H42'],
                            ['HF1','HF2','HF3','H44','H43',]
                           ],
            },
           'NSFC2012SPRING-ES':
           {
               'PSAL2':[],
               'TEMP':[],
               'Profiles':[
                           ['DH3-01','DH3-02','DH3-03','DH3-04','DH3-05','DH3-06','DH3-07','DH3-08',],
                           ['RB-10','RB-16'],
                           ['DH4-01','DH4-02','DH4-03','DH4-04','DH4-05','DH4-06',],
                           ['ZA-05','ZA-03','ZA-01'],
                           ['GC-05','ZB-07','ZB-09','ZB-11','DH5-01','GC-06','ZB-12','DH5-02','DH5-03','DH5-04','DH5-05','DH5-06'],
                           ['DH6-05','DH6-04','DH6-03','DH6-02','ZC-17','DH6-01','ZC-15','ZC-13'],
                           ['ZD-19','ZD-21','DH7-01','ZD-24','DH7-02','DH7-03','DH7-04','DH7-05'],
                           ['DH8-05','DH8-04','DH8-03','DH8-02','ZE-27','DH8-01','ZE-25'],
                           ['GC-07', 'DH9-01','DH9-02','DH9-03','DH9-04','DH9-05','DH9-06']
                           ]
            },
           'IOCAS2012SPRING':
           {
                'PSAL2':['3600-06','3800-01','3800-02','B-03','DH6-2'],
                'TEMP':['DH6-2',],
                'Profiles': [
                             ['DH9-1','DH9-2','DH9-3','DH9-4','DH9-5','DH9-6'],
                             ['DH8-5','DH8-4','DH8-3','DH8-2','DH8-1'],
                             ['DH7-1','DH7-2','DH7-3','DH7-4','DH7-5'],
                             ['DH6-5','DH6-4','DH6-3','DH6-2','DH6-1'],
                             ['DH5-1','DH5-2','DH5-3','DH5-4','DH5-5','DH5-6'],
                             ['DH4-6','DH4-5','DH4-4','DH4-3','DH4-2','DH4-1'],
                             ['DH3-1','DH3-2','DH3-3','DH3-4','DH3-5','DH3-6','DH3-7','DH3-8'],
                             ['DH2-7','DH2-6','DH2-5','DH2-4','DH2-3','DH2-2','DH2-1'],
                             ['DH1-7','DH1-8'],
                             ['CJ-1','CJ-2','CJ-3','CJ-4','CJ-5','CJ-6'],
                             ['3300-06','3300-05','3300-04','3300-03','3300-02','3300-01'],
                             ['3400-01','3400-02','3400-03','3400-04','3400-05','3400-06','3400-07','3400-08'],
                             ['3500-10','3500-09','3500-08','3500-07','3500-06','3500-05','3500-04','3500-03','3500-02','3500-01'],
                             ['3600-01','3600-02','3600-03','3600-04','3600-05','3600-06','3600-07','3600-08'],
                             ['B-08','B-07','B-06','B-05','B-04','B-03','B-02','B-01'],
                             ['3875-01','3875-02','3875-03','3875-04','3875-05'],
                             ]
            },
           'NSFC2012AUTUMN-ES':{'PSAL2':[],'TEMP':[]},
           'IOCAS2012AUTUMN':{'PSAL2':['3500-09',],'TEMP':['3500-09','DH4-1']},
           'NSFC2013SPRING-CJ':{'PSAL2':['A01-03','A01-09','A02-04','A02-08','A03-06','A03-09','A06-05','A08-06','A13-01'],'TEMP':[]},
           'IOCAS2013SPRING':{'PSAL2':['3600-05','3800-04','B-04'],'TEMP':[]}
           }


#for dataset in list(datasets.keys()):

