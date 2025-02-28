import xarray as xr
from dask.distributed import Client
import time
import datetime as dt
import warnings
warnings.filterwarnings('ignore')
import sys
import numpy as np
sys.setrecursionlimit(100000)

class roms_workspace():
    def __init__(self, my_root, wrk_name, prefix='ocean_ecs', file_type='his'):
        self.my_root     = my_root
        self.my_wrk      = wrk_name
        self.my_datapath = my_root+ 

