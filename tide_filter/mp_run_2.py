from concurrent import futures
import os

def mp_filt(pid):
    os.system('python -u mp_filt_roms.py %s > mp_out.%s'%(pid, pid))


with futures.ProcessPoolExecutor(max_workers=25) as executor:
    for i in range(1,6):
        for j in range(1,6):
            executor.submit(mp_filt, str(i).rjust(2,'0')+str(j).rjust(2,'0'))



