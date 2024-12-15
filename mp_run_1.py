from multiprocessing import Pool
import os


def mp_filt(pid):
    os.system('python -u mp_filt_roms.py %s > mp_out.%s'%(pid, pid))



if __name__  == '__main__':
    p = Pool(8)

    for i in range(1,3):
        for j in range(1,5):
            p.apply_async(mp_filt,args=(str(i).rjust(2,'0')+str(j).rjust(2,'0'),))

    p.close()
    p.join()
