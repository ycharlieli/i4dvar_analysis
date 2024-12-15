import multiprocessing as mp
import os


def fun0101():
    os.system('python -u mp_filt_roms.py 0101 > mp_out.0101')

def fun0102():
    os.system('python -u mp_filt_roms.py 0102 > mp_out.0102')

def fun0103():
    os.system('python -u mp_filt_roms.py 0103 > mp_out.0103')

def fun0104():
    os.system('python -u mp_filt_roms.py 0104 > mp_out.0104')

def fun0105():
    os.system('python -u mp_filt_roms.py 0105 > mp_out.0105')

def fun0201():
    os.system('python -u mp_filt_roms.py 0201 > mp_out.0201')

def fun0202():
    os.system('python -u mp_filt_roms.py 0202 > mp_out.0202')

def fun0203():
    os.system('python -u mp_filt_roms.py 0203 > mp_out.0203')

def fun0204():
    os.system('python -u mp_filt_roms.py 0204 > mp_out.0204')

def fun0205():
    os.system('python -u mp_filt_roms.py 0205 > mp_out.0205')

def fun0301():
    os.system('python -u mp_filt_roms.py 0301 > mp_out.0301')

def fun0302():
    os.system('python -u mp_filt_roms.py 0302 > mp_out.0302')

def fun0303():
    os.system('python -u mp_filt_roms.py 0303 > mp_out.0303')

def fun0304():
    os.system('python -u mp_filt_roms.py 0304 > mp_out.0304')

def fun0305():
    os.system('python -u mp_filt_roms.py 0305 > mp_out.0305')

def fun0401():
    os.system('python -u mp_filt_roms.py 0401 > mp_out.0401')

def fun0402():
    os.system('python -u mp_filt_roms.py 0402 > mp_out.0402')

def fun0403():
    os.system('python -u mp_filt_roms.py 0403 > mp_out.0403')

def fun0404():
    os.system('python -u mp_filt_roms.py 0404 > mp_out.0404')

def fun0405():
    os.system('python -u mp_filt_roms.py 0405 > mp_out.0405')

def fun0501():
    os.system('python -u mp_filt_roms.py 0501 > mp_out.0501')

def fun0502():
    os.system('python -u mp_filt_roms.py 0502 > mp_out.0502')

def fun0503():
    os.system('python -u mp_filt_roms.py 0503 > mp_out.0503')

def fun0504():
    os.system('python -u mp_filt_roms.py 0504 > mp_out.0504')

def fun0505():
    os.system('python -u mp_filt_roms.py 0505 > mp_out.0505')

if __name__ == '__main__':
    process0101 = mp.Process(target=fun0101)
    process0102 = mp.Process(target=fun0102)
    process0103 = mp.Process(target=fun0103)
    process0104 = mp.Process(target=fun0104)
    process0105 = mp.Process(target=fun0105)
    process0201 = mp.Process(target=fun0201)
    process0202 = mp.Process(target=fun0202)
    process0203 = mp.Process(target=fun0203)
    process0204 = mp.Process(target=fun0204)
    process0205 = mp.Process(target=fun0205)
    process0301 = mp.Process(target=fun0301)
    process0302 = mp.Process(target=fun0302)
    process0303 = mp.Process(target=fun0303)
    process0304 = mp.Process(target=fun0304)
    process0305 = mp.Process(target=fun0305)
    process0401 = mp.Process(target=fun0401)
    process0402 = mp.Process(target=fun0402)
    process0403 = mp.Process(target=fun0403)
    process0404 = mp.Process(target=fun0404)
    process0405 = mp.Process(target=fun0405)
    process0501 = mp.Process(target=fun0501)
    process0502 = mp.Process(target=fun0502)
    process0503 = mp.Process(target=fun0503)
    process0504 = mp.Process(target=fun0504)
    process0505 = mp.Process(target=fun0505)

    process0101.start()
    process0102.start()
    process0103.start()
    process0104.start()
    process0105.start()
    process0201.start()
    process0202.start()
    process0203.start()
    process0204.start()
    process0205.start()
    process0301.start()
    process0302.start()
    process0303.start()
    process0304.start()
    process0305.start()
    process0401.start()
    process0402.start()
    process0403.start()
    process0404.start()
    process0405.start()
    process0501.start()
    process0502.start()
    process0503.start()
    process0504.start()
    process0505.start()


    process0101.join()
    process0102.join()
    process0103.join()
    process0104.join()
    process0105.join()
    process0201.join()
    process0202.join()
    process0203.join()
    process0204.join()
    process0205.join()
    process0301.join()
    process0302.join()
    process0303.join()
    process0304.join()
    process0305.join()
    process0401.join()
    process0402.join()
    process0403.join()
    process0404.join()
    process0405.join()
    process0501.join()
    process0502.join()
    process0503.join()
    process0504.join()
    process0505.join()

