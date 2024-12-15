import numpy as np
import scipy as sp
def lowpass_cosine_filter_coef(Cf,M):
    coef = Cf*np.r_[1 ,np.sin(np.pi*np.arange(1,M+1)*Cf)/(np.pi*np.arange(1,M+1)*Cf)]
    
    return coef

def lanczos_filter_coef(Cf,M):
    hkcs = lowpass_cosine_filter_coef(Cf,M)
    sigma = np.r_[1 ,np.sin(np.pi*np.arange(1,M+1)/M)/(np.pi*np.arange(1,M+1)/M)]
    hkB = hkcs*sigma
    hkA = -hkB; hkA[0] = hkA[0]+1
    coef = np.c_[hkB,hkA]
    
    return coef
    
def spectral_window(coef,N):
    Ff = np.arange(0,1+2/N,2/N);
    window = np.zeros([len(Ff),])
    for i in range(len(Ff)):
        window[i] = coef[0] + 2*np.sum(coef[1:]*np.cos(np.arange(1,len(coef))*np.pi*Ff[i]))
        
    return window, Ff


def spectral_filtering(x,window):
    Nx = len(x);
    Cx = sp.fft.fft(x)
    # print(Cx.shape)
    Cx = Cx[:int(np.floor(Nx/2)+1)]
    # print(Cx.shape)
    CxH = Cx*window
    # print(CxH.shape)
    CxH = np.r_[CxH,np.conj(CxH[Nx-len(CxH):0:-1])]
    # print(CxH.shape)
    y = np.real(sp.fft.ifft(CxH))
    
    return y,Cx



def lanczos_filter(data,Cf,Nf,M):
    coef = sp.signal.firwin(M+1, Cf/Nf,width = 2/len(data),window='lanczos', pass_zero='lowpass')
    return sp.signal.filtfilt(coef,1.0,data)

