import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import imageio
import napari
import math
from PIL import Image

def linear(x,a,b):
    return a*x+b

def gaussian(x,a,s,m,b):
    return a*np.exp(-(x-m)**2/(2*s**2))+b

def lorentzian(x, a, s, m, b):
    return a/(2*np.pi)*s/((x-m)**2+(1/2*s)**2)+b

def format_string(directory,name, i):
    #function to create file names for multiple frames automatically
    rem=len(str(i))
    output=directory+name+'frame'+(5-rem)*'0'+str(i)+'.tif'
    return output

def reshape(spec):
    #function to reassemble 4tap convoluted spectrum into one single spectrum
    b=len(spec)//4
    spec1=np.zeros(len(spec))
    tap1=spec[0:b]
    tap2=spec[b:2*b]
    tap3=spec[2*b:3*b]
    tap4=spec[3*b:4*b]
    spec1[0::4]=tap1
    spec1[1::4]=tap2
    spec1[2::4]=tap3
    spec1[3::4]=tap4
    return spec1

def reshape2D(spec, n=4, axis=1):
    #only call if len(spec)%n==0
    if axis==0:
        a=len(spec[0])
        b=len(spec)//n
    if axis==1:
        a=len(spec)
        b=len(spec[0])//n
    spec1=np.zeros((a, b*n))
    for i in range(n):
        if axis==0:
            tap=spec[i*b:(i+1)*b]
            spec1[i::n]=tap
        if axis==1:
            tap=spec[:,i*b:(i+1)*b]
            spec1[:,i::n]=tap
    return spec1

def formAxis(gratingType=0, gratingCentre=600, width=2048):
    #Calibration function of wavelength for different gratings
    axis=np.arange(0, 2048)
    #find correct calibration
    if gratingType==0: #300 l/mm
        a=-6.42990401e-02
        b0=6.65503505e+02
        l0=600
    elif gratingType==1: #1200 l/mm
        a=-1.34444474e-02
        b0=7.13788374e+02
        l0=700
    else: #no valid grating type: return pixel values
        a=1
        b0=0
        l0=0
    xaxis=a*axis+b0+(gratingCentre-l0)
    if width==2048:
        return xaxis
    # smaller camera width settings
    crop=(2048-width)//2 #take care of different camera sizes
    return xaxis[crop:-crop]

def formatString(directory,name, i):
    rem=len(str(i))
    output=directory+name+'frame'+(5-rem)*'0'+str(i)+'.tif'
    return output

def fitSpectrumWL(spec, wL,tol, tol2, height):
    peaks, props=find_peaks(spec, distance=tol2, height=height)
    heights=props['peak_heights']
    peakFit=[]
    heightsCorr=[]
    for i in range(len(peaks)):
        peak=peaks[i]
        try:
            popt, pcov = curve_fit(gaussian, wL[math.floor(peak)-tol:math.ceil(peak)+tol], spec[math.floor(peak)-tol:math.ceil(peak)+tol], p0=[250, 0.1, wL[peak], 4])
            peakFit.append(popt[2])  
            heightsCorr.append(popt[0])
        except:
            popt=[-1, -1, -1, -1] #fit failed
    return peakFit, heightsCorr

def construct2D(path, xdim, ydim, xaxis, s1, s2, threshold=20, specImg=True):
    data=imageio.imread(path)
    Imap=np.zeros((xdim, ydim))
    Maxmap=np.zeros((xdim, ydim))
    for x in range (0,xdim):
        for y in range (0,ydim):
            spec=data[y*xdim+x]
            I=np.max(spec)
            if I>threshold and specImg:
                specCrop=reshape(spec)[s1:s2]
                lamCrop=xaxis[s1:s2]
                lam=lamCrop[np.min(np.where(specCrop==np.max(specCrop)))]
                Maxmap[x][y]=lam
            Imap[x][y]=I
    return Imap, Maxmap

def show2DMaps(imap0, maxmap0, min0, max0, scaleVal=200):
    fig, ax = plt.subplots(1,2, figsize=(14, 7))
    ax1, ax2 = ax
    map1=ax1.imshow(imap0, vmax=scaleVal, cmap='magma')
    map2=ax2.imshow(maxmap0, cmap='nipy_spectral', vmin=min0, vmax=max0)
    fig.colorbar(map1, ax=ax1)
    fig.colorbar(map2, ax=ax2)
    return fig

def VAtoRGBA(v0, a, start, end):
    v = v0 - start
    Hue0 = v / (end-start)
    Hue= Hue0 * (Hue0>=0) * (Hue0<=1) #clip values below 0 and above 1
    #find bool arrays of which math to be applied
    caseI=Hue<1/3
    caseII=(Hue>=1/3)*(Hue<2/3)
    caseIII=(Hue>=2/3)
    #find Hn matrix
    Hn=Hue-1/3*caseII-2/3*caseIII
    #find different color components
    compI=a*(1+np.cos(2*np.pi*Hn)/np.cos(2*np.pi*(1/6-Hn)))
    compII=0
    compIII=3*a-compI-compII
    #sum these correctly
    R=compI*caseI+compII*caseII+compIII*caseIII
    G=compII*caseI+compIII*caseII+compI*caseIII
    B=compIII*caseI+compI*caseII+compII*caseIII
    return R, G, B

def normArr(r, g, b):
    myA=np.array([r, g, b])
    myA+=np.min(myA)
    myA/=np.max(myA)
    r=myA[0]
    g=myA[1]
    b=myA[2]
    return r, g, b

def rescale(a, lowLim, upLim):
    #inverts and rescales 16bit images into 8bit unsigned integer for tracking algorithm use
    a-=lowLim
    a*=255
    a/=(upLim-lowLim)
    a=(a<0)*0+(a>=0)*a
    a=(a<=255)*a+(a>255)*255
    return np.uint8(a)

def binStack(stack, nBin):
    stackBinned=np.zeros(np.shape(stack)//nBin)
    for i in range(nBin[0]):
        for j in range(nBin[1]):
            for k in range(nBin[2]):
                stackBinned+=(stack[i::nBin[0], j::nBin[1], k::nBin[2]])/(np.average(nBin))
    return stackBinned