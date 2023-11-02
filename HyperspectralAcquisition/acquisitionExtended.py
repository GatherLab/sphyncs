import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import imageio
import napari
import math
import utils
from IPython.display import display, HTML
from PIL import Image


class acquisition:
    '''Acquisition Data Class to store relevant parameters and make them accessible for processing'''
    
    def __init__(self, dataDirectory, dataName, x0=200, y0=200, z0=20, t0=50, dX=1, dY=1, dZ=1, 
                 dT=60, grating=0, wLC=600, wL0=2048):
        self.name=dataName
        self.directory=dataDirectory
        self.xDim=x0
        self.yDim=y0
        self.zDim=z0
        self.tDim=t0
        self.xStep=dX
        self.yStep=dY
        self.zStep=dZ
        self.tStep=dT
        self.xaxis=utils.formAxis(grating, wLC)
        self.wLDim=wL0
        
    def autoCalibrate(self, logFileName, FOV50=1733.85):
        '''Reads Logfile and extracts acquisition parameters to store in Acquisition Object'''
        logFile=open(self.directory+logFileName)
        acqParams=logFile.readlines()
        self.xDim=int(acqParams[-11].split('\t')[0])
        self.yDim=int(acqParams[-10].split('\t')[0])
        self.zDim=int(acqParams[-7].split('\t')[0])
        self.tDim=(len(acqParams)-16)//self.zDim
        if self.tDim==0: self.tDim=1
        nullDat=imageio.imread(utils.formatString(self.directory, self.name, 0))
        self.wLDim=len(nullDat[0])
        self.xStep=FOV50/100*(int(acqParams[-14].split('\t')[0])-int(acqParams[-16].split('\t')[0]))/self.xDim
        self.yStep=FOV50/100*(int(acqParams[-13].split('\t')[0])-int(acqParams[-15].split('\t')[0]))/self.yDim
        self.zStep=-int(acqParams[-8].split('\t')[0])/1000
        try: self.tStep=1/float(acqParams[-9].split('\t')[0])*self.zDim/60
        except: self.tStep=1 #calculating time per z-stack in min from buffer frame rate
        grating=int(acqParams[-3].split('\t')[0])
        wLC=int(acqParams[-2].split('\t')[0])
        self.xaxis=utils.formAxis(grating, wLC, self.wLDim)
        
                
    def autoCalibrateCCF(self, logFileName, FOV50=1733.85):
        '''Reads Logfile and extracts acquisition parameters to store in Acquisition Object'''
        logFile=open(logFileName)
        acqParams=logFile.readlines()
        self.xDim=int(acqParams[-11].split('\t')[0])
        self.zDim=1
        self.yDim=int(acqParams[-7].split('\t')[0])
        self.tDim=(len(acqParams)-16)//(self.yDim*self.zDim)
        nullDat=imageio.imread(utils.formatString(self.directory, self.name, 0))
        self.wLDim=len(nullDat[0])
        self.xStep=FOV50/100*(int(acqParams[-14].split('\t')[0])-int(acqParams[-16].split('\t')[0]))/self.xDim
        self.yStep=FOV50/100*(int(acqParams[-13].split('\t')[0])-int(acqParams[-15].split('\t')[0]))/self.yDim
        self.zStep=-int(acqParams[-8].split('\t')[0])/1000
        try: self.tStep=1/float(acqParams[-9].split('\t')[0])*self.yDim*self.zDim/60
        except: self.tStep=1 #calculating time per z-stack in min from buffer frame rate
        grating=int(acqParams[-3].split('\t')[0])
        wLC=int(acqParams[-2].split('\t')[0])
        self.xaxis=utils.formAxis(grating, wLC, self.wLDim)
    #============================ Image Processing =========================================    
    def batchProcessND(self, nBin=[5,5,5,10], upLim=500, lowLim=15):
        for t in range(self.tDim):
            print('processing stack {} out of {}'.format(t+1, self.tDim))
            for z in range(0, self.zDim, nBin[3]):
                image0=(imageio.imread(utils.format_string(self.directory, self.name, t*self.zDim+z)))/nBin[3]
                for z2 in range(z+1, z+nBin[3]):
                    image0+=(imageio.imread(utils.format_string(self.directory, self.name, t*self.zDim+z2))/nBin[3])
                im2=utils.reshape2D(image0, n=4, axis=1)
                stack=np.reshape(im2, (self.xDim, self.yDim, self.wLDim))
                nBin0=nBin[:3]
                stackBinned=utils.binStack(stack, nBin0)
                stackNew=utils.rescale(stackBinned, lowLim, upLim)
                for w in range(np.shape(stackNew)[2]):
                    Image.fromarray(255-stackNew[:,:,w]).save(self.directory+'processedND\\'+self.name+'Imap_t{}_w{}_z{}.tif'.format(t,w,z//nBin[3]))
                    
    def batchProcessNDCCF(self, nBin=[5,5,5,10], upLim=500, lowLim=15, show=False):
        for t in range(self.tDim):
            print('processing stack {} out of {}'.format(t+1, self.tDim))
            for z in range(0, self.zDim, nBin[3]):
                image0=np.zeros((self.xDim*self.yDim, self.wLDim))
                for y in range(0, self.yDim):
                    image00=(imageio.imread(utils.format_string(self.directory, self.name, t*self.zDim*self.yDim+self.yDim*z+y)))/nBin[3]
                    image0[y*self.xDim:(y+1)*self.xDim]=image00
                for z2 in range(z+1, z+nBin[3]):
                    image0+=(imageio.imread(utils.format_string(self.directory, self.name, t*self.zDim*self.yDim+self.yDim*z2+y))/nBin[3])
                #im2=utils.reshape2D(image0, n=4, axis=1)
                stack=np.reshape(image0, (self.xDim, self.yDim, self.wLDim))
                nBin0=nBin[:3]
                stackBinned=utils.binStack(stack, nBin0)
                print(np.shape(stackBinned))
                if show: 
                    plt.figure()
                    plt.imshow(stackBinned[:,:,0])
                    plt.savefig(self.directory+'\\test.png')
                    plt.show()
                stackNew=utils.rescale(stackBinned, lowLim, upLim)
                for w in range(np.shape(stackNew)[2]):
                    Image.fromarray(255-stackNew[:,:,w]).save(self.directory+'processedND\\'+self.name+'Imap_t{}_w{}_z{}.tif'.format(t,w,z//nBin[3]))
    
    def construct2D(self, z=0, t=0,  s1=400, s2=600, threshold=20, specImg=False, show=False, save=False, scaleVal=200):
        '''constructor to create intensity and spectral maps'''
        path=utils.formatString(self.directory, self.name, t*self.zDim+z)
        Imap, Maxmap=utils.construct2D(path, self.xDim, self.yDim, self.xaxis, s1, s2, threshold, specImg)
        if show:
            fig=utils.show2DMaps(Imap, Maxmap, self.xaxis[s2],  self.xaxis[s1], scaleVal=scaleVal)
            if save:
                fig.savefig(self.directory+'view2D_t{}_z{}'.format(t,z)+self.name)
        if specImg:
            return Imap, Maxmap
        return Imap
    
    def constructSpec2DCCF(self, z=0, t=0,  s1=400, s2=600, threshold=20, specImg=False, show=False, save=False, scaleVal=200):
        '''constructor to create intensity and spectral maps'''
        Imap = np.zeros((self.xDim, self.yDim))
        Maxmap = np.zeros((self.xDim, self.yDim))
        for y in range(self.yDim):
            path=utils.formatString(self.directory, self.name, t*self.zDim*self.yDim+self.yDim*z+y)
            ImapY, MaxmapY=utils.construct2D(path, self.xDim, 1, self.xaxis, s1, s2, threshold, specImg)
            Imap[:,y]=ImapY[:,0]
            Maxmap[:,y]=MaxmapY[:,0]
        if show:
            fig=utils.show2DMaps(Imap, Maxmap, self.xaxis[s2],  self.xaxis[s1], scaleVal=scaleVal)
            if save:
                fig.savefig(self.directory+'view2D_t{}_z{}'.format(t,z)+self.name)
        if specImg:
            return Imap, Maxmap
        return Imap
    
    def construct3D(self, t=0, s1=400, s2=600, threshold=20, save=False):
        imaps=[]
        for z in range(self.zDim):
            imap=self.construct2D( z, t,  s1=s1, s2=s2, threshold=20, specImg=False, show=False, save=False)
            imaps.append(imap)
        return np.array(imaps)
    
    def constructT(self, z=None, s1=400, s2=600, threshold=20, save=False):
        imaps=[]
        for i in range(self.tDim):
            if z==None:
                #process entire stack at each timepoint
                imap=self.construct3D(t=i, s1=s1, s2=s2, threshold=20, save=save)
            else:
                #process only specified z-plane
                imap=self.construct2D(z=z, t=i,  s1=s1, s2=s2, threshold=20, specImg=False, show=False, save=False)
            imaps.append(imap)
        return np.array(imaps)
    
    def constructSpec3D(self, t=0, s1=400, s2=600, threshold=20, save=False, show=False):
        imaps=[]
        maxmaps=[]
        for z in range(self.zDim):
            imap, maxmap=self.construct2D(z, t,  s1=s1, s2=s2, threshold=20, specImg=True, show=show, save=False)
            imaps.append(imap)
            maxmaps.append(maxmap)
        return np.array(imaps), np.array(maxmaps)
    
    def constructSpecT(self, z=None, s1=400, s2=600, threshold=20, save=False):
        imaps=[]
        maxmaps=[]
        for i in range(self.tDim):
            if z==None:
                #process entire stack at each timepoint
                imap, maxmap=self.constructSpec3D(t=i, s1=s1, s2=s2, threshold=20, save=save)
            else:
                #process only specified z-plane
                imap, maxmap=self.construct2D(z=z, t=i,  s1=s1, s2=s2, threshold=20, specImg=True, show=False, save=False)
            imaps.append(imap)
            maxmaps.append(maxmap)
        return np.array(imaps), np.array(maxmaps)
    #=========================== Fitting functions ================================
    def fitSpectrum(self, spec, tol=50, tol2=100, height=15, plot=True):
        peaks, heights=utils.fitSpectrumWL(spec, self.xaxis, tol=tol, tol2=tol2, height=height)
        if plot:
            plt.figure(figsize=(6, 2))
            plt.plot(self.xaxis, spec)
            plt.scatter(peaks, heights, color='red')
            plt.xlim(min(self.xaxis), max(self.xaxis))
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Counts (a.u.)')
            plt.show()
        return peaks, heights
    
    def getSingleSpectrum(self, x, y, z, t, save=False, pathSave=None):
        path=utils.formatString(self.directory, self.name, t*self.zDim+z)
        rawDat=imageio.imread(path)
        spec=utils.reshape(rawDat[y*self.xDim+x])
        if save:
            if pathSave==None: print('Please specify the file path')
            else: np.savetxt(pathSave, np.transpose(np.array([self.xaxis, spec])))
        return spec
    
    def getSingleSpectrumIntegrated(self, xc, yc, z, t, tol=5, thresh=15, save=False, pathSave=None):
        path=utils.formatString(self.directory, self.name, t*self.zDim+z)
        rawDat=imageio.imread(path)
        specTotal=np.zeros(len(self.xaxis))
        count=0
        for x in range(xc-tol, xc+tol+1):
            for y in range(yc-tol, yc+tol+1):
                spec=utils.reshape(rawDat[y*self.xDim+x])
                if np.max(spec)>thresh:
                    specTotal+=spec
                    count+=1
        specTotal/=count
        if save:
            if pathSave==None: print('Please specify the file path')
            else: np.savetxt(pathSave, np.transpose(np.array([self.xaxis, specTotal])))
        return specTotal
    
    def getSingleSpectrumCCF(self, x, y, z, t, save=False, pathSave=None):
        path=utils.formatString(self.directory, self.name, t*self.zDim*self.yDim+self.yDim*z+y)
        rawDat=imageio.imread(path)
        spec=utils.reshape(rawDat[x])
        if save:
            if pathSave==None: print('Please specify the file path')
            else: np.savetxt(pathSave, np.transpose(np.array([self.xaxis, spec])))
        return spec

    def fitSingleSpectrum(self, x, y, z, t, tol=50, height=15, plot=True):
        spec=self.getSingleSpectrum(x, y, z, t)
        peaks, heights=self.fitSpectrum(spec, tol=tol, height=height, plot=plot)
        return peaks
    
    def fitSpectraToFiles(self, tol=100, tol2=100, height=10, thresh=70, noPeaks=4, t0=0, t1=None, z0=0, z1=None):
        if t1==None: t1=self.tDim
        if z1==None: z1=self.zDim
        for t in range(t0, t1):
            for z in range(z0, z1):
                path=utils.formatString(self.directory, self.name, t*self.zDim+z)
                rawDat=imageio.imread(path)
                frameOutput=[]
                for y in range(self.yDim):
                    for x in range(self.xDim):
                        spec=utils.reshape(rawDat[y*self.xDim+x])
                        if max(spec)>thresh:
                            peaks, heights=self.fitSpectrum(spec, tol=tol, tol2=tol2, height=height, plot=False)
                            #peaks2=sorted(peaks[0:noPeaks])
                            if len(peaks)>=noPeaks:   
                                peaks.insert(0, x)
                                peaks.insert(1, y)
                                peaks.insert(2, max(heights))
                                frameOutput.append(peaks[0:noPeaks+3])
                np.savetxt(self.directory+'fitted\\Fit_t{}_z{}_'.format(t, z)+self.name+'.txt', np.transpose(np.array(frameOutput)))
                
                
    def fitSpectraToFilesASC(self, tol=100, tol2=100, height=10, thresh=70, noPeaks=4, t0=0, t1=None, z0=0, z1=None, plot=False):
        if t0==None: t0=1
        if t1==None: t1=self.tDim+1
        if z1==None: z1=self.zDim
        rawDat=np.loadtxt(self.directory+self.name)
        frameOutput=[]
        for t in range(t0, t1):
            spec=rawDat[:,t]
            if max(spec)>thresh:
                peaks, heights=self.fitSpectrum(spec, tol=tol, tol2=tol2, height=height, plot=plot)
                #peaks2=sorted(peaks[0:noPeaks])
                if len(peaks)>=noPeaks:   
                    peaks.insert(0, t-1)
                    peaks.insert(1, t-1)
                    peaks.insert(2, max(heights))
                    frameOutput.append(peaks[0:noPeaks+3])
        np.savetxt(self.directory+'fitted\\Fit_t{}_z{}_'.format(0, 0)+self.name+'.txt', np.transpose(np.array(frameOutput)))

    def fitSinglePeaks(self, thresh=25, show=True):
        output=[]
        if show:
            plt.figure(figsize=(6,12))
        for t in range(self.tDim):
            for z in range(self.zDim):
                for y in range(self.yDim):
                    path=utils.formatString(self.directory, self.name, t*self.zDim*self.yDim+self.yDim*z+y)
                    data=imageio.imread(path)
                    for x in range(self.xDim):
                        spec=utils.reshape(data[x])
                        I=np.max(spec)
                        if I>thresh and I<4050:
                            result=np.zeros(7)
                            p00=self.xaxis[np.min(np.where(spec==max(spec)))]
                            popt, pcov = curve_fit(utils.gaussian, self.xaxis, spec, p0=[100, 0.1, p00, 15])
                            if show: plt.scatter(t*60*self.tStep, popt[2], color='pink', s=1)
                            result[0]=t
                            result[1]=z
                            result[2]=x
                            result[3]=y
                            result[4]=(np.max(data[x]))
                            result[6]=popt[1]
                            result[5]=popt[2]
                            output.append(result)
        if show:
            plt.ylim(680, 725)
            plt.xlabel('time (s)')
            plt.title(self.name)
            plt.show()
        return np.array(output)
            
    
    def fitSpectraToFilesCCF(self, tol=100, tol2=100, height=10, thresh=70, noPeaks=4, t0=0, t1=None, z0=0, z1=None):
        if t1==None: t1=self.tDim
        if z1==None: z1=self.zDim
        for t in range(t0, t1):
            for z in range(z0, z1):
                frameOutput=[]
                for y in range(self.yDim):
                    path=utils.formatString(self.directory, self.name,  t*self.zDim*self.yDim+self.yDim*z+y)
                    rawDat=imageio.imread(path)
                    for x in range(self.xDim):
                        spec=utils.reshape(rawDat[x])
                        if max(spec)>thresh:
                            peaks, heights=self.fitSpectrum(spec, tol=tol, tol2=tol2, height=height, plot=False)
                            peaks2=sorted(peaks[0:noPeaks])
                            if len(peaks)>=noPeaks:   
                                peaks2.insert(0, x)
                                peaks2.insert(1, y)
                                peaks2.insert(2, max(heights))
                                frameOutput.append(peaks2[0:noPeaks+3])
                np.savetxt(self.directory+'fitted\\Fit_t{}_z{}_'.format(t, z)+self.name+'.txt', np.transpose(np.array(frameOutput)))
    
    #=========================== Displaying functions ================================           
    def cmapResultsFromFiles(self,  axLim0, axLim1, maxRes=10**(-11), param='n', t0=0, t1=None, z0=0, z1=None, fitted=True):
        if t1==None: t1=self.tDim
        if z1==None: z1=self.zDim
        outputMap, iMap = self.ResultsFromFiles(maxRes=maxRes, param=param, t0=t0, t1=t1, z0=z0, z1=z1, fitted=fitted)
        redAll, greenAll, blueAll = utils.VAtoRGBA(outputMap, iMap, axLim0, axLim1)
        return np.array(redAll), np.array(greenAll), np.array(blueAll)  
    
    def ResultsFromFiles(self, maxRes=10**(-11), param='n', t0=0, t1=None, z0=0, z1=None, fitted=True):
        if t1==None: t1=self.tDim
        if z1==None: z1=self.zDim
        outputMap=np.zeros((t1-t0, z1-z0, self.xDim, self.yDim))
        iMap=np.zeros((t1-t0, z1-z0, self.xDim, self.yDim))
        for t in range(t0, t1):
            for z in range(z0, z1):
                try:
                    if fitted: results=np.transpose(np.loadtxt(self.directory+'fitted\\AE_Fit_t{}_z{}_'.format(t, z)+self.name+'.txt'))
                    if not fitted: results=np.transpose(np.loadtxt(self.directory+'fitted\\Fit_t{}_z{}_'.format(t, z)+self.name+'.txt'))
                except:
                    print('Fit_t{}_z{}_'.format(t, z)+self.name+' could not be read')
                if len(np.shape(results))==1:
                    results=[results]
                for result in results:
                    if not fitted or np.abs(result[-1])<maxRes:
                        x=int(result[0])
                        y=int(result[1])
                        I=result[2]
                        iMap[t][z][x][y]=I
                        if not fitted:
                            outputMap[t][z][x][y]=result[3]
                        elif param=='n':
                            outputMap[t][z][x][y]=result[-2]
                        elif param=='d':
                            outputMap[t][z][x][y]=result[-3]
                        elif param=='peak':
                            outputMap[t][z][x][y]=result[-6]
                            iMap[t][z][x][y]=result[-5]
        return outputMap, iMap        
    
    def ResultsFromFilesASC(self, maxRes=10**(-11), param='n', t0=0, t1=None, z0=0, z1=None, fitted=True):
        if t1==None: t1=self.tDim
        if z1==None: z1=self.zDim
        outputMap=np.zeros(self.tDim)
        iMap=np.zeros(self.tDim)
        for t in range(t0, t1):
            for z in range(z0, z1):
                try:
                    if fitted: results=np.transpose(np.loadtxt(self.directory+'fitted\\AE_Fit_t{}_z{}_'.format(t, z)+self.name+'.txt'))
                    if not fitted: results=np.transpose(np.loadtxt(self.directory+'fitted\\Fit_t{}_z{}_'.format(t, z)+self.name+'.txt'))
                except:
                    print('Fit_t{}_z{}_'.format(t, z)+self.name+' could not be read')
                if len(np.shape(results))==1:
                    results=[results]
                for result in results:
                    if not fitted or np.abs(result[-1])<maxRes:
                        t=int(result[1])
                        I=result[2]
                        iMap[t]=I
                        if not fitted:
                            outputMap[t]=result[3]
                        elif param=='n':
                            outputMap[t]=result[-2]
                        elif param=='d':
                            outputMap[t]=result[-3]
                        elif param=='peak':
                            outputMap[t]=result[-6]
                            iMap[t]=result[-5]
                        if param=='peak' and not fitted:
                            outputMap[t]=result[-2]
                            iMap[t]=result[-1]
        return outputMap, iMap   