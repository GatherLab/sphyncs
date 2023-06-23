from acquisitionExtended import * #processing data class and dependencies
from pandas import DataFrame, Series, read_csv  # for convenience
import pims
import trackpy as tp
from zipfile import ZipFile
import os

class disk: #container to store trajectory information of disks
    def __init__(self, timesteps, partNo, tStep):
        self.timesteps=timesteps
        self.partNo=partNo
        self.x=np.zeros(timesteps)
        self.y=np.zeros(timesteps)
        self.z=np.zeros(timesteps)
        self.peaks=np.zeros(timesteps)
        self.time=np.linspace(0, timesteps*tStep, timesteps)
    
    def fitAllWL(self, acq , LPsize=[1,3,3], plot=True, thresh=20): #xdim, zdim, timesteps,
        allpeaks=[]
        spec0=np.zeros(len(xaxis))
        if plot:
            plt.figure(figsize=(12,4))
        for t in range(acq.tDim):
            diskX=int(round(self.x[t], 0))
            diskY=int(round(self.y[t], 0))
            diskZ=int(self.z[t])
            if diskX+diskY > 0: #means there is an actual disk
                for Z in range(diskZ-LPsize[0], diskZ+1+LPsize[0]):
                    path=utils.formatString(acq.directory, acq.name, t*acq.zDim+Z)
                    try:
                        rawdata=imageio.imread(path)
                        for X in range(diskX-LPsize[1], diskX+LPsize[1]):
                            for Y in range(diskY-LPsize[2], diskY+LPsize[2]):
                                spec=reshape(rawdata[acq.xDim*Y+X])
                                if np.max(spec)>thresh:
                                    spec0+=spec
                        del rawdata
                    except:
                        peak=-1
                if plot:
                    r, g, b = utils.VAtoRGBA(t, 1/3, 0,-acq.tDim[-1]*0.214, acq.tDim[-1]+(acq.tDim[-1]*0.385))
                    r, g, b = utils.normArr(r, g, b)
                    plt.plot(xaxis, spec0/np.max(spec0), color=(r, g, b))
                if np.max(spec0)>thresh:
                    peaks, heights=utils.fitSpectrumWL(spec, acq.xaxis, height=0.8*np.max(spec)) #only brightest peak (singlemode)
                    peak=peaks[0]
                    if peak!=-1:
                        allpeaks.append(peak)
                    self.peaks[t]=peak
        if plot:
            peakAvg=np.average(allpeaks)
            plt.title(self.partNo)
            plt.xlabel(r'$\lambda$ (nm)')
            plt.xlim(peakAvg-2, peakAvg+2)
            plt.show()
        
    def plotTraj(self, acq, ax0='x'):
        validData=self.x!=0
        if ax0=='x':
            myX=self.x[validData]
            upLim=acq.xDim
        elif ax0=='y':
            myX=self.y[validData]
            upLim=acq.yDim
        else:
            myX=self.z[validData]
            upLim=acq.zDim
        plt.figure()
        time=self.time[validData]
        for i in range(1, len(time)): 
            r, g, b = utils.VAtoRGBA(time[i], 1/2, -self.time[-1]*0.3, self.time[-1]+(self.time[-1]*0.4))
            r, g, b = utils.normArr(r, g, b) #color re-scaled 
            plt.plot([time[i-1], time[i]], [myX[i-1], myX[i]], color=(r,g,b))
        plt.xlabel('time')
        plt.ylabel(ax0)
        plt.ylim(0, upLim)
        plt.show()
        
    def plotXY(self, acq, colorCode='wL', nBin=[1,1,1,1]):
        validData=self.x!=0
        x=(self.x[validData]) 
        y=(self.y[validData])
        z=(self.z[validData])
        peaks=(self.peaks[validData])
        time=self.time[validData]
        for i in range(1, len(time)): 
            if colorCode=='time':
                r, g, b = utils.VAtoRGBA(time[i], 1/2, -self.time[-1]*0.3, self.time[-1]+(self.time[-1]*0.4)) #color re-scaled 
                r, g, b = utils.normArr(r, g, b)
            elif colorCode=='wL':
                r, g, b = utils.VAtoRGBA(peaks[i], 1/2, min(acq.xaxis), max(acq.xaxis))
                r, g, b = utils.normArr(r, g, b)
            else: ###huh
                r, g, b = VAtoRGBA((zcorr[i]+zcorr[i-1])/2, 1, 8, 23)
                r, g, b = utils.normArr(r, g, b)
            try:
                plt.plot([x[i-1], x[i]], [y[i-1], y[i]], color=(r,g,b))
            except:
                print(time[i])
        plt.xlim(0, acq.xDim/nBin[0])
        plt.ylim(0, acq.yDim/nBin[1]) 

class trackDataND:    
    def __init__(self, linkedD, acq, nBin=[1,1,1,4],name='A0-'):#, xaxis=np.linspace(719.95, 681.43, 150), dim=[375, 375, 40, 48]):
        self.disks={}
        self.xdim=acq.xDim
        self.ydim=acq.yDim
        self.zdim=acq.zDim
        self.timesteps=acq.tDim*acq.tStep    
        self.nBin=nBin
        for i in range(len(linkedD)):
            if type(linkedD)==np.ndarray:
            #construction from processed array
                part=linkedD[i][0]
                x=linkedD[i, 0*self.timesteps+1:1*self.timesteps+1]        
                y=linkedD[i, 1*self.timesteps+1:2*self.timesteps+1]        
                z=linkedD[i, 2*self.timesteps+1:3*self.timesteps+1] 
                t=linkedD[i, 3*self.timesteps+1:4*self.timesteps+1]
                peaks=linkedD[i, 4*self.timesteps+1:5*self.timesteps+1]
                self.disks[part]=disk(dim[3], part)
                (self.disks[part].x)=x
                (self.disks[part].y)=y            
                (self.disks[part].z)=z
                (self.disks[part].t)=t
                (self.disks[part].peaks)=peaks
            elif type(linkedD)==pd.core.frame.DataFrame:
            #construction from tracking dataframe
                part=name+str(linkedD.at[i, 'particle'])
                x=linkedD.at[i, 'x2']        
                y=linkedD.at[i, 'x3']        
                z=linkedD.at[i, 'x1']
                peak=acq.xaxis[int(linkedD.at[i, 'x0'])*self.nBin[3]]
                t=linkedD.at[i, 'frame']
                if (part in self.disks.keys()):
                    (self.disks[part].x)[t]=x
                    (self.disks[part].y)[t]=y            
                    (self.disks[part].z)[t]=z
                    (self.disks[part].peaks)[t]=peak
                else:
                    self.disks[part]=disk(acq.tDim, part, acq.tStep)
                    (self.disks[part].x)[t]=x
                    (self.disks[part].y)[t]=y            
                    (self.disks[part].z)[t]=z
                    (self.disks[part].peaks)[t]=peak
                
    def mergeWith(self, disks2):
        for key in disks2.disks:
            self.disks[key]=disks2.disks[key]
            
    def fitAllDisks(self, LPsize=[1,3,3], plot0=True, minLen=30, thresh=20):
        for key in self.disks:
            xval=list(filter(lambda num: num != 0, self.disks[key].x))
            if len(xval)>minLen:
                self.disks[key].fitAllWL(self.xdim, self.zdim, self.timesteps, LPsize=[1,3,3], plot=plot0, thresh=thresh)
            
    def plotTrajs(self, acq, minLen=47, leg=False, colorCode='wL' , path='temp.png', nBin=[1,1,1,1]):
        count=0
        plt.figure(figsize=(5, 5))
        for key in self.disks:
            xval=list(filter(lambda num: num != 0, self.disks[key].x))
            if len(xval)>minLen:
                (self.disks[key]).plotXY(acq, colorCode=colorCode)
        plt.xlim(0, acq.xDim/nBin[0])
        plt.ylim(0, acq.yDim/nBin[1]) 
        plt.xticks([])
        plt.yticks([])
        plt.savefig(path)
        plt.show()