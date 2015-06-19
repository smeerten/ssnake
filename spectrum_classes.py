import matplotlib
import numpy as np
import scipy.optimize
import scipy.signal
import scipy.ndimage
import copy
from spectrumFrame import Plot1DFrame

#########################################################################
#the generic data class
class Spectrum(object):
    def __init__(self, data, freq, sw , spec=None, wholeEcho=None, ref=None):
        self.dim = len(data.shape)                    #number of dimensions
        self.data = np.array(data,dtype=complex)      #data of dimension dim
        self.freq = freq                              #array of center frequency (length is dim, MHz)
        self.sw = sw                                  #array of sweepwidths
        if spec is None:
            self.spec=[0]*self.dim
        else:
            self.spec = spec                              #int array of length dim where 0 = time domain, 1 = complex spectral domain and 2= real spectral domain
        if wholeEcho is None:
            self.wholeEcho = [False]*self.dim
        else:
            self.wholeEcho = wholeEcho                    #boolean array of length dim where True indicates a full Echo
        if ref is None:
            self.ref = self.freq
        else:
            self.ref = ref
        self.xaxArray = [[] for i in range(self.dim)]
        self.resetXax()

    def resetXax(self,axes=None):
        if axes is not None:
            val=[axes]
        else:
            val=range(self.dim)
        for i in val:
            if self.spec[i]==0:
                self.xaxArray[i]=np.arange(self.data.shape[i])/(self.sw[i])
            elif self.spec[i]==1:
                self.xaxArray[i]=np.fft.fftshift(np.fft.fftfreq(self.data.shape[i],1.0/self.sw[i]))

    def real(self):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.real())
        self.data = np.real(self.data)
        return returnValue

    def imag(self):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.imag())
        self.data = np.imag(self.data)
        return returnValue

    def abs(self):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.abs())
        self.data = np.abs(self.data)
        return returnValue

    def matrixManip(self, pos1, pos2, axes, which):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.matrixManip(pos1,pos2,axes,which))
        minPos = min(pos1,pos2)
        maxPos = max(pos1,pos2)
        slicing = (slice(None),) * axes + (slice(minPos,maxPos),) + (slice(None),)*(self.dim-1-axes)
        if which == 0:
            self.data = np.sum(self.data[slicing],axis=axes)
        elif which == 1:
            self.data = np.amax(self.data[slicing],axis=axes)
        elif which == 2:
            self.data = np.amin(self.data[slicing],axis=axes)
        self.dim = self.dim - 1
        self.freq = np.delete(self.freq,axes)
        self.sw = np.delete(self.sw,axes)
        self.spec = np.delete(self.spec,axes)
        self.wholeEcho = np.delete(self.wholeEcho,axes)
        self.xaxArray = np.delete(self.xaxArray,axes)
        return returnValue

    def getRegion(self, pos1, pos2,axes):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.getRegion(axes,pos1,pos2))
        minPos = min(pos1,pos2)
        maxPos = max(pos1,pos2)
        slicing = (slice(None),) * axes + (slice(minPos,maxPos),) + (slice(None),)*(self.dim-1-axes)
        self.data = self.data[slicing]
        self.xaxArray[axes] = self.xaxArray[axes][slice(minPos,maxPos)] #what to do with sw?
        return returnValue

    def flipLR(self, axes):
        slicing = (slice(None),) * axes + (slice(None,None,-1),) + (slice(None),)*(self.dim-1-axes)
        self.data = self.data[slicing]
        return lambda self: self.flipLR(axes)
    
    def hilbert(self,axes):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.hilbert(axes))
        self.data = scipy.signal.hilbert(np.real(self.data), axis=axes)
        return returnValue

    def setPhase(self, phase0, phase1, axes):
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.data.shape[axes],1.0/self.sw[axes]))*phase1*1j)
        if self.spec[axes]==0:
            self.fourier(axes,tmp=True)
            self.data=self.data*np.exp(phase0*1j)
            for i in range(self.data.shape[axes]):
                slicing = (slice(None),) * axes + (i,) + (slice(None),)*(self.dim-1-axes)
                self.data[slicing]=self.data[slicing]*vector[i]
            self.fourier(axes,tmp=True)
        else:
            self.data=self.data*np.exp(phase0*1j)
            for i in range(self.data.shape[axes]):
                slicing = (slice(None),) * axes + (i,) + (slice(None),)*(self.dim-1-axes)
                self.data[slicing]=self.data[slicing]*vector[i]
        return lambda self: self.setPhase(-phase0,-phase1,axes)

    def apodize(self,lor,gauss, cos2, hamming, shift, shifting, shiftingAxes, axes):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.apodize(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes,axes))
        axLen = self.data.shape[axes]
        t=np.arange(0,axLen)/self.sw[axes]
        if shifting != 0.0:
            for j in range(self.data.shape[shiftingAxes]):
                shift1 = shift + shifting*j
                t2 = t - shift1
                x=np.ones(axLen)
                if lor is not None:
                    x=x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x=x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x=x*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw[axes]/axLen+np.linspace(0,0.5*np.pi,axLen)))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw[axes]/axLen+np.linspace(0,np.pi,axLen))))

                if self.spec[axes] > 0:
                    self.fourier(axes,tmp=True)
                    for i in range(self.data.shape[axes]):
                        if axes < shiftingAxes:
                            slicing = (slice(None),) * axes + (i,) + (slice(None),)*(shiftingAxes-1-axes) + (j,) + (slice(None),)*(self.dim-2-shiftingAxes)
                        else:
                            slicing = (slice(None),) * shiftingAxes + (j,) + (slice(None),)*(axes-1-shiftingAxes) + (i,) + (slice(None),)*(self.dim-2-axes)
                        self.data[slicing]=self.data[slicing]*x[i]
                    self.fourier(axes,tmp=True)
                else:
                    for i in range(self.data.shape[axes]):
                        if axes < shiftingAxes:
                            slicing = (slice(None),) * axes + (i,) + (slice(None),)*(shiftingAxes-1-axes) + (j,) + (slice(None),)*(self.dim-2-shiftingAxes)
                        else:
                            slicing = (slice(None),) * shiftingAxes + (j,) + (slice(None),)*(axes-1-shiftingAxes) + (i,) + (slice(None),)*(self.dim-2-axes)
                        self.data[slicing]=self.data[slicing]*x[i]
        else:
            t2 = t - shift
            x=np.ones(axLen)
            if lor is not None:
                x=x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x=x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw[axes]/axLen+np.linspace(0,0.5*np.pi,axLen)))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw[axes]/axLen+np.linspace(0,np.pi,axLen))))
            if self.spec[axes] > 0:
                self.fourier(axes,tmp=True)
                for i in range(self.data.shape[axes]):
                    slicing = (slice(None),) * axes + (i,) + (slice(None),)*(self.dim-1-axes)
                    self.data[slicing]=self.data[slicing]*x[i]
                self.fourier(axes,tmp=True)
            else:
                for i in range(self.data.shape[axes]):
                    slicing = (slice(None),) * axes + (i,) + (slice(None),)*(self.dim-1-axes)
                    self.data[slicing]=self.data[slicing]*x[i]
        return returnValue

    def setFreq(self,freq,sw,axes):
        oldFreq = self.freq[axes]
        oldSw = self.sw[axes]
        self.freq[axes]=freq
        self.sw[axes]=sw
        self.resetXax(axes)
        return lambda self: self.setFreq(oldFreq,oldSw,axes)

    def setSize(self,size,axes):
        if axes>self.dim:
            print("axes bigger than dim in setsize")
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.setSize(size,axes))
        if size > self.data.shape[axes]:
            if self.wholeEcho[axes]:
                tmpdata = np.array_split(self.data,2,axes)
                self.data = np.concatenate((np.pad(tmpdata[0],[(0,0)]*axes+[(0,size-self.data.shape[axes])]+[(0,0)]*(self.dim-axes-1),'constant',constant_values=0),tmpdata[1]),axes)
            else:
                self.data = np.pad(self.data,[(0,0)]*axes+[(0,size-self.data.shape[axes])]+[(0,0)]*(self.dim-axes-1),'constant',constant_values=0)
        else:
            if self.wholeEcho[axes]:
                slicing1  = (slice(None),) * axes + (slice(0,np.ceil(size/2.0)),) + (slice(None),)*(self.dim-1-axes)
                slicing2  = (slice(None),) * axes + (slice(size/2,None),) + (slice(None),)*(self.dim-1-axes)
                self.data = np.concatenate((self.data[slicing1],self.data[slicing2]),axes)
            else:
                slicing = (slice(None),) * axes + (slice(0,size),) + (slice(None),)*(self.dim-1-axes)
                self.data = self.data[slicing]
        self.dim = len(self.data.shape)
        self.resetXax(axes)
        return returnValue

    def changeSpec(self,val,axes):
        oldVal = self.spec[axes]
        self.spec[axes] = val
        self.resetXax(axes)
        return lambda self: self.changeSpec(oldVal,axes)

    def swapEcho(self,idx,axes):
        slicing1=(slice(None),) * axes + (slice(None,idx),) + (slice(None),)*(self.dim-1-axes)
        slicing2=(slice(None),) * axes + (slice(idx,None),) + (slice(None),)*(self.dim-1-axes)
        self.data = np.concatenate((self.data[slicing2],self.data[slicing1]),axes)
        self.wholeEcho[axes] = not self.wholeEcho[axes]
        return lambda self: self.swapEcho(-idx,axes)
            
    def shiftData(self,shift,axes):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.shiftData(shift,axes))
        self.data = np.roll(self.data,shift,axes)
        if shift < 0:
            slicing = (slice(None),) * axes + (slice(shift,None),) + (slice(None),)*(self.dim-1-axes)
        else:
            slicing = (slice(None),) * axes + (slice(None,shift),) + (slice(None),)*(self.dim-1-axes)
        self.data[slicing]=self.data[slicing]*0
        return returnValue

    def dcOffset(self,offset):
        self.data = self.data+offset
        return lambda self: self.dcOffset(-offset)

    def fourier(self, axes,tmp=False):
        if axes>self.dim:
            print("axes bigger than dim in fourier")
        if self.spec[axes]==0:
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None),) * axes + (0,) + (slice(None),)*(self.dim-1-axes)
                self.data[slicing]= self.data[slicing]*0.5
            self.data=np.fft.fftshift(np.fft.fftn(self.data,axes=[axes]),axes=axes)
            self.spec[axes]=1
            
        else:
            self.data=np.fft.ifftn(np.fft.ifftshift(self.data,axes=axes),axes=[axes])
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None),) * axes + (0,) + (slice(None),)*(self.dim-1-axes)
                self.data[slicing]= self.data[slicing]*2.0
            self.spec[axes]=0
        self.resetXax(axes)
        return lambda self: self.fourier(axes)

    def fftshift(self, axes, inv=False):
        if axes>self.dim:
            print("axes bigger than dim in fourier")
        if inv:
            self.data=np.fft.ifftshift(self.data,axes=[axes])
        else:
            self.data=np.fft.fftshift(self.data,axes=axes)
        return lambda self: self.fftshift(axes,not(inv))

    def shear(self, shear, axes, axes2):
        if self.dim < 2:
            print("The data does not have enough dimensions for a shearing transformation")
            return None
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.shear(shear,axes, axes2))
        shearMatrix = np.identity(self.dim)
        shearMatrix[axes,axes2]=shear
        self.data = scipy.ndimage.interpolation.affine_transform(np.real(self.data),shearMatrix,mode='wrap') + 1j*scipy.ndimage.interpolation.affine_transform(np.imag(self.data),shearMatrix,mode='wrap')
        return returnValue
    
    def getSlice(self,axes,locList):
        return (self.data[tuple(locList[:axes])+(slice(None),)+tuple(locList[axes:])],self.freq[axes],self.sw[axes],self.spec[axes],self.wholeEcho[axes],self.xaxArray[axes],self.ref[axes])

    def getBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None):
        stackSlice = slice(stackBegin, stackEnd, stackStep)
        if axes == axes2:
            print("First and second axes are the same")
            return
        elif axes < axes2:
            return (np.transpose(self.data[tuple(locList[:axes])+(slice(None),)+tuple(locList[axes:axes2-1])+(stackSlice,)+tuple(locList[axes2-1:])]),self.freq[axes],self.sw[axes],self.spec[axes],self.wholeEcho[axes],self.xaxArray[axes],self.xaxArray[axes2],self.ref[axes],self.ref[axes2])
        elif axes > axes2:
            return (self.data[tuple(locList[:axes2])+(stackSlice,)+tuple(locList[axes2:axes-1])+(slice(None),)+tuple(locList[axes-1:])],self.freq[axes],self.sw[axes],self.spec[axes],self.wholeEcho[axes],self.xaxArray[axes],self.xaxArray[axes2],self.ref[axes],self.ref[axes2])

    def restoreData(self,copyData,returnValue): # restore data from an old copy for undo purposes
        self.data = copyData.data
        self.dim = len(self.data.shape)                    #number of dimensions
        self.freq = copyData.freq                              #array of center frequency (length is dim, MHz)
        self.sw = copyData.sw                                  #array of sweepwidths
        self.spec = copyData.spec
        self.wholeEcho = copyData.wholeEcho
        self.xaxArray = copyData.xaxArray
        self.ref = copyData.ref
        return returnValue


#########################################################################################################
#the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing
class Current1D(Plot1DFrame):
    def __init__(self, root, data, axes=None, locList=None, plotType=0, axType=1,ppm=False):
        Plot1DFrame.__init__(self,root)
        self.xax = None               #x-axis
        self.data = data              #the actual spectrum instance
        self.freq = None              #frequency of the slice 
        self.sw = None                #x-data display
        self.data1D = None            #the data1D
        self.spec = None              #boolean where False=time domain and True=spectral domain
        self.wholeEcho = None
        self.ref = None               #reference frequency
        self.ppm = ppm                #display frequency as ppm
        if axes is None:
            self.axes = len(self.data.data.shape)-1
        else:
            self.axes = axes              #dimension of which the data is displayed
        if locList is None:
            self.resetLocList()
        else:
            self.locList = locList    
        self.plotType = plotType
        self.axType = axType
        self.upd()   #get the first slice of data
        self.plotReset() #reset the axes limits
        self.showFid() #plot the data
        
    def upd(self): #get new data from the data instance
        if (len(self.locList)+1) != self.data.dim:
            self.resetLocList()
        updateVar = self.data.getSlice(self.axes,self.locList)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.sw = updateVar[2]
        self.spec = updateVar[3]
        self.wholeEcho = updateVar[4]
        self.xax=updateVar[5]
        self.ref=updateVar[6]

    def setSlice(self,axes,locList): #change the slice 
        axesSame = True
        if self.axes != axes:
            axesSame = False
            self.axes = axes
        self.locList = locList
        self.upd()
        if not axesSame:
            self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-1)
        
    def setPhaseInter(self, phase0in, phase1in): #interactive changing the phase without editing the actual data
        phase0=float(phase0in)
        phase1=float(phase1in)
        if self.spec==0:
            tmpdata=self.fourierLocal(self.data1D,0)
            tmpdata=tmpdata*np.exp(phase0*1j)
            tmpdata=tmpdata*np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata),1.0/self.sw))*phase1*1j)
            tmpdata=self.fourierLocal(tmpdata,1)
        else:
            tmpdata=self.data1D*np.exp(phase0*1j)
            tmpdata=tmpdata*np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata),1.0/self.sw))*phase1*1j)
        self.showFid(tmpdata)

    def applyPhase(self, phase0, phase1):# apply the phase to the actual data
        phase0=float(phase0)
        phase1=float(phase1)
        returnValue = self.data.setPhase(phase0,phase1,self.axes)
        self.upd()
        self.showFid()
        return returnValue

    def fourier(self): #fourier the actual data and replot
        returnValue = self.data.fourier(self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        return returnValue

    def fftshift(self,inv=False): #fftshift the actual data and replot
        returnValue = self.data.fftshift(self.axes,inv)
        self.upd()
        self.showFid()
        return returnValue

    def fourierLocal(self, fourData, spec): #fourier the local data for other functions
        ax = len(self.data1D.shape)-1
        if spec==0:
            fourData=np.fft.fftshift(np.fft.fftn(fourData,axes=[ax]),axes=ax)
        else:
            fourData=np.fft.ifftn(np.fft.ifftshift(fourData,axes=ax),axes=[ax])
        return fourData

    def apodPreview(self,lor=None,gauss=None, cos2=None, hamming=None ,shift=0.0,shifting=0.0,shiftingAxes=None): #display the 1D data including the apodization function
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                print('shiftingAxes cannot be axes')
                return
            elif shiftingAxes < self.axes:
                shift += shifting*self.locList[shiftingAxes]
            else:
                shift += shifting*self.locList[shiftingAxes-1]
        t=np.arange(0,len(self.data1D))/(self.sw)
        t2=t-shift
        x=np.ones(len(self.data1D))
        if lor is not None:
            x=x*np.exp(-lor*abs(t2))
        if gauss is not None:
            x=x*np.exp(-(gauss*t2)**2)
        if cos2 is not None:
            x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,0.5*np.pi,len(self.data1D))))**2)
        if hamming is not None:
            alpha = 0.53836 # constant for hamming window
            x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
        if self.wholeEcho:
            x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
        a=self.fig.gca()
        a.cla()
        y = self.data1D
        if self.spec ==1:
            y=np.fft.ifft(np.fft.ifftshift(y))
            y= y*x
            y=np.fft.fftshift(np.fft.fft(y))
        else:
            y= y*x
        if self.spec==0:
            if self.plotType==0:
                self.showFid(y,[t],[x*max(np.real(self.data1D))],['g'],old=True)
            elif self.plotType==1:
                self.showFid(y,[t],[x*max(np.imag(self.data1D))],['g'],old=True)
            elif self.plotType==2:
                self.showFid(y,[t],[x*max(max(np.real(self.data1D)),max(np.imag(self.data1D)))],['g'],old=True)
            elif self.plotType==3:
                self.showFid(y,[t],[x*max(np.abs(self.data1D))],['g'],old=True)
        else:
            self.showFid(y)

    def applyApod(self,lor=None,gauss=None,cos2=None,hamming=None,shift=0.0,shifting=0.0,shiftingAxes=0): #apply the apodization to the actual data
        returnValue = self.data.apodize(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes,self.axes)
        self.upd() 
        self.showFid()     
        return returnValue

    def setFreq(self,freq,sw): #set the frequency of the actual data
        returnValue = self.data.setFreq(freq,sw,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        return returnValue

    def setSizePreview(self,size): #set size only on local data
        if size > len(self.data1D):
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2)
                self.data1D = np.concatenate((np.pad(tmpdata[0],(0,size-len(self.data1D)),'constant',constant_values=0),tmpdata[1]))
            else:
                self.data1D = np.pad(self.data1D,(0,size-len(self.data1D)),'constant',constant_values=0)
        else:
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2)
                self.data1D = np.concatenate((tmpdata[:np.ceil(size/2.0)],tmpdata[size/2:]))
            else:
                self.data1D = self.data1D[:size]
        if self.spec==0:
            self.xax=np.arange(len(self.data1D))/self.sw
        elif self.spec==1:
            self.xax=np.fft.fftshift(np.fft.fftfreq(len(self.data1D),1.0/self.sw)) 
        self.plotReset()
        self.showFid()
        self.upd()

    def applySize(self,size): #set size to the actual data
        returnValue = self.data.setSize(size,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        return returnValue

    def changeSpec(self,val): #change from time to freq domain of the actual data
        returnValue = self.data.changeSpec(val,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        return returnValue

    def applySwapEcho(self,idx):
        returnValue = self.data.swapEcho(idx,self.axes)
        self.upd()
        self.showFid()
        return returnValue

    def setSwapEchoPreview(self,idx):
        self.data1D = np.concatenate((self.data1D[idx:],self.data1D[:idx]))
        self.plotReset()
        self.showFid()
        self.upd()

    def setWholeEcho(self, value):
        if value == 0:
            self.data.wholeEcho[self.axes]=False
            self.wholeEcho = False
        else:
            self.data.wholeEcho[self.axes]=True
            self.wholeEcho = True

    def applyShift(self,shift):
        returnValue = self.data.shiftData(shift,self.axes)
        self.upd()
        self.showFid()
        return returnValue
        
    def setShiftPreview(self,shift):
        tmpData = np.roll(self.data1D,shift)
        if shift<0:
            tmpData[shift:] = tmpData[shift:]*0
        else:
            tmpData[:shift] = tmpData[:shift]*0
        self.showFid(tmpData)

    def dcOffset(self,pos1,pos2):
        minPos = int(min(pos1,pos2))
        maxPos = int(max(pos1,pos2))
        if minPos != maxPos:
            self.showFid(self.data1D-np.mean(self.data1D[minPos:maxPos]))
            
    def applydcOffset(self,pos1,pos2):
        minPos = int(min(pos1,pos2))
        maxPos = int(max(pos1,pos2))
        returnValue = self.data.dcOffset(-np.mean(self.data1D[minPos:maxPos]))
        self.upd()
        self.showFid()
        return returnValue
    
    def integrate(self,pos1,pos2):
        return self.data.matrixManip(pos1,pos2,self.axes,0)

    def maxMatrix(self,pos1,pos2):
        return self.data.matrixManip(pos1,pos2,self.axes,1)
    
    def minMatrix(self,pos1,pos2):
        return self.data.matrixManip(pos1,pos2,self.axes,2)

    def flipLR(self):
        returnValue = self.data.flipLR(self.axes)
        self.upd()
        self.showFid()
        return returnValue
    
    def getRegion(self,pos1,pos2): #set the frequency of the actual data
        returnValue = self.data.getRegion(pos1,pos2,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        return returnValue
    
    def ACMEentropy(self,phaseIn,phaseAll=True):
        phase0=phaseIn[0]
        if phaseAll:
            phase1=phaseIn[1]
        else:
            phase1=0.0
        L = len(self.data1D)
        if self.spec==1:
            x=np.fft.fftshift(np.fft.fftfreq(L,1.0/self.sw))
        if self.spec>0:
            s0 = self.data1D*np.exp(1j*(phase0+phase1*x))
        else:
            s0 = np.fft.fftshift(np.fft.fft(self.data1D))*np.exp(1j*(phase0+phase1*x))
        s2 = np.real(s0)
        ds1 = np.abs((s2[3:L]-s2[1:L-2])/2.0)
        p1 = ds1/sum(ds1)
        p1[np.where(p1 == 0)] = 1
        h1  = -p1*np.log(p1)
        H1  = sum(h1)
        Pfun = 0.0
        as1 = s2 - np.abs(s2)
        sumas   = sum(as1)
        if (np.real(sumas) < 0): 
            Pfun = Pfun + sum(as1**2)/4/L**2
        return H1+1000*Pfun 

    def autoPhase(self,phaseNum):
        if phaseNum == 0:
            phases = scipy.optimize.fmin(func=self.ACMEentropy,x0=[0],args=(False,))
        elif phaseNum == 1:
            phases = scipy.optimize.fmin(func=self.ACMEentropy,x0=[0,0])
        return phases

    def setXaxPreview(self,xax):
        self.xax = xax
        self.plotReset()
        self.showFid()
        self.upd()

    def setXax(self,xax):
        self.data.xaxArray[self.axes]= xax 
        self.upd()
        self.plotReset()
        self.showFid()
        #for now the changing of the axis cannot be undone, because of problems in case of a non-linear axis on operations such as fourier transform etc.
        #doing one of these operations will result in a return to the default axis defined by sw

    def setAxType(self, val):
        oldAxAdd = 0
        if self.spec == 1:
            if self.ppm:
                oldAxAdd = (self.freq-self.ref)/self.ref*1e6
                oldAxMult = 1e6/self.ref
            else:
                oldAxMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            oldAxMult = 1000.0**self.axType
        newAxAdd = 0
        if self.spec == 1:
            if val == 'ppm':
                newAxAdd = (self.freq-self.ref)/self.ref*1e6
                newAxMult = 1e6/self.ref
            else:
                newAxMult = 1.0/(1000.0**val)
        elif self.spec == 0:
            newAxMult = 1000.0**val 
        if val == 'ppm':
            self.ppm = True
        else:
            self.ppm = False
            self.axType = val
        self.xminlim = self.xminlim * newAxMult / oldAxMult + newAxAdd - oldAxAdd 
        self.xmaxlim = self.xmaxlim * newAxMult / oldAxMult + newAxAdd - oldAxAdd 
        self.showFid()

    def hilbert(self):
        returnValue = self.data.hilbert(self.axes)
        self.upd()
        self.showFid()
        return returnValue
    
    def getDisplayedData(self):
        if self.plotType==0:
            return np.real(self.data1D)
        elif self.plotType==1:
            return np.imag(self.data1D)
        elif self.plotType==2:
            return np.real(self.data1D)
        elif self.plotType==3:
            return np.abs(self.data1D)      

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        a=self.fig.gca()
        a.cla()
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if old:
            if (self.plotType==0):
                a.plot(self.xax*axMult+axAdd,np.real(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==1):
                a.plot(self.xax*axMult+axAdd,np.imag(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==2):
                a.plot(self.xax*axMult+axAdd,np.real(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==3):
                a.plot(self.xax*axMult+axAdd,np.abs(self.data1D),c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraX)):
                a.plot(extraX[num]*axMult+axAdd,extraY[num],c=extraColor[num])
        if (self.plotType==0):
            self.line = a.plot(self.xax*axMult+axAdd,np.real(tmpdata),c='b')
        elif(self.plotType==1):
            self.line = a.plot(self.xax*axMult+axAdd,np.imag(tmpdata),c='b')
        elif(self.plotType==2):
            a.plot(self.xax*axMult+axAdd,np.imag(tmpdata),c='r')
            self.line = a.plot(self.xax*axMult+axAdd,np.real(tmpdata),c='b')
        elif(self.plotType==3):
            self.line = a.plot(self.xax*axMult+axAdd,np.abs(tmpdata),c='b')
        a.set_title("TD"+str(self.axes+1))
        if self.spec==0:
            if self.axType == 0:
                a.set_xlabel('Time [s]')
            elif self.axType == 1:
                a.set_xlabel('Time [ms]')
            elif self.axType == 2:
                a.set_xlabel(r'Time [$\mu$s]')
            else:
                a.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                a.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    a.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    a.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    a.set_xlabel('Frequency [MHz]')
                else:
                    a.set_xlabel('User defined')
        else:
            a.set_xlabel('')
        a.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        a.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
        a=self.fig.gca()
        if self.plotType==0:
            miny = min(np.real(self.data1D))
            maxy = max(np.real(self.data1D))
        elif self.plotType==1:
            miny = min(np.imag(self.data1D))
            maxy = max(np.imag(self.data1D))
        elif self.plotType==2:
            miny = min(min(np.real(self.data1D)),min(np.imag(self.data1D)))
            maxy = max(max(np.real(self.data1D)),max(np.imag(self.data1D)))
        elif self.plotType==3:
            miny = min(np.abs(self.data1D))
            maxy = max(np.abs(self.data1D))
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=max(self.xax*axMult+axAdd)
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)

#########################################################################################################
#the class from which the stacked data is displayed, the operations which only edit the content of this class are for previewing
class CurrentStacked(Current1D):
    def __init__(self, root, data, axes=None, axes2=None, locList=None, plotType=0, axType=1,ppm=False ,stackBegin=None, stackEnd=None, stackStep=None):
        self.data = data
        if axes2 is None:
            self.axes2 = len(self.data.data.shape)-2
        else:
            self.axes2 = axes2            #dimension from which the spectra are stacked
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if locList is None:
            self.resetLocList()
        self.spacing = 0
        Current1D.__init__(self,root, data, axes, locList, plotType, axType,ppm)
        self.resetSpacing()
        self.plotReset()
        self.showFid()
        
    def upd(self): #get new data from the data instance
        updateVar = self.data.getBlock(self.axes,self.axes2,self.locList,self.stackBegin, self.stackEnd, self.stackStep)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.sw = updateVar[2]
        self.spec = updateVar[3]
        self.wholeEcho = updateVar[4]
        self.xax=updateVar[5]
        self.xax2=updateVar[6]
        self.ref=updateVar[7]
        self.ref2=updateVar[8]

    def setBlock(self,axes,axes2,locList,stackBegin=None,stackEnd=None,stackStep=None): #change the slice 
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.locList = locList
        self.upd()
        self.resetSpacing()
        self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)
        
    def stackSelect(self,stackBegin, stackEnd, stackStep):
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.upd()
        self.showFid()

    def setPhaseInter(self, phase0in, phase1in): #interactive changing the phase without editing the actual data
        phase0=float(phase0in)
        phase1=float(phase1in)
        if self.spec==0:
            tmpdata=self.fourierLocal(self.data1D,0)
            tmpdata=tmpdata*np.exp(phase0*1j)
            tmpdata=tmpdata*np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata[0]),1.0/self.sw))*phase1*1j)],len(tmpdata),axis=0)
            tmpdata=self.fourierLocal(tmpdata,1)
        else:
            tmpdata=self.data1D*np.exp(phase0*1j)
            tmpdata=tmpdata*np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata[0]),1.0/self.sw))*phase1*1j)],len(tmpdata),axis=0)
        self.showFid(tmpdata)

    def apodPreview(self,lor=None,gauss=None, cos2=None, hamming=None,shift=0.0,shifting=0.0,shiftingAxes=None): #display the 1D data including the apodization function
        t=np.arange(0,len(self.data1D[0]))/(self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                print('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.data.shape[self.axes2])[slice(self.stackBegin,self.stackEnd,self.stackStep)]
                x=np.ones((len(ar),len(self.data1D[0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting*ar[i]
                    t2 = t - shift1
                    x2=np.ones(len(self.data1D[0]))
                    if lor is not None:
                        x2=x2*np.exp(-lor*abs(t2))
                    if gauss is not None:
                        x2=x2*np.exp(-(gauss*t2)**2)
                    if cos2 is not None:
                        x2=x2*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
                    if hamming is not None:
                        alpha = 0.53836 # constant for hamming window
                        x2=x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
                    if self.wholeEcho:
                        x2[2-1:-(len(x2)/2+1):-1]=x2[:len(x2)/2]
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting*self.locList[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting*self.locList[shiftingAxes-2]
                else:
                    shift += shifting*self.locList[shiftingAxes-1]
                t2 = t - shift
                x=np.ones(len(self.data1D[0]))
                if lor is not None:
                    x=x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x=x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
                if self.wholeEcho:
                    x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
                x = np.repeat([x],len(self.data1D),axis=0)
        else:
            t2 = t - shift
            x=np.ones(len(self.data1D[0]))
            if lor is not None:
                x=x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x=x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
            x = np.repeat([x],len(self.data1D),axis=0)
        y = self.data1D
        a=self.fig.gca()
        a.cla()
        if self.spec ==1:
            y=np.fft.ifftn(np.fft.ifftshift(y,axes=1),axes=[1])
            y= y*x
            y=np.fft.fftshift(np.fft.fftn(y,axes=[1]),axes=1)
        else:
            y= y*x
        if self.spec==0:
            if self.plotType==0:
                self.showFid(y,[t],x*np.amax(np.real(self.data1D)),['g'],old=True)
            elif self.plotType==1:
                self.showFid(y,[t],x*np.amax(np.imag(self.data1D)),['g'],old=True)
            elif self.plotType==2:
                self.showFid(y,[t],x*np.amax(np.amax(np.real(self.data1D)),np.amax(np.imag(self.data1D))),['g'],old=True)
            elif self.plotType==3:
                self.showFid(y,[t],x*np.amax(np.abs(self.data1D)),['g'],old=True)
        else:
            self.showFid(y)

    def setSizePreview(self,size): #set size only on local data
        if size > len(self.data1D[0]):
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=1)
                self.data1D = np.concatenate((np.pad(tmpdata[0],((0,0),(0,size-len(self.data1D[0]))),'constant',constant_values=0),tmpdata[1]),axis=1)
            else:
                self.data1D = np.pad(self.data1D,((0,0),(0,size-len(self.data1D[0]))),'constant',constant_values=0)
        else:
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=1)
                self.data1D = np.concatenate((tmpdata[:,:np.ceil(size/2.0)],tmpdata[:,size/2:]),axis=1)
            else:
                self.data1D = self.data1D[:,:size]
        if self.spec==0:
            self.xax=np.arange(len(self.data1D[0]))/self.sw
        elif self.spec==1:
            self.xax=np.fft.fftshift(np.fft.fftfreq(len(self.data1D[0]),1.0/self.sw)) 
        self.plotReset()
        self.showFid()
        self.upd()

    def setSwapEchoPreview(self,idx):
        self.data1D = np.concatenate((self.data1D[:,idx:],self.data1D[:,:idx]),axis=1)
        self.plotReset()
        self.showFid()
        self.upd()
        
    def setShiftPreview(self,shift):
        tmpData = np.roll(self.data1D,shift)
        if shift<0:
            tmpData[:,shift:] = tmpData[:,shift:]*0
        else:
            tmpData[:,:shift] = tmpData[:,:shift]*0
        self.showFid(tmpData)

    def dcOffset(self,pos1,pos2):
        minPos = int(min(pos1,pos2))
        maxPos = int(max(pos1,pos2))
        if minPos != maxPos:
            self.showFid(self.data1D-np.mean(self.data1D[:,minPos:maxPos]))
            
    def applydcOffset(self,pos1,pos2):
        minPos = int(min(pos1,pos2))
        maxPos = int(max(pos1,pos2))
        returnValue = self.data.dcOffset(-np.mean(self.data1D[:,minPos:maxPos]))
        self.upd()
        self.showFid()
        return returnValue

    def shearing(self,shear,axes,axes2):
        returnValue = self.data.shear(shear,axes,axes2)
        self.upd()
        self.showFid()
        return returnValue
    
    def setSpacing(self, spacing):
        self.spacing = spacing
        self.showFid()

    def resetSpacing(self):
        difference = np.diff(self.data1D,axis=0)
        if self.plotType==0:
            difference = np.amax(np.real(difference))
            amp = np.amax(np.real(self.data1D))-np.amin(np.real(self.data1D))
        elif self.plotType==1:
            difference = np.amax(np.imag(difference))
            amp = np.amax(np.imag(self.data1D))-np.amin(np.imag(self.data1D))
        elif self.plotType==2:
            difference = np.amax((np.real(difference),np.imag(difference)))
            amp = np.amax((np.real(self.data1D),np.imag(self.data1D)))-np.amin((np.real(self.data1D),np.imag(self.data1D)))
        elif self.plotType==3:
            difference = np.amax(np.abs(difference))
            amp = np.amax(np.abs(self.data1D))-np.amin(np.abs(self.data1D))
        self.spacing = difference + 0.1*amp

    def getDisplayedData(self):
        if self.plotType==0:
            return np.real(self.data1D[0])
        elif self.plotType==1:
            return np.imag(self.data1D[0])
        elif self.plotType==2:
            return np.real(self.data1D[0])
        elif self.plotType==3:
            return np.abs(self.data1D[0])      

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        a=self.fig.gca()
        a.cla()
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if old:
            if (self.plotType==0):
                for num in range(len(self.data1D)):
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==1):
                for num in range(len(self.data1D)):
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.imag(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==2):
                for num in range(len(self.data1D)):
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==3):
                for num in range(len(self.data1D)):
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.abs(self.data1D[num]),c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraY)):
                a.plot(extraX[0]*axMult+axAdd,num*self.spacing+extraY[num],c=extraColor[0])
        if (self.plotType==0):
            for num in range(len(tmpdata)):
                if num is 0:
                    self.line = a.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')
                else:
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')
        elif(self.plotType==1):
            for num in range(len(tmpdata)):
                if num is 0:
                    self.line = a.plot(self.xax*axMult+axAdd,num*self.spacing+np.imag(tmpdata[num]),c='b')
                else:
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.imag(tmpdata[num]),c='b')
        elif(self.plotType==2):
            for num in range(len(tmpdata)):
                a.plot(self.xax*axMult+axAdd,num*self.spacing+np.imag(tmpdata[num]),c='r')
                if num is 0:
                    self.line = a.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')
                else:
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')
        elif(self.plotType==3):
            for num in range(len(tmpdata)):
                if num is 0:
                    self.line = a.plot(self.xax*axMult+axAdd,num*self.spacing+np.abs(tmpdata[num]),c='b')
                else:
                    a.plot(self.xax*axMult+axAdd,num*self.spacing+np.abs(tmpdata[num]),c='b')
        a.set_title("TD"+str(self.axes+1))
        #a.set_xlabel('X axis label')
        if self.spec==0:
            if self.axType == 0:
                a.set_xlabel('Time [s]')
            elif self.axType == 1:
                a.set_xlabel('Time [ms]')
            elif self.axType == 2:
                a.set_xlabel(r'Time [$\mu$s]')
            else:
                a.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                a.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    a.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    a.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    a.set_xlabel('Frequency [MHz]')
                else:
                    a.set_xlabel('User defined')
        else:
            a.set_xlabel('')
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        a.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        a.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
        a=self.fig.gca()
        incr = np.repeat(np.arange(len(self.data1D)).reshape((len(self.data1D),1)),len(self.data1D[0]),axis=1)*self.spacing
        if self.plotType==0:
            miny = np.amin(np.real(self.data1D)+incr)
            maxy = np.amax(np.real(self.data1D)+incr)
        elif self.plotType==1:
            miny = np.amin(np.imag(self.data1D)+incr)
            maxy = np.amax(np.imag(self.data1D)+incr)
        elif self.plotType==2:
            miny = np.amin((np.amin(np.real(self.data1D)+incr),np.amin(np.imag(self.data1D)+incr)))
            maxy = np.amax((np.amax(np.real(self.data1D)+incr),np.amax(np.imag(self.data1D)+incr)))
        elif self.plotType==3:
            miny = np.amin(np.abs(self.data1D)+incr)
            maxy = np.amax(np.abs(self.data1D)+incr)
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=max(self.xax*axMult+axAdd)
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)

#########################################################################################################
#the class from which the arrayed data is displayed, the operations which only edit the content of this class are for previewing
class CurrentArrayed(Current1D):
    def __init__(self, root, data, axes=None, axes2=None, locList=None, plotType=0, axType=1, ppm=False, stackBegin=None, stackEnd=None, stackStep=None):
        self.data = data
        if axes2 is None:
            self.axes2 = len(data.data.shape)-2
        else:
            self.axes2 = axes2            #dimension from which the spectra are stacked
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if locList is None:
            self.resetLocList()
        self.spacing = 0
        Current1D.__init__(self, root, data, axes, locList, plotType, axType, ppm)
        self.resetSpacing()
        self.plotReset()
        self.showFid()
        
    def upd(self): #get new data from the data instance
        updateVar = self.data.getBlock(self.axes,self.axes2,self.locList,self.stackBegin, self.stackEnd, self.stackStep)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.sw = updateVar[2]
        self.spec = updateVar[3]
        self.wholeEcho = updateVar[4]
        self.xax=updateVar[5]
        self.xax2=updateVar[6]
        self.ref=updateVar[7]
        self.ref2=updateVar[8]
 
    def setBlock(self,axes,axes2,locList,stackBegin=None,stackEnd=None,stackStep=None): #change the slice 
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.locList = locList
        self.upd()
        self.resetSpacing()
        self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)
        
    def stackSelect(self,stackBegin, stackEnd, stackStep):
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.upd()
        self.showFid()

    def setPhaseInter(self, phase0in, phase1in): #interactive changing the phase without editing the actual data
        phase0=float(phase0in)
        phase1=float(phase1in)
        if self.spec==0:
            tmpdata=self.fourierLocal(self.data1D,0)
            tmpdata=tmpdata*np.exp(phase0*1j)
            tmpdata=tmpdata*np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata[0]),1.0/self.sw))*phase1*1j)],len(tmpdata),axis=0)
            tmpdata=self.fourierLocal(tmpdata,1)
        else:
            tmpdata=self.data1D*np.exp(phase0*1j)
            tmpdata=tmpdata*np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata[0]),1.0/self.sw))*phase1*1j)],len(tmpdata),axis=0)
        self.showFid(tmpdata)


    def apodPreview(self,lor=None,gauss=None, cos2=None, hamming=None,shift=0.0,shifting=0.0,shiftingAxes=None): #display the 1D data including the apodization function
        t=np.arange(0,len(self.data1D[0]))/(self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                print('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.data.shape[self.axes2])[slice(self.stackBegin,self.stackEnd,self.stackStep)]
                x=np.ones((len(ar),len(self.data1D[0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting*ar[i]
                    t2 = t - shift1
                    x2=np.ones(len(self.data1D[0]))
                    if lor is not None:
                        x2=x2*np.exp(-lor*abs(t2))
                    if gauss is not None:
                        x2=x2*np.exp(-(gauss*t2)**2)
                    if cos2 is not None:
                        x2=x2*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
                    if hamming is not None:
                        alpha = 0.53836 # constant for hamming window
                        x2=x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
                    if self.wholeEcho:
                        x2[2-1:-(len(x2)/2+1):-1]=x2[:len(x2)/2]
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting*self.locList[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting*self.locList[shiftingAxes-2]
                else:
                    shift += shifting*self.locList[shiftingAxes-1]
                t2 = t - shift
                x=np.ones(len(self.data1D[0]))
                if lor is not None:
                    x=x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x=x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
                if self.wholeEcho:
                    x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
                x = np.repeat([x],len(self.data1D),axis=0)
        else:
            t2 = t - shift
            x=np.ones(len(self.data1D[0]))
            if lor is not None:
                x=x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x=x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
            x = np.repeat([x],len(self.data1D),axis=0)
        y = self.data1D
        a=self.fig.gca()
        a.cla()
        if self.spec ==1:
            y=np.fft.ifftn(np.fft.ifftshift(y,axes=1),axes=[1])
            y= y*x
            y=np.fft.fftshift(np.fft.fftn(y,axes=[1]),axes=1)
        else:
            y= y*x
        if self.spec==0:
            if self.plotType==0:
                self.showFid(y,[t],x*np.amax(np.real(self.data1D)),['g'],old=True)
            elif self.plotType==1:
                self.showFid(y,[t],x*np.amax(np.imag(self.data1D)),['g'],old=True)
            elif self.plotType==2:
                self.showFid(y,[t],x*np.amax(np.amax(np.real(self.data1D)),np.amax(np.imag(self.data1D))),['g'],old=True)
            elif self.plotType==3:
                self.showFid(y,[t],x*np.amax(np.abs(self.data1D)),['g'],old=True)
        else:
            self.showFid(y)
                 
    def setSizePreview(self,size): #set size only on local data
        if size > len(self.data1D[0]):
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=1)
                self.data1D = np.concatenate((np.pad(tmpdata[0],((0,0),(0,size-len(self.data1D[0]))),'constant',constant_values=0),tmpdata[1]),axis=1)
            else:
                self.data1D = np.pad(self.data1D,((0,0),(0,size-len(self.data1D[0]))),'constant',constant_values=0)
        else:
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=1)
                self.data1D = np.concatenate((tmpdata[:,:np.ceil(size/2.0)],tmpdata[:,size/2:]),axis=1)
            else:
                self.data1D = self.data1D[:,:size]
        if self.spec==0:
            self.xax=np.arange(len(self.data1D[0]))/self.sw
        elif self.spec==1:
            self.xax=np.fft.fftshift(np.fft.fftfreq(len(self.data1D[0]),1.0/self.sw)) 
        self.plotReset()
        self.showFid()
        self.upd()

    def setSwapEchoPreview(self,idx):
        self.data1D = np.concatenate((self.data1D[:,idx:],self.data1D[:,:idx]),axis=1)
        self.plotReset()
        self.showFid()
        self.upd()
        
    def setShiftPreview(self,shift):
        tmpData = np.roll(self.data1D,shift)
        if shift<0:
            tmpData[:,shift:] = tmpData[:,shift:]*0
        else:
            tmpData[:,:shift] = tmpData[:,:shift]*0
        self.showFid(tmpData)

    def dcOffset(self,pos1,pos2):
        minPos = int(min(pos1,pos2))
        maxPos = int(max(pos1,pos2))
        if minPos != maxPos:
            self.showFid(self.data1D-np.mean(self.data1D[:,minPos:maxPos]))
    
    def ACMEentropy(self,phaseIn,phaseAll=True):
        phase0=phaseIn[0]
        if phaseAll:
            phase1=phaseIn[1]
        else:
            phase1=0.0
        L = len(self.data1D)
        if self.spec==1:
            x=np.fft.fftshift(np.fft.fftfreq(L,1.0/self.sw))
        if self.spec>0:
            s0 = self.data1D*np.exp(1j*(phase0+phase1*x))
        else:
            s0 = np.fft.fftshift(np.fft.fft(self.data1D))*np.exp(1j*(phase0+phase1*x))
        s2 = np.real(s0)
        ds1 = np.abs((s2[3:L]-s2[1:L-2])/2.0)
        p1 = ds1/sum(ds1)
        p1[np.where(p1 == 0)] = 1
        h1  = -p1*np.log(p1)
        H1  = sum(h1)
        Pfun = 0.0
        as1 = s2 - np.abs(s2)
        sumas   = sum(as1)
        if (np.real(sumas) < 0): 
            Pfun = Pfun + sum(as1**2)/4/L**2
        return H1+1000*Pfun 

    def autoPhase(self,phaseNum):
        if phaseNum == 0:
            phases = scipy.optimize.fmin(func=self.ACMEentropy,x0=[0],args=(False,))
        elif phaseNum == 1:
            phases = scipy.optimize.fmin(func=self.ACMEentropy,x0=[0,0])
        return phases

    def setSpacing(self, spacing):
        self.spacing = spacing
        self.showFid()

    def resetSpacing(self):
        self.spacing = (self.xax[-1]-self.xax[0])*1.1

    def getDisplayedData(self):
        if self.plotType==0:
            return np.real(self.data1D[0])
        elif self.plotType==1:
            return np.imag(self.data1D[0])
        elif self.plotType==2:
            return np.real(self.data1D[0])
        elif self.plotType==3:
            return np.abs(self.data1D[0])      

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        a=self.fig.gca()
        a.cla()
        if self.spec > 0:
            direc = slice(None,None,-1)
        else:
            direc = slice(None,None,1)
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if old:
            if (self.plotType==0):
                for num in range(len(self.data1D)):
                    a.plot((num*self.spacing+self.xax)*axMult+axAdd,np.real(self.data1D[num])[direc],c='k',alpha=0.2)
            elif(self.plotType==1):
                for num in range(len(self.data1D)):
                    a.plot((num*self.spacing+self.xax)*axMult+axAdd,np.imag(self.data1D[num])[direc],c='k',alpha=0.2)
            elif(self.plotType==2):
                for num in range(len(self.data1D)):
                    a.plot((num*self.spacing+self.xax)*axMult+axAdd,np.real(self.data1D[num])[direc],c='k',alpha=0.2)
            elif(self.plotType==3):
                for num in range(len(self.data1D)):
                    a.plot((num*self.spacing+self.xax)*axMult+axAdd,np.abs(self.data1D[num])[direc],c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraY)):
                a.plot((num*self.spacing+extraX[0])*axMult+axAdd,extraY[num][direc],c=extraColor[0])

        self.line = []
        if (self.plotType==0):
            for num in range(len(tmpdata)):
                self.line.append(a.plot((num*self.spacing+self.xax)*axMult,np.real(tmpdata[num])[direc],c='b')[0])
        elif(self.plotType==1):
            for num in range(len(tmpdata)):
                self.line.append(a.plot((num*self.spacing+self.xax)*axMult,np.imag(tmpdata[num])[direc],c='b')[0])
        elif(self.plotType==2):
            for num in range(len(tmpdata)):
                a.plot((num*self.spacing+self.xax)*axMult,np.imag(tmpdata[num])[direc],c='r')
                self.line.append(a.plot((num*self.spacing+self.xax)*axMult,np.real(tmpdata[num])[direc],c='b')[0])
        elif(self.plotType==3):
            for num in range(len(tmpdata)):
                self.line.append(a.plot((num*self.spacing+self.xax)*axMult,np.abs(tmpdata[num])[direc],c='b')[0])
        a.set_title("TD"+str(self.axes+1))
        if self.spec==0:
            if self.axType == 0:
                a.set_xlabel('Time [s]')
            elif self.axType == 1:
                a.set_xlabel('Time [ms]')
            elif self.axType == 2:
                a.set_xlabel(r'Time [$\mu$s]')
            else:
                a.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                a.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    a.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    a.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    a.set_xlabel('Frequency [MHz]')
                else:
                    a.set_xlabel('User defined')
        else:
            a.set_xlabel('')
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        a.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        a.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
        a=self.fig.gca()
        if self.plotType==0:
            miny = np.amin(np.real(self.data1D))
            maxy = np.amax(np.real(self.data1D))
        elif self.plotType==1:
            miny = np.amin(np.imag(self.data1D))
            maxy = np.amax(np.imag(self.data1D))
        elif self.plotType==2:
            miny = np.amin((np.amin(np.real(self.data1D)),np.amin(np.imag(self.data1D))))
            maxy = np.amax((np.amax(np.real(self.data1D)),np.amax(np.imag(self.data1D))))
        elif self.plotType==3:
            miny = np.amin(np.abs(self.data1D))
            maxy = np.amax(np.abs(self.data1D))
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=(max(self.xax)+(len(self.data1D)-1)*self.spacing)*axMult+axAdd
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)

#########################################################################################################
#the class from which the contour data is displayed, the operations which only edit the content of this class are for previewing
class CurrentContour(Current1D):
    def __init__(self, root, data, axes=None, axes2=None, locList=None, plotType=0, axType=1, ppm=False):
        self.data = data
        if axes2 is None:
            self.axes2 = len(data.data.shape)-2
        else:
            self.axes2 = axes2            
        if locList is None:
            self.resetLocList()
        self.spacing = 0
        Current1D.__init__(self, root, data, axes, locList, plotType, axType, ppm)
        self.plotReset()
        self.showFid()
        
    def upd(self): #get new data from the data instance
        updateVar = self.data.getBlock(self.axes,self.axes2,self.locList)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.sw = updateVar[2]
        self.spec = updateVar[3]
        self.wholeEcho = updateVar[4]
        self.xax=updateVar[5]
        self.xax2=updateVar[6]
        self.ref=updateVar[7]
        self.ref2=updateVar[8]

    def setBlock(self,axes,axes2,locList): #change the slice 
        self.axes = axes
        self.axes2 = axes2
        self.locList = locList
        self.upd()
        self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)

    def setPhaseInter(self, phase0in, phase1in): #interactive changing the phase without editing the actual data
        phase0=float(phase0in)
        phase1=float(phase1in)
        if self.spec==0:
            tmpdata=self.fourierLocal(self.data1D,0)
            tmpdata=tmpdata*np.exp(phase0*1j)
            tmpdata=tmpdata*np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata[0]),1.0/self.sw))*phase1*1j)],len(tmpdata),axis=0)
            tmpdata=self.fourierLocal(tmpdata,1)
        else:
            tmpdata=self.data1D*np.exp(phase0*1j)
            tmpdata=tmpdata*np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata[0]),1.0/self.sw))*phase1*1j)],len(tmpdata),axis=0)
        self.showFid(tmpdata)


    def apodPreview(self,lor=None,gauss=None, cos2=None, hamming=None,shift=0.0,shifting=0.0,shiftingAxes=None): #display the 1D data including the apodization function
        t=np.arange(0,len(self.data1D[0]))/(self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                print('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.data.shape[self.axes2])
                x=np.ones((len(ar),len(self.data1D[0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting*ar[i]
                    t2 = t - shift1
                    x2=np.ones(len(self.data1D[0]))
                    if lor is not None:
                        x2=x2*np.exp(-lor*abs(t2))
                    if gauss is not None:
                        x2=x2*np.exp(-(gauss*t2)**2)
                    if cos2 is not None:
                        x2=x2*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
                    if hamming is not None:
                        alpha = 0.53836 # constant for hamming window
                        x2=x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
                    if self.wholeEcho:
                        x2[2-1:-(len(x2)/2+1):-1]=x2[:len(x2)/2]
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting*self.locList[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting*self.locList[shiftingAxes-2]
                else:
                    shift += shifting*self.locList[shiftingAxes-1]
                t2 = t - shift
                x=np.ones(len(self.data1D[0]))
                if lor is not None:
                    x=x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x=x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
                if self.wholeEcho:
                    x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
                x = np.repeat([x],len(self.data1D),axis=0)
        else:
            t2 = t - shift
            x=np.ones(len(self.data1D[0]))
            if lor is not None:
                x=x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x=x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,0.5*np.pi,len(self.data1D[0]))))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0,np.pi,len(self.data1D)))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
            x = np.repeat([x],len(self.data1D),axis=0)
        y = self.data1D
        a=self.fig.gca()
        a.cla()
        if self.spec ==1:
            y=np.fft.ifftn(np.fft.ifftshift(y,axes=1),axes=[1])
            y= y*x
            y=np.fft.fftshift(np.fft.fftn(y,axes=[1]),axes=1)
        else:
            y= y*x
        # if self.spec==0:
        #     if self.plotType==0:
        #         self.showFid(y,[t],x*np.amax(np.real(self.data1D)),['g'],old=True)
        #     elif self.plotType==1:
        #         self.showFid(y,[t],x*np.amax(np.imag(self.data1D)),['g'],old=True)
        #     elif self.plotType==2:
        #         self.showFid(y,[t],x*np.amax(np.amax(np.real(self.data1D)),np.amax(np.imag(self.data1D))),['g'],old=True)
        #     elif self.plotType==3:
        #         self.showFid(y,[t],x*np.amax(np.abs(self.data1D)),['g'],old=True)
        # else:
        #     self.showFid(y)
        self.showFid(y)
                 
    def setSizePreview(self,size): #set size only on local data
        if size > len(self.data1D[0]):
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=1)
                self.data1D = np.concatenate((np.pad(tmpdata[0],((0,0),(0,size-len(self.data1D[0]))),'constant',constant_values=0),tmpdata[1]),axis=1)
            else:
                self.data1D = np.pad(self.data1D,((0,0),(0,size-len(self.data1D[0]))),'constant',constant_values=0)
        else:
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=1)
                self.data1D = np.concatenate((tmpdata[:,:np.ceil(size/2.0)],tmpdata[:,size/2:]),axis=1)
            else:
                self.data1D = self.data1D[:,:size]
        if self.spec==0:
            self.xax=np.arange(len(self.data1D[0]))/self.sw
        elif self.spec==1:
            self.xax=np.fft.fftshift(np.fft.fftfreq(len(self.data1D[0]),1.0/self.sw)) 
        self.plotReset()
        self.showFid()
        self.upd()

    def setSwapEchoPreview(self,idx):
        self.data1D = np.concatenate((self.data1D[:,idx:],self.data1D[:,:idx]),axis=1)
        self.plotReset()
        self.showFid()
        self.upd()
        
    def setShiftPreview(self,shift):
        tmpData = np.roll(self.data1D,shift)
        if shift<0:
            tmpData[:,shift:] = tmpData[:,shift:]*0
        else:
            tmpData[:,:shift] = tmpData[:,:shift]*0
        self.showFid(tmpData)

    def dcOffset(self,pos1,pos2):
        minPos = int(min(pos1,pos2))
        maxPos = int(max(pos1,pos2))
        if minPos != maxPos:
            self.showFid(self.data1D-np.mean(self.data1D[:,minPos:maxPos]))
    
    def ACMEentropy(self,phaseIn,phaseAll=True):
        phase0=phaseIn[0]
        if phaseAll:
            phase1=phaseIn[1]
        else:
            phase1=0.0
        L = len(self.data1D)
        if self.spec==1:
            x=np.fft.fftshift(np.fft.fftfreq(L,1.0/self.sw))
        if self.spec>0:
            s0 = self.data1D*np.exp(1j*(phase0+phase1*x))
        else:
            s0 = np.fft.fftshift(np.fft.fft(self.data1D))*np.exp(1j*(phase0+phase1*x))
        s2 = np.real(s0)
        ds1 = np.abs((s2[3:L]-s2[1:L-2])/2.0)
        p1 = ds1/sum(ds1)
        p1[np.where(p1 == 0)] = 1
        h1  = -p1*np.log(p1)
        H1  = sum(h1)
        Pfun = 0.0
        as1 = s2 - np.abs(s2)
        sumas   = sum(as1)
        if (np.real(sumas) < 0): 
            Pfun = Pfun + sum(as1**2)/4/L**2
        return H1+1000*Pfun 

    def autoPhase(self,phaseNum):
        if phaseNum == 0:
            phases = scipy.optimize.fmin(func=self.ACMEentropy,x0=[0],args=(False,))
        elif phaseNum == 1:
            phases = scipy.optimize.fmin(func=self.ACMEentropy,x0=[0,0])
        return phases

    def getDisplayedData(self):
        if self.plotType==0:
            return np.real(self.data1D[0])
        elif self.plotType==1:
            return np.imag(self.data1D[0])
        elif self.plotType==2:
            return np.real(self.data1D[0])
        elif self.plotType==3:
            return np.abs(self.data1D[0])      

    def showFid(self, tmpdata=None): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        a=self.fig.gca()
        a.cla()
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        X, Y = np.meshgrid(self.xax*axMult+axAdd,self.xax2)
        self.line = []
        if (self.plotType==0):
            self.line.append(a.contour(X, Y, np.real(tmpdata),c='b'))
        elif(self.plotType==1):
            self.line.append(a.contour(X, Y, np.imag(tmpdata),c='b'))
        elif(self.plotType==2):
            print('type not supported')
            self.line.append(a.contour(X, Y, np.real(tmpdata),c='b'))
        elif(self.plotType==3):
            self.line.append(a.contour(X, Y, np.abs(tmpdata),c='b'))
        a.set_title("TD"+str(self.axes+1))
        if self.spec==0:
            if self.axType == 0:
                a.set_xlabel('Time [s]')
            elif self.axType == 1:
                a.set_xlabel('Time [ms]')
            elif self.axType == 2:
                a.set_xlabel(r'Time [$\mu$s]')
            else:
                a.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                a.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    a.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    a.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    a.set_xlabel('Frequency [MHz]')
                else:
                    a.set_xlabel('User defined')
        else:
            a.set_xlabel('')
        if self.spec:
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        a.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        a.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
        a=self.fig.gca()
        self.yminlim=min(self.xax2)
        self.ymaxlim=max(self.xax2)
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=max(self.xax*axMult+axAdd)
        if self.spec:
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)

    #The peakpicking function needs to be changed for contour plots
    def buttonRelease(self,event):
        a=self.fig.gca()
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0]=None
                    self.peakPick = False
                    axAdd = 0
                    if self.spec == 1:
                        if self.ppm:
                            axAdd = (self.freq-self.ref)/self.ref*1e6
                            axMult = 1e6/self.ref
                        else:
                            axMult = 1.0/(1000.0**self.axType)
                    elif self.spec == 0:
                        axMult = 1000.0**self.axType
                    xdata = self.xax*axMult+axAdd
                    ydata = self.xax2
                    idx = np.argmin(np.abs(xdata-event.xdata))
                    if self.peakPickFunc is not None:
                        #self.peakPickFunc((idx,xdata[idx],ydata[idx]))
                        self.peakPickFunc((idx,xdata[idx],0))
                    if not self.peakPick: #check if peakpicking is still required
                        self.peakPickFunc = None
            else:
                self.leftMouse = False
                if self.rect[0] is not None:
                    self.rect[0].remove()
                if self.rect[1] is not None:
                    self.rect[1].remove()
                if self.rect[2] is not None:
                    self.rect[2].remove()
                if self.rect[3] is not None:
                    self.rect[3].remove()
                self.rect=[None,None,None,None]
                if self.zoomX2 is not None and self.zoomY2 is not None:
                    self.xminlim=min([self.zoomX1,self.zoomX2])
                    self.xmaxlim=max([self.zoomX1,self.zoomX2])
                    self.yminlim=min([self.zoomY1,self.zoomY2])
                    self.ymaxlim=max([self.zoomY1,self.zoomY2])
                    if self.spec > 0:
                        a.set_xlim(self.xmaxlim,self.xminlim)
                    else:
                        a.set_xlim(self.xminlim,self.xmaxlim)
                    a.set_ylim(self.yminlim,self.ymaxlim)
                self.zoomX1=None
                self.zoomX2=None #WF: should also be cleared, memory of old zoom
                self.zoomY1=None
                self.zoomY2=None #WF: should also be cleared, memory of old zoom
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()
