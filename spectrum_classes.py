import matplotlib
import numpy as np
import scipy.optimize
from scipy.interpolate import UnivariateSpline
import scipy.signal
import scipy.ndimage
import copy
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from spectrumFrame import Plot1DFrame
import sys

#########################################################################
#the generic data class
class Spectrum(object):
    def __init__(self, data, loadFunc ,freq , sw , spec=None, wholeEcho=None, ref=None, xaxArray=None):
        self.dim = len(data.shape)                    #number of dimensions
        self.data = np.array(data,dtype=complex)      #data of dimension dim
        self.loadFunc = loadFunc
        self.freq = np.array(freq)                              #array of center frequency (length is dim, MHz)
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
            self.ref = self.dim*[None]
        else:
            self.ref = ref
        if xaxArray is None:
            self.xaxArray = [[] for i in range(self.dim)]
            self.resetXax()
        else:
            self.xaxArray = xaxArray

    def reload(self,mainProgram):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.reload(mainProgram))
        self.restoreData(self.loadFunc(mainProgram),None)
        return returnValue
            
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

    def setXax(self,xax,axes):
        oldXax = self.xaxArray[axes]
        self.xaxArray[axes]=xax
        return lambda self: self.setXax(oldXax,axes)
                
    def insert(self,data,pos,axes):
        self.data = np.insert(self.data,[pos],data,axis=axes)
        self.resetXax(axes)
        return lambda self: self.remove(range(pos,pos+data.shape[axes]),axes)

    def remove(self,pos,axes):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.remove(pos,axes))
        tmpdata = np.delete(self.data,pos,axes)
        if (np.array(tmpdata.shape) != 0).all():
            self.data = tmpdata
            self.xaxArray[axes] = np.delete(self.xaxArray[axes],pos)
            return returnValue
        else:
            print('Cannot delete all data')
            return None

    def add(self,data):
        self.data = self.data + data
        return lambda self: self.subtract(data)
        
    def subtract(self,data):
        self.data = self.data - data
        return lambda self: self.add(data)

    def concatenate(self,axes):
        splitVal = self.data.shape[axes]
        self.data = np.concatenate(self.data,axis=axes)
        self.dim = len(self.data.shape)
        self.freq = np.delete(self.freq, axes)
        self.sw = np.delete(self.sw, axes)
        self.spec = np.delete(self.spec, axes)
        self.wholeEcho =  np.delete(self.wholeEcho, axes)
        self.ref = np.delete(self.ref, axes)
        del self.xaxArray[0]
        #self.resetXax(axes)
        self.resetXax()
        return lambda self: self.split(splitVal,axes)
    
    def split(self,sections,axes):
        self.data = np.array(np.split(self.data,sections,axis=axes))
        self.dim = len(self.data.shape)
        self.freq = np.insert(self.freq, 0, self.freq[axes])
        self.sw = np.insert(self.sw, 0 , self.sw[axes])
        self.spec = np.insert(self.spec, 0 , self.spec[axes])
        self.wholeEcho =  np.insert(self.wholeEcho, 0 , self.wholeEcho[axes])
        self.ref = np.insert(self.ref, 0 , self.ref[axes])
        self.xaxArray.insert(0,[])
        self.resetXax(0)
        self.resetXax(axes+1)
        return lambda self: self.concatenate(axes)
    
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
    
    def states(self,axes):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.states(axes))
        if self.data.shape[axes]%2 != 0:
            print("data has to be even for States")
            return None
        tmpdata = np.real(self.data)
        slicing1 = (slice(None),) * axes + (slice(None,None,2),) + (slice(None),)*(self.dim-1-axes)
        slicing2 = (slice(None),) * axes + (slice(1,None,2),) + (slice(None),)*(self.dim-1-axes)
        tmpdata = tmpdata[slicing1]+1j*tmpdata[slicing2]
        self.data = tmpdata
        self.resetXax(axes)
        return returnValue
    
    def statesTPPI(self,axes):
        copyData=copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.statesTPPI(axes))
        if self.data.shape[axes]%2 != 0:
            print("data has to be even for States-TPPI")
            return None
        tmpdata = np.real(self.data)
        slicing1 = (slice(None),) * axes + (slice(None,None,2),) + (slice(None),)*(self.dim-1-axes)
        slicing2 = (slice(None),) * axes + (slice(1,None,2),) + (slice(None),)*(self.dim-1-axes)
        tmpdata = tmpdata[slicing1]+1j*tmpdata[slicing2]
        tmpdata[slicing2] = -1*tmpdata[slicing2]
        self.data = tmpdata
        self.resetXax(axes)
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
        del self.xaxArray[axes]
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
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.data.shape[axes],1.0/self.sw[axes]))/self.sw[axes]*phase1*1j)
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
        if shiftingAxes==None:
            shiftingAxes = 0
            shifting = 0
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
                if self.wholeEcho[axes]:
                    x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
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
            if self.wholeEcho[axes]:
                x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
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
    
    def setRef(self,ref,axes):
        oldRef = self.ref[axes]
        self.ref[axes]=ref
        return lambda self: self.setRef(oldRef,axes)

    def setWholeEcho(self,val,axes):
        self.wholeEcho[axes]=val
        return lambda self: self.setWholeEcho(not val,axes)
    
    def setSize(self,size,axes):
        if axes>self.dim:
            print("axes bigger than dim in setsize")
            return None
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
        if self.spec[axes] > 0:
            self.fourier(axes,tmp=True)
            self.data = np.roll(self.data,shift,axes)
            if shift < 0:
                slicing = (slice(None),) * axes + (slice(shift,None),) + (slice(None),)*(self.dim-1-axes)
            else:
                slicing = (slice(None),) * axes + (slice(None,shift),) + (slice(None),)*(self.dim-1-axes)
            self.data[slicing]=self.data[slicing]*0
            self.fourier(axes,tmp=True)
        else:
            self.data = np.roll(self.data,shift,axes)
            if shift < 0:
                slicing = (slice(None),) * axes + (slice(shift,None),) + (slice(None),)*(self.dim-1-axes)
            else:
                slicing = (slice(None),) * axes + (slice(None,shift),) + (slice(None),)*(self.dim-1-axes)
            self.data[slicing]=self.data[slicing]*0
        return returnValue
    
    def LPSVD(self,axes):
        M = 1 #Number of frequencies
        y = self.data[axes][0:100]
        N=y.shape[0]						# # of complex data points in FID
        L=math.floor(N*3/4)						# linear prediction order L = 3/4*N
        A=scipy.linalg.hankel(np.conj(y[1:N-L+1]),np.conj(y[N-L:N]))	# backward prediction data matrix
        h=np.conj(y[0:N-L])					# backward prediction data vector
        U,S,V = np.linalg.svd(A)					# singular value decomposition
        S = np.diag(S)
        bias=np.mean(S[M:np.min([N-L-1,L])])	# bias compensation
        PolyCoef=np.dot(-V[:,0:M],np.dot(np.diag(1/(S[0:M]-bias)),np.dot(np.transpose(U[:,0:M]),h)))	# prediction polynomial coefficients
        s=np.conj(np.log(np.roots(np.append(PolyCoef[::-1],1))))		# polynomial rooting
        s = s[np.where(np.real(s)<0)[0]]
        a=1
        
    def dcOffset(self,offset):
        self.data = self.data-offset
        return lambda self: self.dcOffset(-offset)

    def fourier(self, axes,tmp=False):
        if axes>self.dim:
            print("axes bigger than dim in fourier")
            return None
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
            return None
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
        return copy.deepcopy((self.data[tuple(locList[:axes])+(slice(None),)+tuple(locList[axes:])],self.freq[axes],self.sw[axes],self.spec[axes],self.wholeEcho[axes],self.xaxArray[axes],self.ref[axes]))

    def getBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None):
        stackSlice = slice(stackBegin, stackEnd, stackStep)
        if axes == axes2:
            print("First and second axes are the same")
            return None
        elif axes < axes2:
            return copy.deepcopy((np.transpose(self.data[tuple(locList[:axes])+(slice(None),)+tuple(locList[axes:axes2-1])+(stackSlice,)+tuple(locList[axes2-1:])]),self.freq[axes],self.freq[axes2],self.sw[axes],self.sw[axes2],self.spec[axes],self.spec[axes2],self.wholeEcho[axes],self.wholeEcho[axes2],self.xaxArray[axes],self.xaxArray[axes2][stackSlice],self.ref[axes],self.ref[axes2]))
        elif axes > axes2:
            return copy.deepcopy((self.data[tuple(locList[:axes2])+(stackSlice,)+tuple(locList[axes2:axes-1])+(slice(None),)+tuple(locList[axes-1:])],self.freq[axes],self.freq[axes2],self.sw[axes],self.sw[axes2],self.spec[axes],self.spec[axes2],self.wholeEcho[axes],self.wholeEcho[axes2],self.xaxArray[axes],self.xaxArray[axes2][stackSlice],self.ref[axes],self.ref[axes2]))
        
    def restoreData(self,copyData,returnValue): # restore data from an old copy for undo purposes
        self.data = copyData.data
        self.dim = len(self.data.shape)                    #number of dimensions
        self.freq = copyData.freq                              #array of center frequency (length is dim, MHz)
        self.loadFunc = copyData.loadFunc
        self.sw = copyData.sw                                  #array of sweepwidths
        self.spec = copyData.spec
        self.wholeEcho = copyData.wholeEcho
        self.xaxArray = copyData.xaxArray
        self.ref = copyData.ref
        return returnValue


#########################################################################################################
#the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing
class Current1D(Plot1DFrame):
    def __init__(self, root, data,duplicateCurrent=None):
        Plot1DFrame.__init__(self,root)
        self.xax = None               #x-axis
        self.data = data              #the actual spectrum instance
        self.freq = None              #frequency of the slice 
        self.sw = None                #x-data display
        self.data1D = None            #the data1D
        self.spec = None              #boolean where False=time domain and True=spectral domain
        self.wholeEcho = None
        self.ref = None               #reference frequency
        if duplicateCurrent is None:
            self.ppm = False                #display frequency as ppm
            self.axes = len(self.data.data.shape)-1
            self.resetLocList() 
            self.plotType = 0
            self.axType = 1
        else:
            self.ppm = duplicateCurrent.ppm
            self.axes = duplicateCurrent.axes
            if isinstance(self,(CurrentStacked,CurrentArrayed,CurrentContour,CurrentSkewed)):
                if (len(duplicateCurrent.locList) == self.data.dim-2):
                    self.locList = duplicateCurrent.locList
                else:
                    if self.axes < self.axes2:
                        self.locList = np.delete(duplicateCurrent.locList,self.axes2-1)
                    else:
                        self.locList = np.delete(duplicateCurrent.locList,self.axes2)
            else:
                if (len(duplicateCurrent.locList) == self.data.dim-1):
                    self.locList = duplicateCurrent.locList
                else:
                    if self.axes < duplicateCurrent.axes2:
                        self.locList = np.insert(duplicateCurrent.locList,duplicateCurrent.axes2-1,0)
                    else:
                        self.locList = np.insert(duplicateCurrent.locList,duplicateCurrent.axes2,0)
            self.plotType = duplicateCurrent.plotType
            self.axType = duplicateCurrent.axType
        self.upd()   #get the first slice of data
        self.plotReset() #reset the axes limits
        self.showFid() #plot the data
        
    def copyCurrent(self,root, data):
        return Current1D(root,data,self)
        
    def upd(self): #get new data from the data instance
        if self.data.dim <= self.axes:
            self.axes = len(self.data.data.shape)-1
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
        if self.ref == None:
            self.ref = self.freq
        self.single = self.data1D.shape[-1]==1

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
        self.upd()
        if self.spec==0:
            tmpdata=self.fourierLocal(self.data1D,0)
        else:
            tmpdata = self.data1D
        tmpdata=tmpdata*np.exp(phase0*1j)
        if len(self.data1D.shape) > 1:
            mult = np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata),1.0/self.sw))/self.sw*phase1*1j)],len(tmpdata),axis=0)
        else:
            mult = np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata),1.0/self.sw))/self.sw*phase1*1j)
        tmpdata=tmpdata*mult
        if self.spec==0:
            tmpdata=self.fourierLocal(tmpdata,1)
        self.data1D = tmpdata
        self.showFid()

    def applyPhase(self, phase0, phase1):# apply the phase to the actual data
        phase0=float(phase0)
        phase1=float(phase1)
        returnValue = self.data.setPhase(phase0,phase1,self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['phase',(phase0,phase1,self.axes)])
        return returnValue

    def fourier(self): #fourier the actual data and replot
        returnValue = self.data.fourier(self.axes)
        self.upd()
        if isinstance(self,(CurrentStacked,CurrentArrayed)):
            self.resetSpacing()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['fourier',(self.axes,)])
        return returnValue

    def fftshift(self,inv=False): #fftshift the actual data and replot
        returnValue = self.data.fftshift(self.axes,inv)
        self.upd()
        self.showFid()
        self.root.addMacro(['fftshift',(self.axes,inv)])
        return returnValue
    
    def fourierLocal(self, fourData, spec): #fourier the local data for other functions
        ax = len(fourData.shape)-1
        a = fourData
        if spec==0:
            if not self.wholeEcho:
                slicing = (slice(None),) * ax + (0,)
                fourData[slicing]= fourData[slicing]*0.5
            fourData=np.fft.fftshift(np.fft.fftn(fourData,axes=[ax]),axes=ax)
        else:
            fourData=np.fft.ifftn(np.fft.ifftshift(fourData,axes=ax),axes=[ax])
            if not self.wholeEcho:
                slicing = (slice(None),) * ax + (0,) 
                fourData[slicing]= fourData[slicing]*2.0
        return fourData

    def apodPreview(self,lor=None,gauss=None, cos2=None, hamming=None ,shift=0.0,shifting=0.0,shiftingAxes=None): #display the 1D data including the apodization function
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                print('shiftingAxes cannot be equal to axes')
                return
            elif shiftingAxes < self.axes:
                shift += shifting*self.locList[shiftingAxes]
            else:
                shift += shifting*self.locList[shiftingAxes-1]
        length = len(self.data1D)
        t=np.arange(0,length)/(self.sw)
        t2=t-shift
        x=np.ones(length)
        if lor is not None:
            x=x*np.exp(-lor*abs(t2))
        if gauss is not None:
            x=x*np.exp(-(gauss*t2)**2)
        if cos2 is not None:
            x=x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/length+np.linspace(0,0.5*np.pi,len(self.data1D))))**2)
        if hamming is not None:
            alpha = 0.53836 # constant for hamming window
            x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/length+np.linspace(0,np.pi,length))))
        if self.wholeEcho:
            x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
        self.ax.cla()
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
        self.root.addMacro([u'apodize',(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes,self.axes)])
        return returnValue

    def setFreq(self,freq,sw): #set the frequency of the actual data
        returnValue = self.data.setFreq(freq,sw,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['freq',(freq,sw,self.axes)])
        return returnValue

    def setRef(self,ref): #set the frequency of the actual data
        returnValue = self.data.setRef(ref,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['ref',(ref,self.axes)])
        return returnValue

    def SN(self,minNoise,maxNoise,minPeak,maxPeak):
        minN = min(minNoise,maxNoise)
        maxN = max(minNoise,maxNoise)
        minP = min(minPeak,maxPeak)
        maxP = max(minPeak,maxPeak)
        if len(self.data1D.shape) > 1:
            tmpData = self.data1D[0]
        else:
            tmpData = self.data1D
        if (self.plotType==0):
            tmpData = np.real(tmpData)
        elif(self.plotType==1):
            tmpData = np.imag(tmpData)
        elif(self.plotType==2):
            tmpData = np.real(tmpData)
        elif(self.plotType==3):
            tmpData = np.abs(tmpData)
        return (np.amax(tmpData[minP:maxP])/(np.amax(tmpData[minN:maxN])-np.amin(tmpData[minN:maxN])))
    
    def fwhm(self,minPeak,maxPeak):
        minP = min(minPeak,maxPeak)
        maxP = max(minPeak,maxPeak)
        if len(self.data1D.shape) > 1:
            tmpData = self.data1D[0]
        else:
            tmpData = self.data1D
        if (self.plotType==0):
            tmpData = np.real(tmpData)
        elif(self.plotType==1):
            tmpData = np.imag(tmpData)
        elif(self.plotType==2):
            tmpData = np.real(tmpData)
        elif(self.plotType==3):
            tmpData = np.abs(tmpData)
        maxPos = np.argmax(tmpData[minP:maxP])
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        x= self.xax*axMult+axAdd
        maxX = x[minP:maxP][maxPos]
        spline = UnivariateSpline(x, tmpData-tmpData[minP:maxP][maxPos]/2.0, s=0)
        zeroPos = spline.roots()
        left = zeroPos[zeroPos>maxX]
        right = zeroPos[zeroPos<maxX]
        if right.size>0 and left.size>0:
            return abs(left[0]-right[-1])
        else:
            return 0.0
    
    def setSizePreview(self,size): #set size only on local data
        if len(self.data1D.shape) > 1:
            length = len(self.data1D[0])
        else:
            length = len(self.data1D)
        if size > length:
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=(len(self.data1D.shape)-1))
                if len(self.data1D.shape) > 1:
                    self.data1D = np.concatenate((np.pad(tmpdata[0],((0,0),(0,size-length)),'constant',constant_values=0),tmpdata[1]),axis=1)
                else:
                    self.data1D = np.concatenate((np.pad(tmpdata[0],(0,size-length),'constant',constant_values=0),tmpdata[1]))
            else:
                if len(self.data1D.shape) > 1:
                    self.data1D = np.pad(self.data1D,((0,0),(0,size-length)),'constant',constant_values=0)
                else:
                    self.data1D = np.pad(self.data1D,(0,size-length),'constant',constant_values=0)
        else:
            if self.wholeEcho:
                tmpdata = np.array_split(self.data1D,2,axis=(len(self.data1D.shape)-1))
                if len(self.data1D.shape) > 1:
                    self.data1D = np.concatenate((tmpdata[:,:np.ceil(size/2.0)],tmpdata[:,size/2:]),axis=1)
                else:
                    self.data1D = np.concatenate((tmpdata[:np.ceil(size/2.0)],tmpdata[size/2:]))
            else:
                if len(self.data1D.shape) > 1:
                    self.data1D = self.data1D[:,:size]
                else:
                    self.data1D = self.data1D[:size]
        if len(self.data1D.shape) > 1:
            length = len(self.data1D[0])
        else:
            length = len(self.data1D)
        if self.spec==0:
            self.xax=np.arange(length)/self.sw
        elif self.spec==1:
            self.xax=np.fft.fftshift(np.fft.fftfreq(length,1.0/self.sw))
        self.plotReset()
        self.showFid()
        self.upd()

    def applySize(self,size): #set size to the actual data
        returnValue = self.data.setSize(size,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['size',(size,self.axes)])
        return returnValue

    def changeSpec(self,val): #change from time to freq domain of the actual data
        returnValue = self.data.changeSpec(val,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['spec',(val,self.axes)])
        return returnValue
    
    def applySwapEcho(self,idx):
        returnValue = self.data.swapEcho(idx,self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['swapecho',(idx,self.axes)])
        return returnValue

    def setSwapEchoPreview(self,idx):
        if len(self.data1D.shape) > 1:
            self.data1D = np.concatenate((self.data1D[:,idx:],self.data1D[:,:idx]),axis=1)
        else:
            self.data1D = np.concatenate((self.data1D[idx:],self.data1D[:idx]))
        self.plotReset()
        self.showFid()
        self.upd()

    def setWholeEcho(self, value):
        if value == 0:
            returnValue = self.data.setWholeEcho(False,self.axes)
            self.wholeEcho = False
            self.root.addMacro(['wholeEcho',(False,self.axes)])
        else:
            returnValue = self.data.setWholeEcho(True,self.axes)
            self.wholeEcho = True
            self.root.addMacro(['wholeEcho',(True,self.axes)])
        return returnValue

    def applyShift(self,shift):
        returnValue = self.data.shiftData(shift,self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['shift',(shift,self.axes)])
        return returnValue

    def setShiftPreview(self,shift):
        tmpData = self.data1D
        dim = len(self.data1D.shape)
        if self.spec > 0:
            tmpData=self.fourierLocal(tmpData,1)
        tmpData = np.roll(tmpData,shift)
        if shift<0:
            tmpData[(slice(None),)*(dim-1)+(slice(shift,None),)] = tmpData[(slice(None),)*(dim-1)+(slice(shift,None),)]*0
        else:
            tmpData[(slice(None),)*(dim-1)+(slice(None,shift),)] = tmpData[(slice(None),)*(dim-1)+(slice(None,shift),)]*0
        if self.spec > 0:
           tmpData=self.fourierLocal(tmpData,0)
        self.showFid(tmpData)
        
    def getdcOffset(self,pos1,pos2):
        minPos = int(min(pos1,pos2))
        maxPos = int(max(pos1,pos2))
        if minPos != maxPos:
            return np.mean(self.data1D[(len(self.data1D.shape)-1)*(slice(None),)+(slice(minPos,maxPos),)])
        else:
            return 0
            
    def dcOffset(self,offset):
        self.showFid(self.data1D-offset)
            
    def applydcOffset(self,offset):
        returnValue = self.data.dcOffset(offset)
        self.upd()
        self.showFid()
        self.root.addMacro(['offset',(offset,)])
        return returnValue
    
    def applyLPSVD(self):
        returnValue = self.data.LPSVD(self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['lpsvd',(self.axes,)])
        return returnValue
    
    def states(self):
        returnValue = self.data.states(self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['states',(self.axes,)])
        return returnValue
    
    def statesTPPI(self):
        returnValue = self.data.statesTPPI(self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['statesTPPI',(self.axes,)])
        return returnValue
    
    def integrate(self,pos1,pos2):
        self.root.addMacro(['integrate',(pos1,pos2,self.axes,)])
        return self.data.matrixManip(pos1,pos2,self.axes,0)

    def maxMatrix(self,pos1,pos2):
        self.root.addMacro(['max',(pos1,pos2,self.axes,)])
        return self.data.matrixManip(pos1,pos2,self.axes,1)
    
    def minMatrix(self,pos1,pos2):
        self.root.addMacro(['min',(pos1,pos2,self.axes,)])
        return self.data.matrixManip(pos1,pos2,self.axes,2)

    def flipLR(self):
        returnValue = self.data.flipLR(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['fliplr',(self.axes,)])
        return returnValue

    def concatenate(self,axes):
        returnValue = self.data.concatenate(axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['concatenate',(axes,)])
        return returnValue

    def split(self,sections):
        returnValue = self.data.split(sections, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['split',(sections,self.axes)])
        return returnValue
    
    def insert(self,data,pos):
        returnValue = self.data.insert(data,pos,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['insert',(data,pos,self.axes)])
        return returnValue
    
    def delete(self,pos):
        returnValue = self.data.remove(pos,self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['delete',(pos,self.axes)])
        return returnValue

    def deletePreview(self,pos):
        self.data1D = np.delete(self.data1D,pos,axis=len(self.data1D.shape)-1)
        self.xax = np.delete(self.xax,pos)
        if (np.array(self.data1D.shape) != 0).all():
            self.showFid()
        self.upd()
    
    def add(self,data):
        returnValue = self.data.add(data)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['add',(data,)])
        return returnValue
        
    def subtract(self,data):
        returnValue = self.data.subtract(data)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['subtract',(data,)])
        return returnValue
    
    def getRegion(self,pos1,pos2): #set the frequency of the actual data
        returnValue = self.data.getRegion(pos1,pos2,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        return returnValue
    
    def shearing(self,shear,axes,axes2):
        returnValue = self.data.shear(shear,axes,axes2)
        self.upd()
        self.showFid()
        self.root.addMacro(['shear',(shear,axes,axes2)])
        return returnValue
    
    def ACMEentropy(self,phaseIn,phaseAll=True):
        if len(self.data1D.shape) > 1:
            tmp = self.data1D[0]
        else:
            tmp = self.data1D
        phase0=phaseIn[0]
        if phaseAll:
            phase1=phaseIn[1]
        else:
            phase1=0.0
        L = len(tmp)
        if self.spec>0:
            x=np.fft.fftshift(np.fft.fftfreq(L,1.0/self.sw))/self.sw
            s0 = tmp*np.exp(1j*(phase0+phase1*x))
        else:
            s0 = np.fft.fftshift(np.fft.fft(tmp))*np.exp(1j*(phase0+phase1*x))
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
        self.upd()
        if phaseNum == 0:
            phases = scipy.optimize.minimize(self.ACMEentropy,[0],(False,),method='Powell')
            phases = [phases['x']]
        elif phaseNum == 1:
            phases = scipy.optimize.minimize(self.ACMEentropy,[0,0],method='Powell')
            phases = phases['x']
        return phases

    def setXaxPreview(self,xax):
        self.xax = xax
        self.plotReset()
        self.showFid()
        self.upd()

    def setXax(self,xax):
        returnVal = self.data.setXax(xax,self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['setxax',(xax,self.axes)])
        return returnVal

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
        self.xminlim = (self.xminlim - oldAxAdd) * newAxMult / oldAxMult + newAxAdd
        self.xmaxlim = (self.xmaxlim - oldAxAdd) * newAxMult / oldAxMult + newAxAdd
        self.showFid()

    def hilbert(self):
        returnValue = self.data.hilbert(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['hilbert',(self.axes,)])
        return returnValue
    
    def getDisplayedData(self):
        if len(self.data1D.shape) > 1:
            tmp = self.data1D[0]
        else:
            tmp = self.data1D
        if self.plotType==0:
            return np.real(tmp)
        elif self.plotType==1:
            return np.imag(tmp)
        elif self.plotType==2:
            return np.real(tmp)
        elif self.plotType==3:
            return np.abs(tmp)      

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False,output=None): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        self.ax.cla()
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
                if self.single:
                    self.ax.scatter(self.xax*axMult+axAdd,np.real(self.data1D),c='k',alpha=0.2)
                else:
                    self.ax.plot(self.xax*axMult+axAdd,np.real(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==1):
                if self.single:
                    self.ax.scatter(self.xax*axMult+axAdd,np.imag(self.data1D),c='k',alpha=0.2)
                else:
                    self.ax.plot(self.xax*axMult+axAdd,np.imag(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==2):
                if self.single:
                    self.ax.scatter(self.xax*axMult+axAdd,np.real(self.data1D),c='k',alpha=0.2)
                else:
                    self.ax.plot(self.xax*axMult+axAdd,np.real(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==3):
                if self.single:
                    self.ax.scatter(self.xax*axMult+axAdd,np.abs(self.data1D),c='k',alpha=0.2)
                else:
                    self.ax.plot(self.xax*axMult+axAdd,np.abs(self.data1D),c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraX)):
                if self.single:
                    self.ax.scatter(extraX[num]*axMult+axAdd,extraY[num],c=extraColor[num])
                else:
                    self.ax.plot(extraX[num]*axMult+axAdd,extraY[num],c=extraColor[num])
        if (self.plotType==0):
            if self.single:
                self.line = self.ax.scatter(self.xax*axMult+axAdd,np.real(tmpdata),c='b')
            else:
                self.line = self.ax.plot(self.xax*axMult+axAdd,np.real(tmpdata),c='b')
        elif(self.plotType==1):
            if self.single:
                self.line =self.ax.scatter(self.xax*axMult+axAdd,np.imag(tmpdata),c='b')
            else:
                self.line =self.ax.plot(self.xax*axMult+axAdd,np.imag(tmpdata),c='b')
        elif(self.plotType==2):
            if self.single:
                self.ax.scatter(self.xax*axMult+axAdd,np.imag(tmpdata),c='r')
                self.line = self.ax.scatter(self.xax*axMult+axAdd,np.real(tmpdata),c='b')
            else:
                self.ax.plot(self.xax*axMult+axAdd,np.imag(tmpdata),c='r')
                self.line = self.ax.plot(self.xax*axMult+axAdd,np.real(tmpdata),c='b')
        elif(self.plotType==3):
            if self.single:
                self.line =self.ax.scatter(self.xax*axMult+axAdd,np.abs(tmpdata),c='b')
            else:
                self.line =self.ax.plot(self.xax*axMult+axAdd,np.abs(tmpdata),c='b')
        if self.spec==0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        if output is not None:
            self.canvas.print_figure(output)
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
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
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        
#########################################################################################################
#the class from which the stacked data is displayed, the operations which only edit the content of this class are for previewing
class CurrentStacked(Current1D):
    def __init__(self, root, data, duplicateCurrent=None):
        self.data = data
        if hasattr(duplicateCurrent,'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent,'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent,'stackBegin'):
            self.stackBegin = duplicateCurrent.stackBegin
        else:
            self.stackBegin = None
        if hasattr(duplicateCurrent,'stackEnd'):
            self.stackEnd = duplicateCurrent.stackEnd
        else:
            self.stackEnd = None
        if hasattr(duplicateCurrent,'stackStep'):
            self.stackStep = duplicateCurrent.stackStep
        else:
            self.stackStep = None
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.spacing = 0
        Current1D.__init__(self,root, data, duplicateCurrent)
        self.resetSpacing()
        self.plotReset()
        self.showFid()

    def copyCurrent(self,root, data):
        return CurrentStacked(root,data,self)
        
    def upd(self): #get new data from the data instance
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        updateVar = self.data.getBlock(self.axes,self.axes2,self.locList,self.stackBegin, self.stackEnd, self.stackStep)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.freq2 = updateVar[2]
        self.sw = updateVar[3]
        self.sw2 = updateVar[4]
        self.spec = updateVar[5]
        self.spec2 = updateVar[6]
        self.wholeEcho = updateVar[7]
        self.wholeEcho2 = updateVar[8]
        self.xax=updateVar[9]
        self.xax2=updateVar[10]
        self.ref=updateVar[11]
        self.ref2=updateVar[12]
        if self.ref is None:
            self.ref = self.freq
        if self.ref2 is None:
            self.ref2 = self.freq2
        self.single = self.data1D.shape[-1]==1

    def setBlock(self,axes,axes2,locList,stackBegin=None,stackEnd=None,stackStep=None): #change the slice 
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if stackStep == None:
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
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
                        x2=x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
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
                    x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
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
                x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
            x = np.repeat([x],len(self.data1D),axis=0)
        y = self.data1D
        self.ax.cla()
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

    def setSpacing(self, spacing):
        self.spacing = spacing
        self.showFid()

    def resetSpacing(self):
        difference = np.diff(self.data1D,axis=0)
        if difference.size == 0:
            self.spacing = 0
        else:
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

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        self.ax.cla()
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
                    if self.single:
                        self.ax.scatter(self.xax*axMult+axAdd,num*self.spacing+np.real(self.data1D[num]),c='k',alpha=0.2)
                    else:
                        self.ax.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==1):
                for num in range(len(self.data1D)):
                    if self.single:
                        self.ax.scatter(self.xax*axMult+axAdd,num*self.spacing+np.imag(self.data1D[num]),c='k',alpha=0.2)
                    else:
                        self.ax.plot(self.xax*axMult+axAdd,num*self.spacing+np.imag(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==2):
                for num in range(len(self.data1D)):
                    if self.single:
                        self.ax.scatter(self.xax*axMult+axAdd,num*self.spacing+np.real(self.data1D[num]),c='k',alpha=0.2)
                    else:
                        self.ax.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==3):
                for num in range(len(self.data1D)):
                    if self.single:
                        self.ax.scatter(self.xax*axMult+axAdd,num*self.spacing+np.abs(self.data1D[num]),c='k',alpha=0.2)
                    else:
                        self.ax.plot(self.xax*axMult+axAdd,num*self.spacing+np.abs(self.data1D[num]),c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraY)):
                if self.single:
                    self.ax.scatter(extraX[0]*axMult+axAdd,num*self.spacing+extraY[num],c=extraColor[0])
                else:
                    self.ax.plot(extraX[0]*axMult+axAdd,num*self.spacing+extraY[num],c=extraColor[0])
        if (self.plotType==0):
            tmpdata = np.real(tmpdata)
        elif (self.plotType==1):
            tmpdata = np.imag(tmpdata)
        elif(self.plotType==3):
            tmpdata = np.abs(tmpdata)
        if self.single:
            for num in range(len(tmpdata)):
                if (self.plotType==2):
                    self.ax.scatter(self.xax*axMult+axAdd,num*self.spacing+np.imag(tmpdata[num]),c='r')
                if num is 0:
                    self.line = self.ax.scatter(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')
                else:
                    self.ax.scatter(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')
        else:
            for num in range(len(tmpdata)):
                if (self.plotType==2):
                    self.ax.plot(self.xax*axMult+axAdd,num*self.spacing+np.imag(tmpdata[num]),c='r')
                if num is 0:
                    self.line = self.ax.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')
                else:
                    self.ax.plot(self.xax*axMult+axAdd,num*self.spacing+np.real(tmpdata[num]),c='b')  
        if self.spec==0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
        self.ax=self.fig.gca()
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
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)

#########################################################################################################
#the class from which the arrayed data is displayed, the operations which only edit the content of this class are for previewing
class CurrentArrayed(Current1D):
    def __init__(self, root, data, duplicateCurrent = None):
        self.data = data
        if hasattr(duplicateCurrent,'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent,'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent,'stackBegin'):
            self.stackBegin = duplicateCurrent.stackBegin
        else:
            self.stackBegin = None
        if hasattr(duplicateCurrent,'stackEnd'):
            self.stackEnd = duplicateCurrent.stackEnd
        else:
            self.stackEnd = None
        if hasattr(duplicateCurrent,'stackStep'):
            self.stackStep = duplicateCurrent.stackStep
        else:
            self.stackStep = None
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.spacing = 0
        Current1D.__init__(self, root, data, duplicateCurrent)
        self.resetSpacing()
        self.plotReset()
        self.showFid()

    def copyCurrent(self,root, data):
        return CurrentArrayed(root,data,self)
        
    def upd(self): #get new data from the data instance
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        updateVar = self.data.getBlock(self.axes,self.axes2,self.locList,self.stackBegin, self.stackEnd, self.stackStep)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.freq2 = updateVar[2]
        self.sw = updateVar[3]
        self.sw2 = updateVar[4]
        self.spec = updateVar[5]
        self.spec2 = updateVar[6]
        self.wholeEcho = updateVar[7]
        self.wholeEcho2 = updateVar[8]
        self.xax=updateVar[9]
        self.xax2=updateVar[10]
        self.ref=updateVar[11]
        self.ref2=updateVar[12]
        if self.ref is None:
            self.ref = self.freq
        if self.ref2 is None:
            self.ref2 = self.freq2
        self.single = self.data1D.shape[-1]==1
 
    def setBlock(self,axes,axes2,locList,stackBegin=None,stackEnd=None,stackStep=None): #change the slice 
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if self.data.data.shape[self.axes2] > 100:
            self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
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
                        x2=x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
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
        self.ax.cla()
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

    def setSpacing(self, spacing):
        self.spacing = spacing
        self.showFid()

    def resetSpacing(self):
        self.spacing = (self.xax[-1]-self.xax[0])*1.1     

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        self.ax.cla()
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
                oldData = np.real(self.data1D)
            elif(self.plotType==1):
                oldData = np.imag(self.data1D)
            elif(self.plotType==2):
                oldData = np.real(self.data1D)
            elif(self.plotType==3):
                oldData = np.abs(self.data1D)
            for num in range(len(self.data1D)):
                if self.single:
                    self.ax.scatter((num*self.spacing+self.xax)*axMult+axAdd,oldData[num][direc],c='k',alpha=0.2)
                else:
                    self.ax.plot((num*self.spacing+self.xax)*axMult+axAdd,oldData[num][direc],c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraY)):
                if self.single:
                    self.ax.scatter((num*self.spacing+extraX[0])*axMult+axAdd,extraY[num][direc],c=extraColor[0])
                else:
                    self.ax.plot((num*self.spacing+extraX[0])*axMult+axAdd,extraY[num][direc],c=extraColor[0])
        self.line = []
        if (self.plotType==0):
            tmpdata = np.real(tmpdata)
        elif(self.plotType==1):
            tmpdata = np.imag(tmpdata)
        elif(self.plotType==3):
            tmpdata = np.abs(tmpdata)
        if self.single:
            for num in range(len(tmpdata)):
                if (self.plotType==2):
                    self.ax.scatter((num*self.spacing+self.xax)*axMult,np.imag(tmpdata[num])[direc],c='r')
                self.line.append(self.ax.scatter((num*self.spacing+self.xax)*axMult,np.real(tmpdata[num])[direc],c='b')[0])
        else:
            for num in range(len(tmpdata)):
                if (self.plotType==2):
                    self.ax.plot((num*self.spacing+self.xax)*axMult,np.imag(tmpdata[num])[direc],c='r')
                self.line.append(self.ax.plot((num*self.spacing+self.xax)*axMult,np.real(tmpdata[num])[direc],c='b')[0])
        
        if self.spec==0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
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
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)

#########################################################################################################
#the class from which the contour data is displayed, the operations which only edit the content of this class are for previewing
class CurrentContour(Current1D):
    def __init__(self, root, data, duplicateCurrent=None):
        self.data = data
        if hasattr(duplicateCurrent,'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent,'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent,'axType2'):
            self.axType2 = duplicateCurrent.axType2
        else:
            self.axType2 = 1
        if hasattr(duplicateCurrent,'ppm2'):
            self.ppm2 = duplicateCurrent.ppm2
        else:
            self.ppm2 = False
        if hasattr(duplicateCurrent,'numLevels'):
            self.numLevels = duplicateCurrent.numLevels
        else:
            self.numLevels = 20
        if hasattr(duplicateCurrent,'minLevels'):
            self.minLevels = duplicateCurrent.minLevels
        else:
            self.minLevels = 0.0
        if hasattr(duplicateCurrent,'maxLevels'):
            self.maxLevels = duplicateCurrent.maxLevels
        else:
            self.maxLevels = 1.0
        Current1D.__init__(self, root, data, duplicateCurrent)
        self.plotReset()
        self.showFid()

    def copyCurrent(self,root, data):
        return CurrentContour(root,data,self)
        
    def upd(self): #get new data from the data instance
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        updateVar = self.data.getBlock(self.axes,self.axes2,self.locList)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.freq2 = updateVar[2]
        self.sw = updateVar[3]
        self.sw2 = updateVar[4]
        self.spec = updateVar[5]
        self.spec2 = updateVar[6]
        self.wholeEcho = updateVar[7]
        self.wholeEcho2 = updateVar[8]
        self.xax=updateVar[9]
        self.xax2=updateVar[10]
        self.ref=updateVar[11]
        self.ref2=updateVar[12]
        if self.ref is None:
            self.ref = self.freq
        if self.ref2 is None:
            self.ref2 = self.freq2
        self.single = self.data1D.shape[-1]==1

    def setBlock(self,axes,axes2,locList): #change the slice 
        self.axes = axes
        self.axes2 = axes2
        self.locList = locList
        self.upd()
        self.plotReset()
        self.showFid()

    def setLevels(self,numLevels, maxLevels, minLevels):
        self.numLevels = numLevels
        self.maxLevels = maxLevels
        self.minLevels = minLevels
        self.showFid()
        
    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)

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
                        x2=x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
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
                    x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
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
                x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
            x = np.repeat([x],len(self.data1D),axis=0)
        y = self.data1D
        self.ax.cla()
        if self.spec ==1:
            y=np.fft.ifftn(np.fft.ifftshift(y,axes=1),axes=[1])
            y= y*x
            y=np.fft.fftshift(np.fft.fftn(y,axes=[1]),axes=1)
        else:
            y= y*x
        self.showFid(y)

    def showFid(self, tmpdata=None): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        self.ax.cla()
        self.x_ax.cla()
        self.y_ax.cla()
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        axAdd2 = 0
        if self.spec2 == 1:
            if self.ppm2:
                axAdd2 = (self.freq2-self.ref2)/self.ref2*1e6
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        x=self.xax*axMult+axAdd
        y=self.xax2*axMult2+axAdd2
        X, Y = np.meshgrid(x,y)
        self.line = []
        if (self.plotType==0):
            differ=np.amax(np.real(tmpdata))-np.amin(np.real(tmpdata))
            contourLevels = np.linspace(self.minLevels*differ+np.amin(np.real(tmpdata)),self.maxLevels*differ+np.amin(np.real(tmpdata)),self.numLevels)
            self.line.append(self.ax.contour(X, Y, np.real(tmpdata),c='b',levels=contourLevels))
            self.x_ax.plot(x,np.amax(np.real(tmpdata),axis=0),'b')
            self.y_ax.plot(np.amax(np.real(tmpdata),axis=1),y,'b')
        elif(self.plotType==1):
            differ=np.amax(np.imag(tmpdata))-np.amin(np.imag(tmpdata))
            contourLevels = np.linspace(self.minLevels*differ+np.amin(np.imag(tmpdata)),self.maxLevels*differ+np.amin(np.imag(tmpdata)),self.numLevels)
            self.line.append(self.ax.contour(X, Y, np.imag(tmpdata),c='b',levels=contourLevels))
            self.x_ax.plot(x,np.amax(np.imag(tmpdata),axis=0),'b')
            self.y_ax.plot(np.amax(np.imag(tmpdata),axis=1),y,'b')
        elif(self.plotType==2):
            print('type not supported')
            differ=np.amax(np.real(tmpdata))-np.amin(np.real(tmpdata))
            contourLevels = np.linspace(self.minLevels*differ+np.amin(np.real(tmpdata)),self.maxLevels*differ+np.amin(np.real(tmpdata)),self.numLevels)
            self.line.append(self.ax.contour(X, Y, np.real(tmpdata),c='b',levels=contourLevels))
            self.x_ax.plot(x,np.amax(np.real(tmpdata),axis=0),'b')
            self.y_ax.plot(np.amax(np.real(tmpdata),axis=1),y,'b')
        elif(self.plotType==3):
            differ=np.amax(np.abs(tmpdata))-np.amin(np.abs(tmpdata))
            contourLevels = np.linspace(self.minLevels*differ+np.amin(np.abs(tmpdata)),self.maxLevels*differ+np.amin(np.abs(tmpdata)),self.numLevels)
            self.line.append(self.ax.contour(X, Y, np.abs(tmpdata),c='b',levels=contourLevels))
            self.x_ax.plot(x,np.amax(np.abs(tmpdata),axis=0),'b')
            self.y_ax.plot(np.amax(np.abs(tmpdata),axis=1),y,'b')
        if self.spec==0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        if self.spec2==0:
            if self.axType2 == 0:
                self.ax.set_ylabel('Time [s]')
            elif self.axType2 == 1:
                self.ax.set_ylabel('Time [ms]')
            elif self.axType2 == 2:
                self.ax.set_ylabel(r'Time [$\mu$s]')
            else:
                self.ax.set_ylabel('User defined')
        elif self.spec2==1:
            if self.ppm2:
                self.ax.set_ylabel('Frequency [ppm]')
            else:
                if self.axType2 == 0:
                    self.ax.set_ylabel('Frequency [Hz]')
                elif self.axType2 == 1:
                    self.ax.set_ylabel('Frequency [kHz]')
                elif self.axType2 == 2:
                    self.ax.set_ylabel('Frequency [MHz]')
                else:
                    self.ax.set_ylabel('User defined')
        else:
            self.ax.set_ylabel('')
        if self.spec:
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim,self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
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
        axAdd2 = 0
        if self.spec2 == 1:
            if self.ppm2:
                axAdd2 = (self.freq2-self.ref2)/self.ref2*1e6
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        self.yminlim=min(self.xax2*axMult2+axAdd2)
        self.ymaxlim=max(self.xax2*axMult2+axAdd2)
        if self.spec:
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim,self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim,self.ymaxlim)

    #The peakpicking function needs to be changed for contour plots
    def buttonRelease(self,event):
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
                        self.ax.set_xlim(self.xmaxlim,self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim,self.xmaxlim)
                    if self.spec2 > 0:
                        self.ax.set_ylim(self.ymaxlim,self.yminlim)
                    else:
                        self.ax.set_ylim(self.yminlim,self.ymaxlim)
                self.zoomX1=None
                self.zoomX2=None #WF: should also be cleared, memory of old zoom
                self.zoomY1=None
                self.zoomY2=None #WF: should also be cleared, memory of old zoom
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()


#########################################################################################################
#The skewed plot class
class CurrentSkewed(Current1D):
    def __init__(self, root, data, duplicateCurrent=None):
        self.data = data
        if hasattr(duplicateCurrent,'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent,'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent,'stackBegin'):
            self.stackBegin = duplicateCurrent.stackBegin
        else:
            self.stackBegin = None
        if hasattr(duplicateCurrent,'stackEnd'):
            self.stackEnd = duplicateCurrent.stackEnd
        else:
            self.stackEnd = None
        if hasattr(duplicateCurrent,'stackStep'):
            self.stackStep = duplicateCurrent.stackStep
        else:
            self.stackStep = None
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        if hasattr(duplicateCurrent,'axType2'):
            self.axType2 = duplicateCurrent.axType2
        else:
            self.axType2 = 1
        if hasattr(duplicateCurrent,'ppm2'):
            self.ppm2 = duplicateCurrent.ppm2
        else:
            self.ppm2 = False
        Current1D.__init__(self, root, data, duplicateCurrent)
        self.plotReset()
        self.setSkewed(-0.2,70)
        
    def copyCurrent(self,root, data):
        return CurrentSkewed(root,data,self)
    
    def upd(self): #get new data from the data instance
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        updateVar = self.data.getBlock(self.axes,self.axes2,self.locList,self.stackBegin, self.stackEnd, self.stackStep)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.freq2 = updateVar[2]
        self.sw = updateVar[3]
        self.sw2 = updateVar[4]
        self.spec = updateVar[5]
        self.spec2 = updateVar[6]
        self.wholeEcho = updateVar[7]
        self.wholeEcho2 = updateVar[8]
        self.xax=updateVar[9]
        self.xax2=updateVar[10]
        self.ref=updateVar[11]
        self.ref2=updateVar[12]
        if self.ref is None:
            self.ref = self.freq
        if self.ref2 is None:
            self.ref2 = self.freq2
        self.single = self.data1D.shape[-1]==1
 
    def setBlock(self,axes,axes2,locList,stackBegin=None,stackEnd=None,stackStep=None): #change the slice 
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if self.data.data.shape[self.axes2] > 100:
            self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.locList = locList
        self.upd()
        self.plotReset()
        self.showFid()

    def setSkewed(self,skewed,elevation):
        self.skewed = skewed
        self.elevation = elevation
        proj3d.persp_transformation = lambda zfront, zback : np.array([[1,0,0,0],
                                                                       [skewed ,1.0,0,0],
                                                                       [0,0,zfront,0],
                                                                       [0,0,-0.00001,zback]])
        self.ax.view_init(elev=self.elevation, azim=180*np.arctan(skewed/np.sin(np.pi*elevation/180))/np.pi-90)
        self.showFid()
        
    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)
        
    def stackSelect(self,stackBegin, stackEnd, stackStep):
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.upd()
        self.showFid()

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
                        x2=x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
                    if self.wholeEcho:
                        x2[-1:-(len(x2)/2+1):-1]=x2[:len(x2)/2]
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
                    x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
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
                x=x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0,np.pi,len(self.data1D[0])))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
            x = np.repeat([x],len(self.data1D),axis=0)
        y = self.data1D
        self.ax.cla()
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

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        self.ax.cla()
        axAdd = 0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        axAdd2 = 0
        if self.spec2 == 1:
            if self.ppm2:
                axAdd2 = (self.freq2-self.ref2)/self.ref2*1e6
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        x=self.xax*axMult+axAdd
        y=self.xax2*axMult2+axAdd2
        if old:
            if (self.plotType==0):
                for num in range(len(self.data1D)):
                    self.ax.plot(x,y[num]*np.ones(len(x)),np.real(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==1):
                for num in range(len(self.data1D)):
                    self.ax.plot(x,y[num]*np.ones(len(x)),np.imag(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==2):
                for num in range(len(self.data1D)):
                    self.ax.plot(x,y[num]*np.ones(len(x)),np.real(self.data1D[num]),c='k',alpha=0.2)
            elif(self.plotType==3):
                for num in range(len(self.data1D)):
                    self.ax.plot(x,y[num]*np.ones(len(x)),np.abs(self.data1D[num]),c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraY)):
                self.ax.plot(extraX[0]*axMult+axAdd,y[num]*np.ones(len(extraX[0])),extraY[num],c=extraColor[0])

        self.line = []
        if (self.plotType==0):
            for num in range(len(tmpdata)):
                self.line.append(self.ax.plot(x,y[num]*np.ones(len(x)),np.real(tmpdata[num]),c='b'))
        elif(self.plotType==1):
            for num in range(len(tmpdata)):
                self.line.append(self.ax.plot(x,y[num]*np.ones(len(x)),np.imag(tmpdata[num]),c='b'))
        elif(self.plotType==2):
            for num in range(len(tmpdata)):
                self.ax.plot(x,y[num]*np.ones(len(x)),np.imag(tmpdata[num]),c='r')
                self.line.append(self.ax.plot(x,y[num]*np.ones(len(x)),np.real(tmpdata[num]),c='b'))
        elif(self.plotType==3):
            for num in range(len(tmpdata)):
                self.line.append(self.ax.plot(x,y[num]*np.ones(len(x)),np.abs(tmpdata[num]),c='b'))
        if self.spec==0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        if self.spec2==0:
            if self.axType2 == 0:
                self.ax.set_ylabel('Time [s]')
            elif self.axType2 == 1:
                self.ax.set_ylabel('Time [ms]')
            elif self.axType2 == 2:
                self.ax.set_ylabel(r'Time [$\mu$s]')
            else:
                self.ax.set_ylabel('User defined')
        elif self.spec2==1:
            if self.ppm2:
                self.ax.set_ylabel('Frequency [ppm]')
            else:
                if self.axType2 == 0:
                    self.ax.set_ylabel('Frequency [Hz]')
                elif self.axType2 == 1:
                    self.ax.set_ylabel('Frequency [kHz]')
                elif self.axType2 == 2:
                    self.ax.set_ylabel('Frequency [MHz]')
                else:
                    self.ax.set_ylabel('User defined')
        else:
            self.ax.set_ylabel('')
        if self.spec:
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim,self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.w_zaxis.line.set_lw(0.)
        self.ax.set_zticks([])
        self.ax.grid(False)
        self.ax.xaxis.pane.set_edgecolor('white')
        self.ax.yaxis.pane.set_edgecolor('white')
        self.ax.zaxis.pane.set_edgecolor('white')
        self.ax.xaxis.pane.fill = False
        self.ax.yaxis.pane.fill = False
        self.ax.zaxis.pane.fill = False
        self.canvas.draw()

    def plotReset(self): #set the plot limits to min and max values
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
        axAdd2 = 0
        if self.spec2 == 1:
            if self.ppm2:
                axAdd2 = (self.freq2-self.ref2)/self.ref2*1e6
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        self.yminlim=min(self.xax2*axMult2+axAdd2)
        self.ymaxlim=max(self.xax2*axMult2+axAdd2)
        if self.spec:
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim,self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim,self.ymaxlim)
