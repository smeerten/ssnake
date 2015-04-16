import matplotlib
import numpy as np
import scipy.optimize
import copy

#########################################################################
#the generic data class
class Spectrum(object):
    def __init__(self, data, freq, sw , spec=None, wholeEcho=None):
        self.dim = len(data.shape)                    #number of dimensions
        self.data = np.array(data,dtype=complex)      #data of dimension dim
        self.freq = freq                              #array of center frequency (length is dim, MHz)
        self.sw = sw                                  #array of sweepwidths
        if spec is None:
            self.spec=[False]*self.dim
        else:
            self.spec = spec                              #boolean array of length dim where False=time domain and True=spectral domain
        if wholeEcho is None:
            self.wholeEcho = [False]*self.dim
        else:
            self.wholeEcho = wholeEcho                    #boolean array of length dim where True indicates a full Echo


    def setPhase(self, phase0, phase1, axes):
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.data.shape[axes],1.0/self.sw[axes]))*phase1*1j)
        if not(self.spec[axes]):
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

    def apodize(self,lor,gauss, cos2, axes):
        copyData=copy.deepcopy(self)
        axLen = self.data.shape[axes]
        t=np.arange(0,axLen)/(2.0*self.sw[axes])
        x=np.ones(axLen)
        if lor is not None:
            x=x*np.exp(-lor*t)
        if gauss is not None:
            x=x*np.exp(-(gauss*t)**2)
        if cos2 is not None:
            x=x*(np.cos(cos2*np.linspace(0,0.5*np.pi,len(self.data1D)))**2)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.apodize(lor,gauss,axes))
        if self.spec[axes]:
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
        return returnValue

    def changeSpec(self,axes):
        if self.spec[axes]:
            self.spec[axes] = False
        else:
            self.spec[axes] = True
        return lambda self: self.changeSpec(axes)

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
        if not self.spec[axes]:
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None),) * axes + (0,) + (slice(None),)*(self.dim-1-axes)
                self.data[slicing]= self.data[slicing]*0.5
            self.data=np.fft.fftshift(np.fft.fftn(self.data,axes=[axes]),axes=axes)
            self.spec[axes]=True
        else:
            self.data=np.fft.ifftn(np.fft.ifftshift(self.data,axes=axes),axes=[axes])
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None),) * axes + (0,) + (slice(None),)*(self.dim-1-axes)
                self.data[slicing]= self.data[slicing]*2.0
            self.spec[axes]=False
        return lambda self: self.fourier(axes)

    def getSlice(self,axes,locList):
        return (self.data[tuple(locList[:axes])+(slice(None),)+tuple(locList[axes:])],self.freq[axes],self.sw[axes],self.spec[axes],self.wholeEcho[axes])

    def restoreData(self,copyData,returnValue): # restore data from an old copy for undo purposes
        self.data = copyData.data
        self.dim = len(self.data.shape)                    #number of dimensions
        self.freq = copyData.freq                              #array of center frequency (length is dim, MHz)
        self.sw = copyData.sw                                  #array of sweepwidths
        self.spec = copyData.spec
        self.wholeEcho = copyData.wholeEcho
        return returnValue
        


#########################################################################################################
#the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing
class Current1D(object):
    def __init__(self, data, axes, locList, plotType , fig, canvas):
        self.data = data              #the actual spectrum instance
        self.freq = None              #frequency of the slice 
        self.sw = None                #x-data display
        self.data1D = None            #the data1D
        self.spec = None              #boolean where False=time domain and True=spectral domain
        self.wholeEcho = None
        self.axes = axes              #dimension of which the data is displayed
        self.locList = locList        #list of length dim-1 with the matrix coordinates of current spectrum 
        self.fig = fig                #figure
        self.canvas = canvas          #canvas
        self.mainLine = None          #keep the main line of the plot
        self.leftMouse = False        #is the left mouse button currently pressed
        self.panX = None              #start position of dragging the spectrum
        self.panY = None              #start position of dragging the spectrum 
        self.zoomX1 = None            #first corner of the zoombox 
        self.zoomY1 = None            #first corner of the zoombox
        self.zoomX2 = None            #second corner of the zoombox
        self.zoomY2 = None            #second corner of the zoombox
        self.rect=[None,None,None,None]      #lines for zooming or peak picking
        self.rightMouse = False              #is the right mouse button currently pressed
        self.peakPick = False         #currently peakPicking
        self.peakPickFunc = None      #the function that needs to be called after peakPicking
        #connect click events to the canvas
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        self.plotType = plotType      #0=real,1=imag,2=both 3=abs ...
        self.update()   #get the first slice of data
        self.plotReset() #reset the axes limits
        self.showFid() #plot the data
        
    def update(self): #get new data from the data instance
        updateVar = self.data.getSlice(self.axes,self.locList)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.sw = updateVar[2]
        self.spec = updateVar[3]
        self.wholeEcho = updateVar[4]

    def setSlice(self,axes,locList): #change the slice 
        axesSame = True
        if self.axes != axes:
            axesSame = False
            self.axes = axes
        self.locList = locList
        self.update()
        if not axesSame:
            self.plotReset()
        self.showFid()

    def setPhaseInter(self, phase0in, phase1in): #interactive changing the phase without editing the actual data
        phase0=float(phase0in)
        phase1=float(phase1in)
        if not(self.spec):
            tmpdata=self.fourierLocal(self.data1D,self.spec)
            tmpdata=tmpdata*np.exp(phase0*1j)
            tmpdata=tmpdata*np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata),1.0/self.sw))*phase1*1j)
            tmpdata=self.fourierLocal(tmpdata,not self.spec)
        else:
            tmpdata=self.data1D*np.exp(phase0*1j)
            tmpdata=tmpdata*np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata),1.0/self.sw))*phase1*1j)
        self.showFid(tmpdata)

    def applyPhase(self, phase0, phase1):# apply the phase to the actual data
        phase0=float(phase0)
        phase1=float(phase1)
        returnValue = self.data.setPhase(phase0,phase1,self.axes)
        self.update()
        self.showFid()
        return returnValue

    def fourier(self): #fourier the actual data and replot
        returnValue = self.data.fourier(self.axes)
        self.update()
        self.plotReset()
        self.showFid()
        return returnValue

    def fourierLocal(self, fourData, spec): #fourier the local data for other functions
        if not(spec):
            fourData=np.fft.fftshift(np.fft.fft(fourData))
        else:
            fourData=np.fft.ifft(np.fft.ifftshift(fourData))
        return fourData

    def apodPreview(self,lor=None,gauss=None, cos2=None): #display the 1D data including the apodization function
        t=np.arange(0,len(self.data1D))/(2.0*self.sw)
        x=np.ones(len(self.data1D))
        if lor is not None:
            x=x*np.exp(-lor*t)
        if gauss is not None:
            x=x*np.exp(-(gauss*t)**2)
        if cos2 is not None:
            x=x*(np.cos(cos2*np.linspace(0,0.5*np.pi,len(self.data1D)))**2)
        if self.wholeEcho:
            x[-1:-(len(x)/2+1):-1]=x[:len(x)/2]
        a=self.fig.gca()
        a.cla()
        y = self.data1D
        if self.spec:
            y=np.fft.ifft(np.fft.ifftshift(y))
            y= y*x
            y=np.fft.fftshift(np.fft.fft(y))
        else:
            y= y*x
        if not(self.spec):

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

    def applyApod(self,lor=None,gauss=None,cos2=None): #apply the apodization to the actual data
        returnValue = self.data.apodize(lor,gauss,cos2,self.axes)
        self.update() 
        self.showFid()     
        return returnValue

    def setFreq(self,freq,sw): #set the frequency of the actual data
        returnValue = self.data.setFreq(freq,sw,self.axes)
        self.update()
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
        self.plotReset()
        self.showFid()
        self.update()

    def applySize(self,size): #set size to the actual data
        returnValue = self.data.setSize(size,self.axes)
        self.update()
        self.plotReset()
        self.showFid()
        return returnValue

    def changeSpec(self): #change from time to freq domain of the actual data
        returnValue = self.data.changeSpec(self.axes)
        self.update()
        self.plotReset()
        self.showFid()
        return returnValue

    def applySwapEcho(self,idx):
        returnValue = self.data.swapEcho(idx,self.axes)
        self.update()
        self.showFid()
        return returnValue

    def setSwapEchoPreview(self,idx):
        self.data1D = np.concatenate((self.data1D[idx:],self.data1D[:idx]))
        self.plotReset()
        self.showFid()
        self.update()

    def setWholeEcho(self, value):
        if value == 0:
            self.data.wholeEcho[self.axes]=False
            self.wholeEcho = False
        else:
            self.data.wholeEcho[self.axes]=True
            self.wholeEcho = True

    def applyShift(self,shift):
        returnValue = self.data.shiftData(shift,self.axes)
        self.update()
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
        self.update()
        self.showFid()
        return returnValue

    def ACMEentropy(self,phaseIn,phaseAll=True):
        phase0=phaseIn[0]
        if phaseAll:
            phase1=phaseIn[1]
        else:
            phase1=0.0
        L = len(self.data1D)
        x=np.fft.fftshift(np.fft.fftfreq(L,1.0/self.sw))
        if self.spec:
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

    def plotReset(self): #set the plot limits to min and max values
        a=self.fig.gca()
        if self.plotType==0:
            self.yminlim=min(np.real(self.data1D))
            self.ymaxlim=max(np.real(self.data1D))
        elif self.plotType==1:
            self.yminlim=min(np.imag(self.data1D))
            self.ymaxlim=max(np.imag(self.data1D))
        elif self.plotType==2:
            self.yminlim=min(min(np.real(self.data1D)),min(np.imag(self.data1D)))
            self.ymaxlim=max(max(np.real(self.data1D)),max(np.imag(self.data1D)))
        elif self.plotType==3:
            self.yminlim=min(np.abs(self.data1D))
            self.ymaxlim=max(np.abs(self.data1D))
        else:
            self.yminlim=-1
            self.ymaxlim=1
        if self.spec:
            self.xminlim=-self.sw/2.0
            self.xmaxlim=self.sw/2.0
        else:
            self.xminlim=0
            self.xmaxlim=(len(self.data1D)-1)/(2.0*self.sw)
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None,old=False): #display the 1D data
        if tmpdata is None:
            tmpdata=self.data1D
        a=self.fig.gca()
        a.cla()
        if not(self.spec):
            x=np.arange(len(tmpdata))/(2.0*self.sw)
        else:
            x=np.fft.fftshift(np.fft.fftfreq(len(tmpdata),1.0/self.sw))
        if old:
            if (self.plotType==0):
                a.plot(x,np.real(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==1):
                a.plot(x,np.imag(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==2):
                a.plot(x,np.real(self.data1D),c='k',alpha=0.2)
            elif(self.plotType==3):
                a.plot(x,np.abs(self.data1D),c='k',alpha=0.2)
        if (extraX is not None):
            for num in range(len(extraX)):
                a.plot(extraX[num],extraY[num],c=extraColor[num])
        if (self.plotType==0):
            self.line = a.plot(x,np.real(tmpdata),c='b')
        elif(self.plotType==1):
            self.line = a.plot(x,np.imag(tmpdata),c='b')
        elif(self.plotType==2):
            a.plot(x,np.imag(tmpdata),c='r')
            self.line = a.plot(x,np.real(tmpdata),c='b')
        elif(self.plotType==3):
            self.line = a.plot(x,np.abs(tmpdata),c='b')
        a.set_title("TD"+str(self.axes+1))
        a.set_xlabel('X axis label')
        a.set_ylabel('Y label')
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        a.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        a.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        #self.fig.set_tight_layout(True)
        self.canvas.draw()

    ################
    # mouse events #
    ################

    def peakPickReset(self):
        self.rect=[None,None,None,None]
        self.peakPick=False
        self.peakPickFunc = None

    def scroll(self,event):
        a=self.fig.gca()
        if self.rightMouse:
            middle = (self.xmaxlim+self.xminlim)/2.0
            width = self.xmaxlim-self.xminlim
            width = width*0.9**event.step
            self.xmaxlim = middle+width/2.0
            self.xminlim = middle-width/2.0
            a.set_xlim(self.xminlim,self.xmaxlim)
        else:
            middle = (self.ymaxlim+self.yminlim)/2.0
            width = self.ymaxlim-self.yminlim
            width = width*0.9**event.step
            self.ymaxlim = middle+width/2.0
            self.yminlim = middle-width/2.0
            a.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

    def buttonPress(self,event):
        if event.button == 1 and not self.peakPick:
            self.leftMouse = True
            self.zoomX1 = event.xdata
            self.zoomY1 = event.ydata
        elif (event.button == 3) and event.dblclick:
            self.plotReset()
        elif event.button == 3:
            self.rightMouse = True
            self.panX = event.xdata
            self.panY = event.ydata

    def buttonRelease(self,event):
        a=self.fig.gca()
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0]=None
                    self.peakPick = False
                    xdata = self.line[0].get_xdata()
                    ydata = self.line[0].get_ydata()
                    idx =np.argmin(np.abs(xdata-event.xdata))
                    if self.peakPickFunc is not None:
                        self.peakPickFunc((idx,xdata[idx],ydata[idx]))
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
                    a.set_xlim(self.xminlim,self.xmaxlim)
                    a.set_ylim(self.yminlim,self.ymaxlim)
                self.zoomX1=None
                self.zoomX2=None #WF: should also be cleared, memory of old zoom
                self.zoomY1=None
                self.zoomY2=None #WF: should also be cleared, memory of old zoom
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()

    def pan(self,event):
        if self.rightMouse and self.panX is not None and self.panY is not None:
            a=self.fig.gca()
            inv = a.transData.inverted()
            point = inv.transform((event.x,event.y))
            diffx = point[0]-self.panX
            diffy = point[1]-self.panY
            self.xmaxlim = self.xmaxlim-diffx
            self.xminlim = self.xminlim-diffx
            self.ymaxlim = self.ymaxlim-diffy
            self.yminlim = self.yminlim-diffy
            a.set_xlim(self.xminlim,self.xmaxlim)
            a.set_ylim(self.yminlim,self.ymaxlim)
            self.canvas.draw()
        elif self.peakPick:
            a=self.fig.gca()
            if self.rect[0] is not None:
                self.rect[0].remove()
                self.rect[0]=None
            if event.xdata is not None:
                self.rect[0]=a.axvline(event.xdata,c='k',linestyle='--')
            self.canvas.draw()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            a=self.fig.gca()
            inv = a.transData.inverted()
            point = inv.transform((event.x,event.y))
            self.zoomX2 =  point[0]
            self.zoomY2 = point[1]
            if self.rect[0] is not None:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                if self.rect[1] is not None:
                    self.rect[1].remove()
                if self.rect[2] is not None:
                    self.rect[2].remove()
                if self.rect[3] is not None:
                    self.rect[3].remove()
                self.rect=[None,None,None,None]
            self.rect[0],=a.plot([self.zoomX1,self.zoomX2],[self.zoomY2,self.zoomY2],'k',clip_on=False)
            self.rect[1],=a.plot([self.zoomX1,self.zoomX2],[self.zoomY1,self.zoomY1],'k',clip_on=False)
            self.rect[2],=a.plot([self.zoomX1,self.zoomX1],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.rect[3],=a.plot([self.zoomX2,self.zoomX2],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.canvas.draw()
            
