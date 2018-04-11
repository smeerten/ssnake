import numpy as np
import scipy.integrate as SI
import scipy.special as SP
import math
import multiprocessing
try: #If numba exists, use jit, otherwise make a mock decorator
    from numba import jit
except:
    def jit(func):
        return func


@jit
def gammaFunc(gamma,a11pa15,a51pre,a55part,preVal):
    cos2G = math.cos(2 * gamma)
    sin2G = math.sin(2 * gamma)
    a51 = a51pre * cos2G
    a55 = a55part * sin2G
    integral = math.exp( a11pa15 + a51 + a55)
    amp = preVal * integral
    return amp

def alphaFunc(alpha, preCalc,preVal,a11,a15pre,a55pre,a51prepre1,a51prepre2):
    cos2A = math.cos(2 * alpha)
    sin2A = math.sin(2 * alpha)
    a11pa15 = a11 - a15pre * cos2A
    a55part = sin2A * a55pre
    a51pre = a51prepre1 + a51prepre2 * cos2A
    
    Amp = SI.quad(gammaFunc,0, np.pi,args = (a11pa15,a51pre,a55part,preVal),epsrel=5,epsabs=5)
    return Amp[0]

def betaFunc(beta, eta, eta0, preCalc):
    cosB = math.cos(beta)
    cosBS = cosB**2
    sinB = math.sin(beta)
    sinBS = sinB**2
    preVal = preCalc[0] * sinB
    a11 =   - 0.5 * (3 * cosBS - 1) * preCalc[1] * preCalc[2] + preCalc[3]
    a15pre = eta0 * preCalc[1] / 2 * sinBS * preCalc[2]
    a55pre = -cosB * preCalc[4] * preCalc[2]
    a51prepre1 = preCalc[1] / 2 * sinBS * eta * preCalc[2]
    a51prepre2 = 0.5 * (1 + cosBS) * preCalc[4] * preCalc[2]

    Amp = SI.quad(alphaFunc,0, np.pi/2,args = (preCalc,preVal,a11,a15pre,a55pre,a51prepre1,a51prepre2),epsrel=5,epsabs=5)
    return Amp[0]


def Czjzek(inp):
    wq, eta, wq0, eta0, sigma, d = inp
    N1 = wq ** (d - 1) / (sigma ** d) * eta * (1 - eta**2 / 9)
    afterfact = -1.0 / (2 * sigma ** 2)
    preCalc =[N1,np.sqrt(3),  2.0/np.sqrt(3) * wq * wq0 * afterfact, (wq0**2 * (1 + eta0**2 / 3) + wq**2 * (1 + eta**2/3)) * afterfact]
    eta0deveta = eta0/preCalc[1]*eta #eta0 divided by sqrt3 multiply by eta
    preCalc.append(eta0deveta)
    Amp = SI.quad(betaFunc,0, np.pi,args = (eta,eta0,preCalc),epsrel=1.8,epsabs=3)
    return Amp[0]

@jit
def tFunc(t,pre,pre2, fact, eta):
    expo = math.exp(fact * (3*t**2-1) + pre2)
    z = eta * abs(fact) * (1-t**2)
    bessel = SP.iv(0,z)
    return pre * expo * bessel
    

def CzjzekNoEta0(inp):
    wq, eta, wq0, sigma, d = inp
    pre = wq**(d-1)/sigma**d * eta * (1-eta**2/9.0) 
    pre2 = -(wq0**2 + wq**2 * (1+eta**2/3.0))/(2 * sigma**2)
    fact = wq*wq0/(2*sigma**2)
    Amp = SI.quad(tFunc,0, 1,args = (pre,pre2,fact,eta))[0]
    return Amp

def getInts(sigma, d, eta0, wq0, wq, eta):
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    if eta0 != 0.0: #If eta not 0, use hard routine
        fit = pool.map_async(Czjzek, [(wq[i],eta[i],wq0,eta0,sigma,d) for i in range(len(wq))])
    else:
        fit = pool.map_async(CzjzekNoEta0, [(wq[i],eta[i],wq0,sigma,d) for i in range(len(wq))])
    pool.close()
    pool.join()
    Amp = np.array(fit.get())
    if np.sum(Amp) == 0.0: #Protect against divide by zero
        Amp = 0 * Amp
    else:
        Amp = Amp / np.sum(Amp)
    return Amp
