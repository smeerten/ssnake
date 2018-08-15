import numpy as np
import scipy.integrate as SI
import scipy.special as SP
import math
import multiprocessing
try: #If numba exists, use jit, otherwise make a mock decorator
    from numba import jit
except Exception:
    def jit(func):
        return func


@jit
def gammaFunc(gamma, a11pa15, a51pre, a55part, preVal):
    cos2G = math.cos(2 * gamma)
    sin2G = math.sin(2 * gamma)
    a51 = a51pre * cos2G
    a55 = a55part * sin2G
    integral = math.exp(a11pa15 + a51 + a55)
    amp = preVal * integral
    return amp

def alphaFunc(alpha, preCalc, preVal, a11, a15pre, a55pre, a51prepre1, a51prepre2):
    cos2A = math.cos(2 * alpha)
    sin2A = math.sin(2 * alpha)
    a11pa15 = a11 - a15pre * cos2A
    a55part = sin2A * a55pre
    a51pre = a51prepre1 + a51prepre2 * cos2A
    Amp = SI.quad(gammaFunc, 0, np.pi, args=(a11pa15, a51pre, a55part, preVal), epsrel=0.1, epsabs=0)
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
    Amp = SI.quad(alphaFunc, 0, np.pi/2, args=(preCalc,preVal,a11,a15pre,a55pre,a51prepre1,a51prepre2), epsrel=0.1, epsabs=0)
    return Amp[0]

def extendedCzjzek(inp):
    cq, eta, cq0, eta0, sigma, d = inp
    #Prevent overflow error (if this condition is true, the point is far outside the distribution, and is 0 anyway)
    if abs(cq**2*(1+eta**2/3) - cq0**2)/(2 * sigma**2) > 1000:
        return 0.0
    #Prevent overflow error on eta range
    if cq0 / sigma * abs(eta0 - eta) > 10:
        return 0.0
    N1 = cq ** (d - 1) / (sigma ** d) * eta * (1 - eta**2 / 9)
    afterfact = -1.0 / (2 * sigma ** 2)
    preCalc =[N1,np.sqrt(3),  2.0/np.sqrt(3) * cq * cq0 * afterfact, (cq0**2 * (1 + eta0**2 / 3) + cq**2 * (1 + eta**2/3)) * afterfact]
    eta0deveta = eta0 / preCalc[1] * eta # eta0 divided by sqrt3 multiply by eta
    preCalc.append(eta0deveta)
    Amp = SI.quad(betaFunc, 0, np.pi, args=(eta,eta0,preCalc), epsrel=0.1, epsabs=0)
    return Amp[0]

@jit
def tFunc(t, pre, pre2, fact, eta):
    expo = math.exp(fact * (3*t**2-1) + pre2)
    z = eta * abs(fact) * (1-t**2)
    bessel = SP.iv(0,z)
    return expo * bessel

def extendedCzjzekNoEta0(inp):
    cq, eta, cq0, sigma, d = inp
    #Prevent overflow error (if this condition is true, the point is far outside the distribution, and is 0 anyway)
    if abs(cq**2*(1+eta**2/3) - cq0**2)/(2 * sigma**2) > 1000: 
        return 0.0
    #Prevent overflow error on eta range
    if cq0 / sigma * eta > 10:
        return 0.0
    pre = cq**(d-1) / sigma**d * eta * (1 - eta**2 / 9.0) 
    pre2 = -(cq0**2 + cq**2 * (1 + eta**2 / 3.0)) / (2 * sigma**2)
    fact = cq * cq0 / (2*sigma**2)
    Amp = SI.quad(tFunc,0, 1, args=(pre,pre2,fact,eta))[0]
    return Amp * pre

def normalCzjzekFunc(cq, eta, sigma, d):
    return cq**(d - 1) * eta / (np.sqrt(2 * np.pi) * sigma**d) * (1 - eta**2 / 9.0) * np.exp(-cq**2 / (2.0 * sigma**2) * (1 + eta**2 / 3.0))

def czjzekIntensities(sigma, d, cq, eta, cq0=0, eta0=0):
    #Calculates an intensity distribution for a Czjzek library
    #Depending on et0 and cq0, either a regular Czjzek, or an extended Czjzek is used
    #cq: C_q grid (2D flattened)
    #eta: eta grid (2D flattened)
    if cq0 == 0.0 and eta0 == 0.0:
        if sigma == 0.0:  # protect against divide by zero
            czjzek = np.zeros_like(cq)
            czjzek[:, 0] = 1
        else:
            czjzek = normalCzjzekFunc(cq, eta, sigma, d)
    else:
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        if eta0 != 0.0:
            fit = pool.map_async(extendedCzjzek, [(cq[i],eta[i],cq0,eta0,sigma,d) for i in range(len(cq))])
        else:
            fit = pool.map_async(extendedCzjzekNoEta0, [(cq[i],eta[i],cq0,sigma,d) for i in range(len(cq))])
        pool.close()
        pool.join()
        czjzek = np.array(fit.get())
    if np.sum(czjzek) == 0.0: #Protect against divide by zero
        czjzek = 0 * czjzek
    else:
        czjzek = czjzek / np.sum(czjzek)
    return czjzek
