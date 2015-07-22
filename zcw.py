import numpy as np

############# Fibonacci number generator ##################     
def fib(n):
    if n > 50:
        print('number of angles becomes too large')
    start = np.array([[1,1],[1,0]],dtype='int64')
    temp = start[:]
    for i in xrange(n):
        temp = np.dot(start,temp)
    return temp[0,0],temp[0,1],temp[1,1]
                
############# ZCW Angles ##########
def zcw_angles(m,symm=0):
    samples, fib_1, fib_2 = fib(m)
    js = np.arange(samples, dtype='Float64')/samples
    if symm ==0:
        #full
        c = (1.,2.,1.)
    elif symm==1:
        #hemi
        c = (-1.,1.,1.)
    elif symm==2:
        #oct
        c  = (-1.,1.,4.)
    j_samples = fib_2*js
    phi = 2*np.pi/c[2] * np.mod(j_samples,1.0)
    theta = np.arccos(c[0]*(c[1]*np.mod(js,1.0) - 1))
    weight = np.ones(samples)/samples
    return phi,theta,weight
