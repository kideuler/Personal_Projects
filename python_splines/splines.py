import numpy as np
import math
import scipy as sp
from dataclasses import dataclass

@dataclass
class splinestruct:
    nv: int
    degree: int
    x: np.array
    y: np.array

def spline2(xs):
    nv = np.size(xs,0)
    x = xs[:,0]
    y = xs[:,1]

    bx = np.zeros((nv,1),dtype=float)
    by = np.zeros((nv,1),dtype=float)
    A = np.zeros((nv,nv),dtype=float)

    bx[-1] = 2*(x[0]-x[-1])
    by[-1] = 2*(y[0]-y[-1])

    for i in range(nv-1):
        bx[i] = 2*(x[i+1]-x[i])
        by[i] = 2*(y[i+1]-y[i])

    A[-1,0] = 1
    A[-1,-1] = 1
    for i in range(nv-1):
        A[i,i] = 1
        A[i,i+1] = 1

    #Q,R = np.linalg.qr(A,mode='reduced')
    Dx,res,rank,s = np.linalg.lstsq(A,bx)
    Dy,res,rank,s = np.linalg.lstsq(A,by)

    #A_inv = np.linalg.inv(A)
    #Dx = np.dot(A_inv,bx)
    #Dy = np.dot(A_inv,by)

    cx = np.zeros((nv,1),dtype=float)
    cy = np.zeros((nv,1),dtype=float)
    
    cx[-1] = (Dx[0]-Dx[-1])/2
    cy[-1] = (Dy[0]-Dy[-1])/2
    for i in range(nv-1):
        cx[i] = (Dx[i+1]-Dx[i])/2
        cy[i] = (Dy[i+1]-Dy[i])/2

    x = np.array(x)[np.newaxis].T
    y = np.array(y)[np.newaxis].T

    Mx = np.stack([x,Dx,cx],axis=1)
    Mx = Mx.reshape((nv,3))
    My = np.stack([y,Dy,cy],axis=1)
    My = My.reshape((nv,3))
    spline = splinestruct(nv,2,Mx,My)
    return(spline)






def spline_init(xs, degree):
    if degree == 1:
        spline = spline1(xs)
        return(spline)
    if degree == 2:
        spline = spline2(xs)
        return(spline)
    if degree == 3:
        spline = spline3(xs)
        return(spline)



def spline_var(t,spline):
    if t<0:
        t=1-T
    if t>1:
        t = t%1
    
    n = math.ceil(spline.nv*t)-1
    if n == -1:
        n=0

    x_j = (t-(n/(spline.nv)))*spline.nv
    vec = np.zeros((spline.degree+1),dtype=float)

    for i in range(spline.degree+1):
        vec[i] = x_j**(i)
    
    xs = [np.dot(spline.x[n,:],vec), np.dot(spline.y[n,:],vec)]

    
    return(xs)