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

    Dx,res,rank,s = np.linalg.lstsq(A,bx)
    Dy,res,rank,s = np.linalg.lstsq(A,by)

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


def spline3(xs):
    nv = np.size(xs,0)
    x = xs[:,0]
    y = xs[:,1]

    bx = np.zeros((nv,1),dtype=float)
    by = np.zeros((nv,1),dtype=float)
    A = np.zeros((nv,nv),dtype=float)

    bx[0] = 3*(x[1]-x[-1])
    by[0] = 2*(y[1]-y[-1])
    bx[-1] = 3*(x[0]-x[-2])
    by[-1] = 3*(y[0]-y[-2])

    for i in range(1,nv-1):
        bx[i] = 2*(x[i+1]-x[i-1])
        by[i] = 2*(y[i+1]-y[i-1])

    A[0,0] = 4
    A[0,1] = 1
    A[0,-1] = 1
    A[-1,1] = 1
    A[-1,-2] = 1
    A[-1,-1] = 4
    for i in range(1,nv-1):
        A[i,i-1] = 1
        A[i,i] = 4
        A[i,i+1] = 1

    A_inv = np.linalg.inv(A)
    Dx = np.dot(A_inv,bx)
    Dy = np.dot(A_inv,by)

    cx = np.zeros((nv,1),dtype=float)
    cy = np.zeros((nv,1),dtype=float)
    dx = np.zeros((nv,1),dtype=float)
    dy = np.zeros((nv,1),dtype=float)

    cx[-1] = 3*(x[0]-x[-1]) - 2*Dx[-1] - Dx[0]
    cy[-1] = 3*(y[0]-y[-1]) - 2*Dy[-1] - Dy[0]
    dx[-1] = 2*(x[-1]-x[0]) + Dx[-1] + Dx[0]
    dy[-1] = 2*(y[-1]-y[0]) + Dy[-1] + Dy[0]

    for i in range(nv-1):
        cx[i] = 3*(x[i+1]-x[i]) - 2*Dx[i] - Dx[i+1]
        cy[i] = 3*(y[i+1]-y[i]) - 2*Dy[i] - Dy[i+1]
        dx[i] = 2*(x[i]-x[i+1]) + Dx[i] + Dx[i+1]
        dy[i] = 2*(y[i]-y[i+1]) + Dy[i] + Dy[i+1]
    
    x = np.array(x)[np.newaxis].T
    y = np.array(y)[np.newaxis].T

    Mx = np.stack([x,Dx,cx,dx],axis=1)
    Mx = Mx.reshape((nv,4))
    My = np.stack([y,Dy,cy,dy],axis=1)
    My = My.reshape((nv,4))
    spline = splinestruct(nv,3,Mx,My)
    return(spline)

def spline_n(xs,degree):
    global columns
    nv = np.size(xs,0)
    columns = (degree+1)*nv
    x = xs[:,0]
    y = xs[:,1]

    bx = np.zeros((columns,1),dtype=float)
    by = np.zeros((columns,1),dtype=float)
    A = np.zeros((columns,columns),dtype=float)

    t = np.array([0,1])
    V = np.zeros((2*degree,degree+1))
    k=0
    for i in range(degree):
        for j in range(i,degree+1):
            V[k:k+2,j] = (math.factorial(j)/math.factorial(j-i))*t**(j-i)
        k+=2   

    k=0
    for i in range(nv):
        bx[2*(k)] = x[k]
        bx[2*(k)+1] = x[(k+1)%nv]
        by[2*(k)] = y[k]
        by[2*(k)+1] = y[(k+1)%nv]
        k+=1
    

    k=0
    vec = V[0:2,:]
    for i in range(nv):
        A[2*k,((degree+1)*k) : ((degree+1)*(k+1))] = vec[0,:]
        A[2*k+1,((degree+1)*k) : ((degree+1)*(k+1))] = vec[1,:]
        k+=1

    for deg in range(1,degree):
        vec = V[(2*deg):(2*deg+2),:]
        p = 0
        k = nv + deg*nv
        for i in range(nv):
            A[k,((degree+1)*p):((degree+1)*(p+1))] = vec[1,:]
            A[k, ((degree+1)*(p+1))%columns : ((degree+1)*(p+2)-1)%(columns)+1] = -vec[0,:]
            k+=1
            p+=1

    if degree%2 == 1:
        A_inv = np.linalg.inv(A)
        Dx = np.dot(A_inv,bx)
        Dy = np.dot(A_inv,by)
    elif degree%2 == 0:
        Dx,res,rank,s = np.linalg.lstsq(A,bx)
        Dy,res,rank,s = np.linalg.lstsq(A,by)

    Mx = np.zeros((nv,degree+1))
    My = np.zeros((nv,degree+1))

    for i in range(nv):
        Mx[i,:] = Dx[(degree+1)*i : (degree+1)*(i+1),0]
        My[i,:] = Dy[(degree+1)*i : (degree+1)*(i+1),0]

    spline = splinestruct(nv,degree,Mx,My)
    return(spline)

def spline_init(xs, degree):
    if degree == 2:
        spline = spline2(xs)
        return(spline)
    elif degree == 3:
        spline = spline3(xs)
        return(spline)
    else:
        spline = spline_n(xs,degree)
        return(spline)




def spline_var(t,spline,nargout=1):
    if t<0:
        t=1-t
    if t>1:
        t = t%1
    
    n = math.ceil(spline.nv*t)
    if n == 0:
        n=1

    x_j = (t-(n-1)/(spline.nv))*spline.nv
    vec = np.zeros((spline.degree+1),dtype=float)

    for i in range(spline.degree+1):
        vec[i] = x_j**(i)
    
    xs = [np.dot(spline.x[n-1,:],vec), np.dot(spline.y[n-1,:],vec)]
    if nargout == 1:
        return(xs)

    vec = np.zeros((spline.degree+1),dtype=float)
    for i in range(1,spline.degree+1):
        vec[i] = (i)*x_j**(i-1) 
    grad = [np.dot(spline.x[n-1,:],vec), np.dot(spline.y[n-1,:],vec)]
    if nargout == 2:
        return(xs,grad)

    vec = np.zeros((spline.degree+1),dtype=float)
    for i in range(2,spline.degree+1):
        vec[i] = i*(i-1)*x_j**(i-2)
    hess = [np.dot(spline.x[n-1,:],vec), np.dot(spline.y[n-1,:],vec)]
    if nargout == 3:
        return(xs,grad,hess)
    

def spline_proj(xs, spline):
    x_0 = xs[0]
    y_0 = xs[1]

    avex = np.sum(spline.x[:,0])/spline.nv
    avey = np.sum(spline.y[:,0])/spline.nv
    u = (np.arctan2(spline.y[0,0]-avey,spline.x[0,0]-avex)+np.pi)/(2*np.pi)
    t = (np.arctan2(y_0-avey,x_0-avex)+np.pi)/(2*np.pi) - u
    t_init = t

    iter = 0
    f = 1
    grad = np.array([0,0])
    curve = np.array([0,0])
    while iter < 1000:
        normal = np.array([grad[1], -grad[0]])/np.linalg.norm([grad[1],-grad[0]])
        vec = np.array([x_0-curve[0],y_0-curve[1]])
        vec = vec/np.linalg.norm(vec)

        dir = np.dot(normal,vec)
        if abs(f) < 0.01 and abs(dir-1)<0.01:
            break
        
        iter+=1
        curve,grad = spline_var(t,spline,2)
        f = grad[0]*(curve[0]-x_0) + grad[1]*(curve[1]-y_0)
        
        xspp = spline_var((t+0.001)%1,spline)
        xspn = spline_var((t-0.001)%1,spline)
        if np.linalg.norm(xs-xspp) < np.linalg.norm(xs-curve):
            t = (t+0.001)%1
        elif np.linalg.norm(xs-xspn) < np.linalg.norm(xs-curve):
            t = (t-0.001)%1

    if iter == 2000:
        print('max iters reached \n')
    xsp = spline_var(t,spline)
    return(xsp)