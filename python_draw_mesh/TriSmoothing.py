import numpy as np
import math
import triangle as tr
import matplotlib.pyplot as plt
from drawcurve import *

def tri_energy(xs):

    # compute edge lengths
    e12 = np.subtract(xs[1,:],xs[0,:])
    sql12 = np.dot(e12,e12)
    e23 = np.subtract(xs[2,:],xs[1,:])
    sql23 = np.dot(e23,e23)
    e31 = np.subtract(xs[0,:],xs[2,:])
    sql31 = np.dot(e31,e31)

    area2 = e12[0]*e23[1]-e12[1]*e23[0]
    area = 0.5*area2
    e12_orth = [-e12[1], e12[0]]
    e23_orth = [-e23[1], e23[0]]
    e31_orth = [-e31[1], e31[0]]
    
    cts_a = np.zeros((3,1))
    cts_a[:] = (1/math.sqrt(3))/area

    energy = cts_a[2]*sql12 + cts_a[0]*sql23 + cts_a[1]*sql31

    # Gradient
    t = -2*cts_a[2] * e12
    grads1  = t
    grads2 = -t
    t = -2*cts_a[0] * e23
    grads2 = np.add(grads2,t)
    grads3 = -t
    t = -2*cts_a[1] * e31
    grads3 = np.add(grads3,t)
    grads1 = np.add(grads1,-t) 

    energy_a = energy/area2
    grads1 = np.subtract(grads1,energy_a*e23_orth)
    grads2 = np.subtract(grads2,energy_a*e31_orth)
    grads3 = np.subtract(grads3,energy_a*e12_orth)

    grad = np.array([grads1, grads2, grads3])
    grad = np.transpose(grad)
    
    # Hessian
    d = 2*(cts_a[2]+cts_a[1])
    C = d*np.eye(2)

    A = np.outer(np.transpose(grads1/area2),e23_orth)
    hess1 = C - (A + np.transpose(A))

    A = np.outer(np.transpose(grads2/area2),e31_orth)
    hess2 = C - (A + np.transpose(A))

    A = np.outer(np.transpose(grads3/area2),e12_orth)
    hess3 = C - (A + np.transpose(A))

    hess = np.array([hess1, hess2, hess3])
    hess = np.transpose(hess)
    return(grad, hess)


def Mesh_smoothing(Tris, iters, app, canvas, h):
    is_bnd = Tris['vertex_markers']
    nelems = np.size(Tris['triangles'],axis=0)
    nv = np.size(Tris['vertices'],axis=0)

    conn_elem = -np.ones((nv,10),dtype=int)
    loc_node = np.zeros((nv,10),dtype=int)
    nnz = np.zeros(nv,dtype=int)

    for i in range(nelems):
        for j in range(3):
            conn_elem[Tris['triangles'][i,j], nnz[Tris['triangles'][i,j]]] = i
            loc_node[Tris['triangles'][i,j], nnz[Tris['triangles'][i,j]]] = j
            nnz[Tris['triangles'][i,j]] = nnz[Tris['triangles'][i,j]]+1

    for p in range(iters):
        for n in range(nv):
            if is_bnd[n] == True:
                continue
            
            Grad = np.zeros((1,2))
            Hess = np.zeros((2,2))
            for E in range(nnz[n]):
                grad, hess = tri_energy(Tris['vertices'][Tris['triangles'][conn_elem[n,E],:],:])
                Grad = np.add(Grad, grad[:,loc_node[n,E]])
                Hess = np.add(Hess, hess[:,:,loc_node[n,E]])
            
            Grad = Grad/nnz[n]
            Hess = Hess/nnz[n]
            is_nan = False
            Hess = np.transpose(Hess)
            Hess_inv = np.linalg.inv(Hess)
            if np.isnan(Grad).any() or np.isnan(Hess_inv).any():
                is_nan = True

            if is_nan == False:
                r = np.dot(Hess_inv, np.transpose(Grad))
            else:
                r = np.transpose([0,0])

            if np.linalg.norm(r) > 1.5*h:
                r = np.transpose([0,0])

            Tris['vertices'][n,:] = np.add(Tris['vertices'][n,:], np.transpose(-r))
        Draw_mesh(Tris,app,canvas)
    return(app,canvas)

