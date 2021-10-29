import numpy as np
import math
import triangle as tr
import matplotlib.pyplot as plt
from drawpoints import *
from splines import *
import time

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


def Mesh_smoothing(Tris, spline, iters, app, canvas, h):

    # initial angles
    mina,maxa = find_angles(Tris)
    S = "Min angle:"+ "{:10.4f}".format(mina) +"   Max angle:"+ "{:10.4f}".format(maxa)
    canvas.create_text(250,10,text = S)
    app.update()

    is_bnd = Tris['vertex_markers']
    nelems = np.size(Tris['triangles'],axis=0)
    nv = np.size(Tris['vertices'],axis=0)

    conn_elem = -np.ones((nv,10),dtype=int)
    loc_node = np.zeros((nv,10),dtype=int)
    nnz = np.zeros(nv,dtype=int)
    nbnd = 0
    for i in range(nv):
        if is_bnd[i] == True:
            nbnd+=1

    t_curve = -np.ones((nv,1))
    for i in range(nv):
        if is_bnd[i] == True:
            t_curve[i] = float(i/nbnd)

    for i in range(nelems):
        for j in range(3):
            conn_elem[Tris['triangles'][i,j], nnz[Tris['triangles'][i,j]]] = i
            loc_node[Tris['triangles'][i,j], nnz[Tris['triangles'][i,j]]] = j
            nnz[Tris['triangles'][i,j]] = nnz[Tris['triangles'][i,j]]+1

    p=0
    while p < iters and (mina<30 or maxa>120):
        for n in range(nv):
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
            if is_bnd[n] == False:
                if np.isnan(Grad).any() or np.isnan(Hess_inv).any():
                    is_nan = True
                if is_nan == False:
                    r = np.dot(Hess_inv, np.transpose(Grad))
                    flag = Check_Neighbor_Jacobians(n,conn_elem[n,0:nnz[n]], np.transpose(r)[0,:], Tris)
                    if flag == False:
                        Tris['vertices'][n,:] -= np.transpose(r)[0,:]

            elif is_bnd[n] == True: # if veetex on boundary

                # evaluate Tangent line and hessian on spline curve
                t_0 = t_curve[n]
                curve, T, H = spline_var(t_0,spline,3)
                T = T/np.linalg.norm(T)
                H = H/np.linalg.norm(H)

                # reshape
                T = np.reshape(T,(2,1))
                H = np.reshape(H,(2,1))

                # iterations step of newton minimization
                rhs = np.dot(Grad,T)
                lhs = np.dot(np.transpose(T),np.dot(Hess,T)) + np.dot(Grad,H)
                res = rhs/lhs
                
                # attempt to find t that matches x,y values on curve
                t,found = Rootfinding(Tris['vertices'][n,:],spline,t_0)

                if found == True and found==False: # broken right now
                    # if found move vertex
                    t = t - res[0,0]
                    xs_new = spline_var(t,spline)
                    r = Tris['vertices'][n,:]-xs_new
                    t_curve[n] = t

                else:
                    # if root is not found move along tangent
                    #T = np.reshape([0,0],(2,1))
                    r = res[0,0]*np.transpose(T)[0,:]

                # check Jacobians and if good move vertex
                flag = Check_Neighbor_Jacobians(n,conn_elem[n,0:nnz[n]], r, Tris)
                if (flag == False) and (np.linalg.norm(r) < 0.5*h):
                    Tris['vertices'][n,:] -=r

        # draw mesh
        Draw_mesh(Tris,app,canvas)
        # update min/max angles
        mina,maxa = find_angles(Tris)
        S = "Min angle: "+ "{:10.4f}".format(mina) +" Max angle: "+ "{:10.4f}".format(maxa)
        canvas.create_text(250,10,text = S)
        app.update()
        time.sleep(.3)
        p=p+1
    return(app,canvas)

# check whether moving a vertex would cause any neghboring elements to have negative jacobians
# this in result prevents any mesh folding
def Check_Neighbor_Jacobians(n, Conn_elem, r, elems):
    xs = elems['vertices'][n,:]
    elems['vertices'][n,:] -= r
    dphi = np.array([[-1,-1],[1,0],[0,1]])
    nelems = len(Conn_elem)
    for i in range(nelems):
        U = np.transpose(elems['vertices'][elems['triangles'][Conn_elem[i],:],:])
        Jac = np.linalg.det(np.dot(U,dphi))
        if Jac < 0:
            elems['vertices'][n,:] +=r
            return(True)

    elems['vertices'][n,:] +=r
    return(False)

# simple newtons iteration for finding parametric variable on spline  
def Rootfinding(xs,spline,t_0):
    res = abs(2-t_0)
    iter = 0
    t = t_0
    while res > 0.001 and iter < 1000:
        curve, grad = spline_var(t,spline,2)
        grad = grad/np.linalg.norm(grad)
        t_0 = t
        t = t - (curve[0]-xs[0])/(2*grad[0]) - (curve[1]-xs[1])/(2*grad[1])
        iter+=1
        res = abs(t-t_0)
        if res < 0.001:
            return(t,True)
    return(t,False)

def find_angles(Tris):
    min_angle = 180.0
    max_angle = 0.0

    nelems = np.size(Tris['triangles'],0)
    for i in range(nelems):
        xs = Tris['vertices'][Tris['triangles'][i,:],:]
        v1 = xs[1,:]-xs[0,:]
        v2 = xs[2,:]-xs[1,:]
        v3 = xs[0,:]-xs[2,:]

        v1 = v1/np.linalg.norm(v1)
        v2 = v2/np.linalg.norm(v2)
        v3 = v3/np.linalg.norm(v3)

        a1 = np.degrees(np.arccos(np.inner(v1,-v3)))
        a2 = np.degrees(np.arccos(np.inner(-v1,v2)))
        a3 = np.degrees(np.arccos(np.inner(-v2,v3)))

        min_angle = min(min_angle,a1,a2,a3)
        max_angle = max(max_angle,a1,a2,a3)

    return(min_angle,max_angle)