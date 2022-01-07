import numpy as np
import gmsh
from tkinter import *
from math import pi
import time
import pickle
import random
import tensorflow as tf

def Draw_mesh(elems,xs,bdy,app, canvas):
    canvas.delete("all")
    nv = np.size(xs, axis=0)
    nelems = np.size(elems, axis=0)
    nbnd = np.size(bdy,0)
    xs5 = 500*xs
    edge = np.array([[0,1],[1,2],[2,0]])
    for i in range(nelems):
        for j in range(3):
            canvas.create_line((xs5[elems[i,edge[j,0]],0], 500-xs5[elems[i,edge[j,0]],1], 
            xs5[elems[i,edge[j,1]],0], 500-xs5[elems[i,edge[j,1]],1]),
            fill='black', width=1)
    for n in range(nv):
        canvas.create_oval(xs5[n,0]-1,500-xs5[n,1]-1,xs5[n,0]+1,500-xs5[n,1]+1,fill='red')
    for n in range(nbnd):
        canvas.create_line((xs5[bdy[n,0],0], 500-xs5[bdy[n,0],1], 
            xs5[bdy[n,1],0], 500-xs5[bdy[n,1],1]),
            fill='blue', width=2)
    app.update()
    app.update_idletasks()
    return app,canvas

def gmsh_init(xs,h):
    # 2-D assertion
    assert(np.size(xs,1)==2)

    # initialize gmsh
    gmsh.initialize()

    # defining mesh geometry
    nv = np.size(xs,0)
    for n in range(nv):
        gmsh.model.geo.addPoint(xs[n,0], xs[n,1], 0, h, n+1)

    SP = list(range(1,nv+1))
    SP = np.append(SP,1)
    gmsh.model.geo.addSpline(SP,1)
    gmsh.model.geo.addCurveLoop([1],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()
    
    # generate mesh
    gmsh.model.mesh.generate(2)

    # obtaining adjacency matrix and nodal position data
    ps = gmsh.model.mesh.getNodes(-1,-1)
    ps = np.reshape(ps[1],(-1,3))
    ps = ps[:,[0,1]]
    elems = gmsh.model.mesh.getElements(2,-1)
    bdy = gmsh.model.mesh.getElements(1,-1)
    bdy = np.reshape(bdy[2],(-1,2))-1
    elems = np.reshape(elems[2],(-1,3))-1


    gmsh.finalize()
    return elems,ps,bdy


def create_data_lambda_functions():
    app  = Tk()
    app.geometry("500x500")
    app.title("Draw mesh")
    canvas = Canvas(app)
    canvas = Canvas(app, bg="white")
    canvas.pack(anchor='nw', fill='both', expand=1)

    # reading lambda functions from textfile
    lambdas = {}
    f = open('lambdas.txt','r')
    lines = f.read().splitlines()

    # copying lambdas in an array
    nl = len(lines)
    for i in range(nl):
        lambdas[i] = eval(lines[i])

    # creating meshes
    for i in range(nl):
        for n in range(100,301):
            t = np.linspace(0,1,n-1,False)
            ps = np.zeros((n-1,2))

            for k in range(n-1):
                ps[k,:] = lambdas[i](t[k])

            # creating mesh
            elems,xs,bdy = gmsh_init(ps,1/(n-1))

            # plot mesh
            app,canvas = Draw_mesh(elems,xs,bdy,app,canvas)

            # save training data as pickle files
            filename = 'training_data/'+str(i)+'_'+str(n)+'.pickle'
            with open(filename,'wb') as f:
                pickle.dump([elems,bdy,xs],f)

            print(filename)


def create_data_random_meshes():
    app  = Tk()
    app.geometry("500x500")
    app.title("Draw mesh")
    canvas = Canvas(app)
    canvas = Canvas(app, bg="white")
    canvas.pack(anchor='nw', fill='both', expand=1)


    for n in range(20,301):
        for i in range(10):
            ps = np.zeros((n-1,2))
            t = np.linspace(0,1,n-1,False)
            t = t.reshape((len(t),1))
            H = np.random.randint(5,np.floor(n/2))
            rho = np.random.rand(H,1)*np.logspace(-0.5,-2,H)
            phi = np.random.rand(H,1)*2*pi
            r = np.ones((np.size(t,0),1))

            for k in range(H):
                r = r + rho[k]*(np.sin(k*2*pi*t + phi[k])+1)

            for k in range(n-1):
                ps[k,:] = np.array([r[k,0]*np.cos(2*pi*t[k,0])/H +0.5,r[k,0]*np.sin(2*pi*t[k,0])/H+0.5])

            dx = max(ps[:,0])-min(ps[:,0])
            dy = max(ps[:,1])-min(ps[:,1])
            ps[:,0] = 0.7*(ps[:,0]-0.5)/dx + 0.5
            ps[:,1] = 0.7*(ps[:,1]-0.5)/dy + 0.5

            # creating mesh
            elems,xs,bdy = gmsh_init(ps,1/(n-1))
            

            # plot mesh
            app,canvas = Draw_mesh(elems,xs,bdy,app,canvas)
            #time.sleep(0.5)

            # save as a pickle file
            filename = 'training_data/rand_'+str(n)+'_'+str(i)+'.pickle'
            with open(filename,'wb') as f:
                pickle.dump([elems,bdy,xs],f)
            print(filename)

#create_data_lambda_functions()
#create_data_random_meshes()

print(tf.__version__)