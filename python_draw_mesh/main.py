from drawcurve import *
from TriSmoothing import *
import numpy as np
import matplotlib.pyplot as plt
import triangle as tr
import math
from tkinter import *
from time import *
import sys

def main_subroutine(app, in_):
    # get intial points of curve
    xs, app, canvas = draw_on_canvas(app)

    # delete every other node twice and form initial segments
    for i in range(2):
        xs = np.delete(xs, list(range(0, xs.shape[0], 2)), axis=0)
    nv = np.size(xs,0)
    s1 = np.arange(nv)
    segs = np.stack([s1,s1+1],axis=1) % nv
    Draw_curve(xs,segs, app, canvas)



    # loop through each segment and delete nodes too close to one another
    Keep  = np.ones(nv, dtype=bool)
    h = 0.0
    for i in range(nv):
        h = h + math.sqrt((xs[segs[i,1],0]-xs[segs[i,0],0])**2 + (xs[segs[i,1],1]-xs[segs[i,0],1])**2)
    h = h/nv
    for i in range(nv):
        norm = math.sqrt((xs[segs[i,1],0]-xs[segs[i,0],0])**2 + (xs[segs[i,1],1]-xs[segs[i,0],1])**2)
        if norm < 1.5*h:
            xs[segs[i,1],0] = xs[segs[i,0],0]
            xs[segs[i,1],1] = xs[segs[i,0],1]
            Keep[segs[i,1]] = False
    xs = xs[Keep,:]

    # add interior/exterior points
    nv = np.size(xs,0)
    s1 = np.arange(nv)
    segs = np.stack([s1,s1+1],axis=1) % nv
    canvas.delete("all")
    Draw_curve(xs,segs, app, canvas)

    h = 0.0
    for i in range(nv):
        h = h + math.sqrt((xs[segs[i,1],0]-xs[segs[i,0],0])**2 + (xs[segs[i,1],1]-xs[segs[i,0],1])**2)
    h = h/nv
    ii = math.ceil(1/h)
    k = 0
    x = np.linspace(0.0,1.0,num=ii)
    ps = np.zeros((ii**2, 2))
    seg_iter = 0
    box_segs = np.zeros((4*(ii-1), 2))

    for i in range(ii-1):
        ps[k,:] = [x[i],0]
        k+=1
    for i in range(ii-1):
        ps[k,:] = [1,x[i]]
        k+=1
    for i in range(ii-1):
        ps[k,:] = [x[ii-1-i], 1]
        k+=1
    for i in range(ii-1):
        ps[k,:] = [0,x[ii-1-i]]
        k+=1

    for i in range(1,ii-1):
        for j in range(1,ii-1):
            ps[k,:] = [x[i], x[j]] 
            k = k+1 
    mask = find_winding_number(ps,xs,segs)

    if in_ == "1":
        ps = ps[mask,:]
        center = [0.0,0.0]
    elif in_ == "0":
        ps = ps[~mask,:]
        nnv = 4*(ii-1) + nv
        s1 = np.arange(4*ii-4)
        box_segs = np.stack([s1,s1+1],axis=1) % (4*ii-4) + nv
        segs = np.append(segs, box_segs,axis=0)

        center = [0.0,0.0]
        for i in range(nv):
            center += xs[i,:]
        center = center/nv

    xs = np.append(xs,ps,axis=0)

    # perform constrained delaunay triangulation
    xs = np.transpose(xs)
    A = dict(vertices=np.transpose(xs), segments=segs)
    tri = tr.triangulate(A, 'p')
    Draw_mesh(tri,app,canvas)

    # smoothing
    app,canvas = Mesh_smoothing(tri,5,app,canvas,h)
    return(app,canvas)

global app, canvas
def start_main():
    global app, canvas, in_
    canvas.destroy()
    app, canvas = main_subroutine(app, in_)

in_ = sys.argv[1]
app  = Tk()
app.geometry("500x500")
app.title("Draw to mesh")
canvas = Canvas(app)
B = Button(app,text = "Reset",command = start_main)
B.pack()
app.mainloop()