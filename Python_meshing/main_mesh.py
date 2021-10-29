from drawpoints import *
from tkinter import *
from splines import *
import numpy as np
import matplotlib.pyplot as plt
from TriSmoothing import *
import triangle as tr
import sys
import time

def main_subroutine(app, degree, nv, Point):
    if Point == 1:
        xs,app,canvas = draw_on_canvas(app)
    if Point == 2:
        xs, app, canvas = draw_curve_on_canvas(app)
        xs = xs[0::2,:]

    nn = np.size(xs,0)
    canvas.delete("all")
    spline = spline_init(xs, degree)

    #draw curve and normals of the spline
    n = 100
    ps = np.zeros((n,2))
    grad = np.zeros((n,2))
    for i in range(n):
        ps[i,:], grad[i,:] = spline_var(float(i/n),spline,2)
    
    Draw_curve(ps, app, canvas)
    Draw_normal(ps,grad,app,canvas)

    # create boundary nodes based on the inpute nv
    bnd_xs = np.zeros((nv,2))
    for i in range(nv):
        bnd_xs[i,:] = spline_var(float(i/nv),spline)
    Draw_points(bnd_xs,app,canvas)
    
    # create 1D segment elements
    s1 = np.arange(nv)
    segs = np.stack([s1,s1+1],axis=1) % nv

    # find average edge length
    h = 0.0
    for i in range(nv):
        h = h + math.sqrt((bnd_xs[segs[i,1],0]-bnd_xs[segs[i,0],0])**2 + (bnd_xs[segs[i,1],1]-bnd_xs[segs[i,0],1])**2)
    h = h/nv

    # create a bunch of points on grid and determine which ones are inside the curve
    ii = math.ceil(1/h)
    k = 0
    x = np.linspace(0.0,1.0,num=ii)
    ps = np.zeros((ii**2, 2))
    for i in range(1,ii-1):
        for j in range(1,ii-1):
            ps[k,:] = [x[i], x[j]] 
            k = k+1 
    mask = find_winding_number(ps,bnd_xs,segs,h)
    ps = ps[mask,:]
    xs = np.append(bnd_xs,ps,axis=0)

    # perform constrained delaunay triangulation
    xs = np.transpose(xs)
    A = dict(vertices=np.transpose(xs), segments=segs)
    tri = tr.triangulate(A, 'p')
    time.sleep(1)
    Draw_mesh(tri,app,canvas)

    # smoothing
    app,canvas = Mesh_smoothing(tri,spline,10,app,canvas,h)
    return(app,canvas)



global app, canvas
def start_main():
    global app, canvas, degree, n
    canvas.destroy()
    app, canvas = main_subroutine(app, degree, n, Point)

degree = int(sys.argv[1])
n = int(sys.argv[2])
Point = int(sys.argv[3])
app  = Tk()
app.geometry("500x500")
app.title("Draw to mesh")
canvas = Canvas(app)
B = Button(app,text = "Reset",command = start_main)
B.pack()
app.mainloop()