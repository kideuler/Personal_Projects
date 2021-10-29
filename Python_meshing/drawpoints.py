from tkinter import *
import numpy as np
from splines import *

def draw_on_canvas(app):

    global xs 
    global flag
    flag = False
    xs = [[0.0, 0.0]]
    global stop
    stop = False
    def get_x_and_y(event):
        global lasx, lasy, xs
        lasx, lasy = event.x, event.y
        xs = np.append(xs, [[event.x, event.y]], axis=0)
        canvas.create_oval(event.x-5, event.y-5, event.x+5, event.y+5,fill='red')

    def end(event):
        global lasx, lasy, xs, flag
        xs = np.delete(xs,0,0)
        xs[:,1] = 500-xs[:,1]
        xs = xs/500
        flag = True
        

    
    canvas = Canvas(app, bg="white")
    canvas.pack(anchor='nw', fill='both', expand=1)
    canvas.bind("<Button-1>", get_x_and_y)
    canvas.bind("<Button-3>", end)
    
    while flag == False:
        app.update_idletasks()
        app.update()

    app.update_idletasks()
    app.update()

    return(xs, app, canvas)

# find if point is inside or outside curve
def find_winding_number(ps,xs,segs,h):
    nv = np.size(xs,0)
    nnv = np.size(ps,0)
    nsegs = np.size(segs,0)
    inn = np.zeros(nnv, dtype=bool)

    for n in range(nnv):
        wind = 0
        for i in range(nsegs):
            tempxs = np.array([[xs[segs[i,0],0] , xs[segs[i,0],1]], [xs[segs[i,1],0] , xs[segs[i,1],1]]])
            d1 = np.linalg.norm(ps[n,:] - xs[segs[i,0],:])
            d2 = np.linalg.norm(ps[n,:] - xs[segs[i,1],:])
            if ps[n,1] <= max(tempxs[:,1]) and ps[n,1] > min(tempxs[:,1]):
                x3 = ((tempxs[1,0]-tempxs[0,0])/(tempxs[1,1]-tempxs[0,1]))*(ps[n,1]-tempxs[0,1]) + tempxs[0,0]
                if ps[n,0] < x3:
                    wind=wind+1

        if wind%2 == 1:
            inn[n] = True

    return(inn)

def Draw_points(xs,app,canvas):
    nv = np.size(xs,axis=0)
    xs5 = 500*xs
    for n in range(nv-1):
        canvas.create_oval((xs5[n,0]-4, 500-xs5[n,1]-4, xs5[n,0]+4, 500-xs5[n,1]+4),
             fill='red')

    app.update()
    app.update_idletasks()

# draw curve on tkinter
def Draw_curve(xs, app, canvas):
    nv = np.size(xs,axis=0)
    xs5 = 500*xs
    for n in range(nv-1):
        canvas.create_line((xs5[n,0], 500-xs5[n,1], xs5[n+1,0], 500-xs5[n+1,1]),
             fill='black', width=2)

    canvas.create_line((xs5[-1,0], 500-xs5[-1,1], xs5[0,0], 500-xs5[0,1]),
             fill='black', width=2)
    app.update()
    app.update_idletasks()

def Draw_normal(xs,n,app,canvas):
    nv = np.size(xs,axis=0)
    xs5 = 500*xs
    n=n/5
    n5 = 500*n
    for i in range(nv):
        canvas.create_line(xs5[i,0], 500-xs5[i,1], xs5[i,0]+n5[i,1], 500-(xs5[i,1]-n5[i,0]),arrow=LAST)
    app.update()
    app.update_idletasks()


# draw mesh
def Draw_mesh(Tris, app, canvas):
    canvas.delete("all")
    nv = np.size(Tris["vertices"], axis=0)
    nelems = np.size(Tris["triangles"], axis=0)
    xs5 = 500*Tris["vertices"]
    edge = np.array([[0,1],[1,2],[2,0]])
    for i in range(nelems):
        for j in range(3):
            canvas.create_line((xs5[Tris["triangles"][i,edge[j,0]],0], 500-xs5[Tris["triangles"][i,edge[j,0]],1], 
            xs5[Tris["triangles"][i,edge[j,1]],0], 500-xs5[Tris["triangles"][i,edge[j,1]],1]),
            fill='black', width=2)
    for n in range(nv):
        canvas.create_oval(xs5[n,0]-3,500-xs5[n,1]-3,xs5[n,0]+3,500-xs5[n,1]+3,fill='red')
    app.update()
    app.update_idletasks() 


def draw_projections(app, canvas, spline):
    global xs 
    global flag
    flag = False
    xs = [[0.0, 0.0]]
    global stop
    stop = False
    def get_x_and_y(event):
        global lasx, lasy, xs
        lasx, lasy = event.x, event.y
        xs = np.array([event.x/500, (500-event.y)/500])
        canvas.create_oval(event.x-5, event.y-5, event.x+5, event.y+5,fill='blue')
        xsp = spline_proj(xs, spline)
        canvas.create_oval(500*xsp[0]-4, 500*(1-xsp[1])-4, 500*xsp[0]+4, 500*(1-xsp[1])+4,fill='blue')
        canvas.create_line(event.x,event.y,500*xsp[0],500*(1-xsp[1]),fill='black',width=2)
                


    def end(event):
        global lasx, lasy, xs, flag
        flag = True
        

    
    canvas.bind("<Button-1>", get_x_and_y)
    canvas.bind("<Button-3>", end)
    
    while flag == False:
        app.update_idletasks()
        app.update()

    app.update_idletasks()
    app.update()

    return(app, canvas)

def draw_curve_on_canvas(app):

    global xs 
    global flag
    flag = False
    xs = [[0.0, 0.0]]
    global stop
    stop = False
    def get_x_and_y(event):
        global lasx, lasy, xs
        lasx, lasy = event.x, event.y
        xs = np.append(xs, [[event.x, event.y]], axis=0)


    def draw_smth(event):
        global lasx, lasy, xs
        canvas.create_line((lasx, lasy, event.x, event.y), 
                      fill='black', 
                      width=2)
        lasx, lasy = event.x, event.y
        xs = np.append(xs, [[event.x, event.y]], axis=0)

    def end(event):
        global lasx, lasy, xs, flag
        xs = np.delete(xs,0,0)
        xs[:,1] = 500-xs[:,1]
        xs = xs/500
        flag = True
        

    
    canvas = Canvas(app, bg="white")
    canvas.pack(anchor='nw', fill='both', expand=1)
    canvas.bind("<Button-1>", get_x_and_y)
    canvas.bind("<B1-Motion>", draw_smth)
    canvas.bind("<Button-3>", end)
    
    while flag == False:
        app.update_idletasks()
        app.update()

    canvas.delete("all")
    app.update_idletasks()
    app.update()

    return(xs, app, canvas)