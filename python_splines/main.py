from drawpoints import *
from tkinter import *
from splines import *
import numpy as np
import matplotlib.pyplot as plt

def main_subroutine(app, degree):
    xs, app, canvas = draw_on_canvas(app)

    spline = spline_init(xs, degree)

    n = 100
    ps = np.zeros((n,2))
    grad = np.zeros((n,2))
    for i in range(n):
        ps[i,:], grad[i,:] = spline_var(float(i/n),spline)
    
    Draw_curve(ps, app, canvas)
    Draw_normal(ps,grad,app,canvas)

    app,canvas = draw_projections(app,canvas,spline)
    return(app,canvas)



global app, canvas
def start_main():
    global app, canvas, degree
    canvas.destroy()
    app, canvas = main_subroutine(app, degree)

degree = int(sys.argv[1])
app  = Tk()
app.geometry("500x500")
app.title("Draw to mesh")
canvas = Canvas(app)
B = Button(app,text = "Reset",command = start_main)
B.pack()
app.mainloop()