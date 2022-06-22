import numpy as np
import matplotlib.pyplot as plt

# Class for a box to be used later
class Box:
    def __init__(self, cx, cy, w, h):
        self.cx = cx
        self.cy = cy

        self.left = cx - w/2
        self.right = cx + w/2
        self.bottom = cy - h/2
        self.top = cy + h/2

        self.w = w
        self.h = h

    def contains(self, xs, n):
        assert(len(xs[n,:]) == 2)

        return (xs[n,0] >= self.left and
        xs[n,0] < self.right and
        xs[n,1] >= self.bottom and
        xs[n,1] < self.top)

    def intersect(self, other):
        return not (other.left > self.right or
        other.right < self.left or
        other.top < self.bottom or
        other.bottom > self.top)

    def draw(self,ax,c='k',linewidth=1, **kwargs):
        x0 = self.left
        y1 = self.top
        x1 = self.right
        y0 = self.bottom
        ax.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],c=c,lw=linewidth, **kwargs)

class QuadTree:
    def __init__(self, boundary, max_points, depth=0):

        self.boundary = boundary
        self.max_points = max_points
        self.points = []
        self.depth = depth

        self.divided = False

    def divide(self):
        cx,cy = self.boundary.cx, self.boundary.cy
        w,h = self.boundary.w/2, self.boundary.h/2

        self.nw = QuadTree(Box(cx-w/2,cy-h/2,w,h), self.max_points, self.depth + 1)
        self.ne = QuadTree(Box(cx+w/2,cy-h/2,w,h), self.max_points, self.depth + 1)
        self.se = QuadTree(Box(cx+w/2,cy+h/2,w,h), self.max_points, self.depth + 1)
        self.sw = QuadTree(Box(cx-w/2,cy+h/2,w,h), self.max_points, self.depth + 1)

        self.divided = True

    def insert(self, xs, n):

        if not self.boundary.contains(xs,n):
            return False

        if len(self.points) < self.max_points:
            self.points.append(n)
            return True

        if not self.divided:
            self.divide()
        
        return (self.ne.insert(xs,n) or
        self.nw.insert(xs,n) or
        self.se.insert(xs,n) or
        self.sw.insert(xs,n))

    def query_circle(self, boundary, centre, radius, found_points, xs):

        if not self.boundary.intersect(boundary):
            return found_points

        for point in self.points:
            if (boundary.contains(xs,point) and np.linalg.norm(xs[point,:] - centre) <= radius):
                found_points.append(point)

        # Recurse the search into this node's children.
        if self.divided:
            self.nw.query_circle(boundary, centre, radius, found_points, xs)
            self.ne.query_circle(boundary, centre, radius, found_points, xs)
            self.se.query_circle(boundary, centre, radius, found_points, xs)
            self.sw.query_circle(boundary, centre, radius, found_points, xs)
        return found_points

    def query_radius(self, centre, radius, found_points, xs):
        boundary = Box(centre[0], centre[1], 2*radius, 2*radius)
        return self.query_circle(boundary, centre, radius, found_points, xs)

    def find_closest_point(self, centre, h, found_points,xs):
        while not found_points:
            search = self.query_radius(centre,h,found_points,xs)
            h = 2*h
        
        dist = 1000000000.0
        for n in search:
            test_dist = np.linalg.norm(centre-xs[n,:])
            if test_dist < dist:
                dist = test_dist
                closest = n
        
        return closest


    def draw(self,ax):
        self.boundary.draw(ax)
        if self.divided:
            self.nw.draw(ax)
            self.ne.draw(ax)
            self.se.draw(ax)
            self.sw.draw(ax)

def points2quadtree(xs,maxpnts):
    r = max(xs[:,0])
    l = min(xs[:,0])
    b = min(xs[:,1]) 
    t = max(xs[:,1])
    w = r-l
    h = t-b
    cx = (r+l)/2
    cy = (t+b)/2

    Rect = Box(cx,cy,w,h)
    QT = QuadTree(Rect,maxpnts)
    for i in range(np.size(xs,0)):
        QT.insert(xs,i)
        
    return QT

def ForceField(xs,QT,r,Force):
    F = np.zeros((np.size(xs,0),2),dtype='float')
    for ii in range(np.size(xs,0)):
        points_in_circle = QT.query_radius(xs[ii,:],r,[],xs)
        for p in points_in_circle:
            if ii != p:
                vector = xs[p,:] - xs[ii,:]
                radius = np.linalg.norm(vector)
                F[p,:] = F[p,:] + Force(radius)*vector
           
    for ii in range(np.size(xs,0)):
        if xs[ii,0]+F[ii,0] > 0 and xs[ii,0]+F[ii,0] < 1 \
            and xs[ii,1]+F[ii,1] > 0 and xs[ii,1]+F[ii,1] < 1:
                xs[ii,:] = xs[ii,:]+F[ii,:]
        
    return xs
            
            
    

# random points to test certain class functions using matplotlib
nv = 150
for n in range(nv):
    xs = np.random.rand(nv,2) 
        

QT = points2quadtree(xs, 1)
Force = lambda r: 0.01*(1/nv)/(r**2)

fig, ax = plt.subplots()
for _ in range(50):
    ax.cla()
    ax.set(xlim=(0,1), ylim=(0,1))
    ax.scatter(xs[:,0],xs[:,1],5)
    QT = points2quadtree(xs, 4)
    QT.draw(ax)
    xs = ForceField(xs, QT, 0.2, Force)
    plt.pause(0.1)
