#import ocl
import openvoronoi as ovd
import ovdvtk
#import time
import vtk
#import datetime
import math
import random

"""
This script does not use OpenVoronoi, it is used merely for drawing and
verifying solutions etc.
There are functions for:
- voronoi-bisectors and verifying the parametric equations for bisectors
- voronoi-vertex solvers (separator-solver)
- bi-tangent lines/points (for medial-axis pocketing)


"""

def drawVertex(myscreen, p, vertexColor, rad=1):
    myscreen.addActor( ovdvtk.Sphere( center=(p.x,p.y,0), radius=rad, color=vertexColor ) )

def drawEdge(myscreen, e, edgeColor=ovdvtk.yellow):
    p1 = e[0]
    p2 = e[1]
    myscreen.addActor( ovdvtk.Line( p1=( p1.x,p1.y,0), p2=(p2.x,p2.y,0), color=edgeColor ) )

def drawCircle(myscreen, c, circleColor):
    myscreen.addActor( ovdvtk.Circle( center=(c.c.x,c.c.y,c.c.z), radius=c.r, color=circleColor ) )

def drawCircle(myscreen, c, r, circleColor):
    myscreen.addActor( ovdvtk.Circle( center=(c.x,c.y, 0), radius=r, color=circleColor ) )
    
def drawLine(myscreen, p1,p2, lineColor=ovdvtk.yellow):
    myscreen.addActor( ovdvtk.Line( p1=( p1.x,p1.y,0), p2=(p2.x,p2.y,0), color=lineColor ) )
    
"""
def drawArc(myscreen, p1, p2, c, cw, arcColor):
    r = (p1-c).norm()
    
    pass
"""
    
# draw line  a x + b y + c = 0
# draws lines roughly in a 100x100 box (?)
"""
def drawLine(myscreen, l, lineColor):
    #  x = -by-c / a
    #if l.a != 0:
    if (abs(l.a) > abs(l.b)):
        y=100
        p1 = ovd.Point( float(-l.b*y-l.c)/l.a , y )
        y=-100
        p2 = ovd.Point( float(-l.b*y-l.c)/l.a , y )
        myscreen.addActor( ovdvtk.Line( p1=( p1.x,p1.y,0), p2=(p2.x,p2.y,0), color=lineColor ) )
    else:
        x=100
        p1 = ovd.Point( x, float(-l.a*x-l.c)/l.b )
        x=-100
        p2 = ovd.Point( x, float(-l.a*x-l.c)/l.b )
        myscreen.addActor( ovdvtk.Line( p1=( p1.x,p1.y,0), p2=(p2.x,p2.y,0), color=lineColor ) )
"""

# CIRCLE definition
# circle offset is  (x(t) - xc)^2 + (y(t)-yc)^2 = (r+k*t)^2
# POINT is circle with r=1 and k=1 
class Circle:
    def __init__(self,c=ovd.Point(0,0),r=1,cw=1,k=1):
        self.c = c
        self.r = r
        self.cw = cw # CW=1, CCW = -1
        self.k = k # offset direction

# LINE def
# line offset is  a1 x + b1 y + c1 + k1 t = 0 and a*a + b*b = 1 
class Line:
    def __init__(self,a,b,c,k):
        self.a = float(a)
        self.b = float(b)
        det = self.a*self.a+self.b*self.b
        self.c = float(c)
        self.k = float(k) # offset to left or right of line





def arc_in_region(p1,p2,c,cw,p):
    if cw:
        return p.is_right(c,p1) and (not p.is_right(c,p2))
    else:
        return (not p.is_right(c,p1)) and p.is_right(c,p2)

    
def drawArcPredicate():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  100 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    c1 = ovd.Point(0,0)
    r1=20
    
    alfa1 = (float(23)/float(360) )* 2*math.pi
    alfa2 = (float(123)/float(360) )* 2*math.pi
    
    alfa2, alfa1 = alfa1, alfa2 # swap
    
    alfa1_dir = ovd.Point( math.cos(alfa1), math.sin(alfa1) )
    alfa2_dir = ovd.Point( math.cos(alfa2), math.sin(alfa2) )
    p1 = c1 + r1*alfa1_dir
    p2 = c1 + r1*alfa2_dir
    cw = True
    
    drawVertex(myscreen, c1, ovdvtk.orange, rad=1)
    fa1 = ovdvtk.FollowerText(text="c",color=ovdvtk.orange,center=(c1.x+1,c1.y,0),scale=1)
    myscreen.addActor(fa1)
    
    drawVertex(myscreen, p1, ovdvtk.green, rad=1)
    fa2 = ovdvtk.FollowerText(text="p1",color=ovdvtk.green,center=(p1.x+1,p1.y,0),scale=1)
    myscreen.addActor(fa2)
    
    drawVertex(myscreen, p2, ovdvtk.red, rad=1)    
    fa3 = ovdvtk.FollowerText(text="p2",color=ovdvtk.red,center=(p2.x+1,p2.y,0),scale=1)
    myscreen.addActor(fa3)
    
    #drawArc(myscreen, p1, p2, c1, True, ovdvtk.yellow)
    #ovdvtk.drawArc(myscreen, p1, p2, (p1-c1).norm(), c1, True, ovdvtk.yellow, da=0.1)
    ovdvtk.drawArc(myscreen, p1, p2, (p1-c1).norm(), c1, cw, ovdvtk.orange, da=0.1)
    
    Nmax = 5000
    for n in range(Nmax):
        p = 100*ovd.Point( random.random()-0.5, random.random()-0.5 )
        if arc_in_region(p1,p2,c1,cw,p):
            drawVertex(myscreen, p, ovdvtk.lgreen, rad=0.1)
        else:
            drawVertex(myscreen, p, ovdvtk.pink, rad=0.1)

    myscreen.render()
    myscreen.iren.Start()

def closer_endpoint(p1,p2,p):
    if (p1-p).norm() < (p2-p).norm():
        return p1
    else:
        return p2
        
def projection_point(p1,p2,c1,cw,p):
    if p==c1:
        return p1
    else:
        n = (p-c1)
        n.normalize()
        return c1 + (p1-c1).norm()*n
        
def apex_point(p1,p2,c1,cw,p):
    if arc_in_region(p1,p2,c1,cw,p):
        return projection_point(p1,p2,c1,cw,p)
    else:
        return closer_endpoint(p1,p2,p)
            
def drawArcPredicate2():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  100 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    c1 = ovd.Point(0,0)
    r1=20
    
    alfa1 = (float(23)/float(360) )* 2*math.pi
    alfa2 = (float(123)/float(360) )* 2*math.pi
    
    #alfa2, alfa1 = alfa1, alfa2 # swap
    
    alfa1_dir = ovd.Point( math.cos(alfa1), math.sin(alfa1) )
    alfa2_dir = ovd.Point( math.cos(alfa2), math.sin(alfa2) )
    p1 = c1 + r1*alfa1_dir
    p2 = c1 + r1*alfa2_dir
    cw = False
    
    drawVertex(myscreen, c1, ovdvtk.orange, rad=1)
    fa1 = ovdvtk.FollowerText(text="c",color=ovdvtk.orange,center=(c1.x+1,c1.y,0),scale=1)
    myscreen.addActor(fa1)
    
    drawVertex(myscreen, p1, ovdvtk.green, rad=1)
    fa2 = ovdvtk.FollowerText(text="p1",color=ovdvtk.green,center=(p1.x+1,p1.y,0),scale=1)
    myscreen.addActor(fa2)
    
    drawVertex(myscreen, p2, ovdvtk.red, rad=1)    
    fa3 = ovdvtk.FollowerText(text="p2",color=ovdvtk.red,center=(p2.x+1,p2.y,0),scale=1)
    myscreen.addActor(fa3)
    
    #drawArc(myscreen, p1, p2, c1, True, ovdvtk.yellow)
    #ovdvtk.drawArc(myscreen, p1, p2, (p1-c1).norm(), c1, True, ovdvtk.yellow, da=0.1)
    ovdvtk.drawArc(myscreen, p1, p2, (p1-c1).norm(), c1, cw, ovdvtk.orange, da=0.1)
    
    Nmax = 5000
    for n in range(Nmax):
        p = 100*ovd.Point( random.random()-0.5, random.random()-0.5 )
        
        apex = apex_point(p1,p2,c1,cw,p)
        linecolor= ovdvtk.pink
        if arc_in_region(p1,p2,c1,cw,p):
            linecolor = ovdvtk.lgreen
        
        drawLine(myscreen,p,apex,linecolor)

    myscreen.render()
    myscreen.iren.Start()

if __name__ == "__main__":
    #drawArcPredicate()
    drawArcPredicate2()
