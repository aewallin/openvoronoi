#import ocl
import openvoronoi as ovd
import ovdvtk
#import time
import vtk
#import datetime
import math
#import random

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

"""
def drawArc(myscreen, p1, p2, c, cw, arcColor):
    r = (p1-c).norm()
    
    pass
"""
    
# draw line  a x + b y + c = 0
# draws lines roughly in a 100x100 box (?)
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






def drawArcPredicate():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  100 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    c1 = ovd.Point(0,0)
    r1=20
    
    alfa1 = (float(23)/float(360) )* 2*math.pi
    alfa2 = (float(123)/float(360) )* 2*math.pi
    alfa1_dir = ovd.Point( math.cos(alfa1), math.sin(alfa1) )
    alfa2_dir = ovd.Point( math.cos(alfa2), math.sin(alfa2) )
    p1 = c1 + r1*alfa1_dir
    p2 = c1 + r1*alfa2_dir
    
    drawVertex(myscreen, c1, ovdvtk.orange, rad=1)
    drawVertex(myscreen, p1, ovdvtk.green, rad=1)
    drawVertex(myscreen, p2, ovdvtk.red, rad=1)    
    
    #drawArc(myscreen, p1, p2, c1, True, ovdvtk.yellow)
    ovdvtk.drawArc(myscreen, p1, p2, (p1-c1).norm(), c1, True, ovdvtk.yellow, da=0.1)
    ovdvtk.drawArc(myscreen, p1, p2, (p1-c1).norm(), c1, False, ovdvtk.orange, da=0.1)
    #c2 = ovd.Point(6,13)
    
    #r2=23
    # external ma-pocket fails with this input:
    """
    c1 = ovd.Point(0, -25)
    r1 = 5 #16.7033
    c2 = ovd.Point(0, 0.2167)
    r2 = 15.9299
    r2 = 2 #16 # 5
    drawCircle(myscreen, c1, r1, ovdvtk.red)
    drawCircle(myscreen, c2, r2, ovdvtk.green)
    # when machining c2 the maximum cut-width is 
    # w_max = | c2 - c1 | + r2 - r1
    #print "dr = ",dr
    print " cut-width = ", ((c2-c1).norm()+r2-r1)
    #[bd1,bd2] = bitanget_direction(c1,c2,r1,r2)
    # from C, go a distance r along the normal to the line.
    #p1 = c1 + r1*bd1
    #p2 = c1 + r1*bd2
    drawCircle(myscreen,p1,0.5,ovdvtk.pink)
    drawCircle(myscreen,p2,0.5,ovdvtk.red)

    p3 = c2 + r2*bd1
    p4 = c2 + r2*bd2
    drawCircle(myscreen,p3,0.5,ovdvtk.lgreen)
    drawCircle(myscreen,p4,0.5,ovdvtk.green)
    """
    myscreen.render()
    myscreen.iren.Start()
    
if __name__ == "__main__":  
    drawArcPredicate()
