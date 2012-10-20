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

# from Held 1991, page 94->
#
# bisectors are of the form
# line, parabola, ellipse, hyperbola
#     x(t) = x1 - x2 - x3*t +/- x4 sqrt( square(x5+x6*t) - square(x7+x8*t) )
#     y(t) = y1 - y2 - y3*t +/- y4 sqrt( square(y5+y6*t) - square(y7+y8*t) ) 
# line/line: line
# circle/line: parabola
# circle/circle: ellipse/hyperbola
# !only valid if no parallel lines and no concentric arcs
#
# line:   (a, b, c, k)     
#         a1 x + b1 y + c1 + k1 t = 0 and a*a + b*b = 1
#         k= +/- 1 indicates offset to right/left
#
# circle: (xc, yc, r, lambda) 
#         (x(t) - xc1)^2 + (y(t)-yc1)^2 = (r1+k1*t)^2
#        lambda=-1 for CCW arc and +1 otherwise 
#         k +1 enlarging circle, k -1 shrinking circle
#
# for a bisector we store only four parameters (alfa1, alfa2, alfa3, alfa4)
#
# line/line
#    delta = a1*b2 - a2*b1
#    alfa1= (b1*d2-b2*d1)/delta
#    alfa2= (a2*d1-a1*d2)/delta
#    alfa3= b2-b1
#    alfa4= a1-a2 
#   bisector-params:
#    x1 = alfa1, x3 = -alfa3, x2 = x4 = x5 = x6 = x7 = x8 = 0
#    y1 = alfa2, y3 = -alfa4, y2=y4=y5=y6=y7=y8 = 0 
#
# circle/line
#     
#    alfa1= a2
#    alfa2= b2
#    alfa3= a2*xc1 + b2*yc1+d2
#    alfa4= r1
#   params: 
#    x1 = xc1, x2 = alfa1*alfa3, x3 = -alfa1, x3 = alfa2, x5 = alfa4, x6 = lambda1, x7 = alfa3, x8 = -1
#    y1 = yc1, y2 = alfa2*alfa3, y3 = -alfa2, y4 = alfa1, y5 = alfa4, y6 = lambda1, y7 = alfa3, y8 = -1 
#
# circle / circle
#    d= sqrt( square(xc1-xc2) + square(yc1-yc2) ) 
#    alfa1= (xc2-xc1)/d
#    alfa2= (yc2-yc1)/d
#    alfa3= (r2*r2-r1*r1-d*d)/2d
#    alfa4= (lambda2*r2-lambda1*r1)/d 
#   params:
#    x1 = xc1, x2 = alfa1*alfa3, x3 = alfa1*alfa4, x4 = alfa2, x5 = r1, x6 = lambda1, x7 = alfa3, x8 = alfa4
#    y1 = yc1, y2 = alfa2*alfa3, y3 = alfa2*alfa4, y4 = alfa1, y5 = r1, y6 = lambda1, y7 = alfa3, y8 = alfa4 

class LineLine:
    """ line/line bisector is a line """
    def __init__(self,l1,l2):
        self.delta = l1.a*l2.b-l2.a*l1.b
        self.alfa1 = (l1.b*l2.c-l2.b*l1.c) / self.delta # intersection x-coord
        self.alfa2 = (l2.a*l1.c-l1.a*l2.c) / self.delta # intersection y-coord
        self.alfa3 = (l2.b*l1.k-l1.b*l2.k) / ( self.delta ) # x-direction
        self.alfa4 = (l1.a*l2.k-l2.a*l1.k) / ( self.delta ) # y-direction
        #self.alfa3 = l2.b-l1.b
        #self.alfa4 = l1.a-l2.a
        if self.delta>0:
            self.sign = +1
        else:
            self.sign = +1 
        self.tmin = float(0)
        self.tmax = float(20)
    def __str__(self):
        s1 = "alfa1= "+str(self.alfa1)+" alfa3="+str(self.alfa3)
        s2 = "alfa2= "+str(self.alfa2)+" alfa4="+str(self.alfa4)
        return s1+s2

    def getX(self):
        x = []
        x.append( self.alfa1 )
        x.append( 0 )
        x.append( -self.alfa3*self.sign )
        x.append( 0 )
        x.append( 0 )
        x.append( 0 )
        x.append( 0 )
        x.append( 0 )
        return x
    def getY(self):
        y = []
        y.append( self.alfa2 )
        y.append( 0 )
        y.append( -self.alfa4*self.sign )
        y.append( 0 )
        y.append( 0 )
        y.append( 0 )
        y.append( 0 )
        y.append( 0 )
        return y

class CircleCircle:
    # CIRCLE / CIRCLE
    # d= sqrt( square(xc1-xc2) + square(yc1-yc2) ) 
    # cw=-1 for CCW arc and +1 otherwise 
    def __init__(self, c1, c2):
        self.d     = (c1.c-c2.c).norm() 
        self.alfa1 = 0.0
        self.alfa2 = 0.0
        self.alfa3 = 0.0
        self.alfa4 = 0.0            
        if ( self.d > 0.0 ):
            self.alfa1 = (c2.c.x-c1.c.x)/self.d
            self.alfa2 = (c2.c.y-c1.c.y)/self.d
            self.alfa3 = (c2.r*c2.r-c1.r*c1.r-self.d*self.d)/(2*self.d)
            self.alfa4 = (c2.cw*c2.r-c1.cw*c1.r)/self.d # k1 k2 (expand/contract??)
        self.c1 = c1 # store all of c1 also??
        #print " alfa4:",self.alfa4
        self.tmin = ((c1.c-c2.c).norm()-c2.r-c1.r) / 2 
        self.tmax = 100
        
    def getX(self):
        x = []
        x.append( self.c1.c.x )
        x.append( self.alfa1*self.alfa3 )
        x.append( self.alfa1*self.alfa4 )
        x.append( self.alfa2 )
        x.append( self.c1.r  )
        x.append( self.c1.cw ) # expand/contract?
        x.append( self.alfa3 )
        x.append( self.alfa4 )
        return x
    def getY(self):
        y = []
        y.append( self.c1.c.y )
        y.append( self.alfa2*self.alfa3 )
        y.append( self.alfa2*self.alfa4 )
        y.append( self.alfa1 )
        y.append( self.c1.r  )
        y.append( self.c1.cw )
        y.append( self.alfa3 )
        y.append( self.alfa4 )
        return y

# CIRCLE/LINE (same as point-line?)
#    * alfa1= a2
#    * alfa2= b2
#    * alfa3= a2*xc1 + b2*yc1+d2 (c2?)
#    * alfa4= r1 
#       x1 = xc1
#       x2 = alfa1*alfa3
#       x3 = -alfa1, 
#       x3 = alfa2, 
#       x5 = alfa4, 
#       x6 = lambda1, 
#       x7 = alfa3, 
#       x8 = -1
#       y1 = yc1, 
#       y2 = alfa2*alfa3, 
#       y3 = -alfa2, 
#       y4 = alfa1, 
#       y5 = alfa4, 
#       y6 = lambda1, 
#       y7 = alfa3, 
#       y8 = -1 

class CircleLine:
    def __init__(self, c1, l2):
        #self.delta = l2.a*c1.c.x + l2.b*c1.c.y + l2.c 
        self.alfa1 = l2.a
        self.alfa2 = l2.b
        self.alfa3 = l2.a*c1.c.x + l2.b*c1.c.y + l2.c
        self.alfa4 = c1.r
        self.c1 = c1 # store all of c1 also??
        self.k = +1 # positive line-offset
        if self.alfa3 > 0:
            self.k = -1 # negative line-offset
        
        # alfa3 is distance from circle-center to line
        self.tmin = (abs(self.alfa3)-c1.r)/(2.0)
        
        self.tmax = 100
        print " tmax=",self.tmax," tmin=",self.tmin
    def getX(self):
        x = []
        x.append( self.c1.c.x )           # x0 = xc1
        x.append( self.alfa1*self.alfa3 ) # x1 = a2*h
        x.append( self.alfa1*self.k )     # x2 = a2*k2
        x.append( self.alfa2 )            # x3 = b2
        x.append( self.alfa4 )            # x4 = 0 for Point!
        x.append( self.c1.cw )            # x5 = 1 for Point
        x.append( self.alfa3 )            # x6 = h
        x.append(  self.k )               # x7, line-offset
        return x
    def getY(self):
        y = []
        y.append( self.c1.c.y )
        y.append( self.alfa2*self.alfa3 )
        y.append( self.alfa2*self.k )
        y.append( self.alfa1 )
        y.append( self.alfa4  )
        y.append( self.c1.cw )
        y.append( self.alfa3 )
        y.append( self.k ) 
        return y

# this is the general bisector representation that can be used to get points
# and draw all bisectors ( POINTPOINT, POINTLINE, LINELINE, etc)
class Bisector:
    def __init__(self, Bis):
        self.bis = Bis
        self.x= Bis.getX()
        self.y= Bis.getY()
        #print self.y
        #print self.x
    def Point(self, t, k):
        x=self.x
        y=self.y
        detx = ( math.pow((x[4]+x[5]*t),2) - math.pow((x[6]+x[7]*t),2) )
        dety = ( math.pow((y[4]+y[5]*t),2) - math.pow((y[6]+y[7]*t),2) )
        xp = x[0]-x[1]-x[2]*t + x[3]*math.sqrt( detx ) # plus sign
        yp = y[0]-y[1]-y[2]*t - y[3]*math.sqrt( dety )
        xm = x[0]-x[1]-x[2]*t - x[3]*math.sqrt( detx ) # minus sign
        ym = y[0]-y[1]-y[2]*t + y[3]*math.sqrt( dety ) 
        return [ovd.Point(xp,yp), ovd.Point(xm,ym)]
        """
    def minT(self):
        # the minimum t that makes sense sets the sqrt() to zero
        # (x[4]+x[5]*t)^2 - (x[6]+x[7]*t)^2 = 0
        # (x[4]+x[5]*t)^2 = (x[6]+x[7]*t)^2
        # (x[4]+x[5]*t) = (x[6]+x[7]*t)  OR   (x[4]+x[5]*t) = -(x[6]+x[7]*t)
        # (x[5]-x[7])*t = (x[6]-x[4])  OR   (x[5]+x7*t) = x4-x[6]
        # t = x6-x4 / (x5-x7)    or    t = x4-x6 / (x5+x7)
        x = self.x 
        y = self.y 
        t1=0
        t2=0
        t3=0
        t4=0
        if (((x[5]-x[7])!=0) and ((x[5]+x[7])!=0) ):
            t1 = (x[6]-x[4]) / (x[5]-x[7])
            t2 = (-x[4]-x[6]) / (x[5]+x[7])
            t3 = (y[6]-y[4]) / (y[5]-y[7])
            t4 = (-y[4]-y[6]) / (y[5]+y[7])
        else:
            print "NO solutions!"
        print " t1 solution= ",t1
        print " t2 solution= ",t2
        print " t3 solution= ",t3
        print " t4 solution= ",t4
        return t2
        """
    def trange(self,N):
        tmin = self.bis.tmin
        tmax = self.bis.tmax
        dt = (tmax-tmin)/N
        tlist=[]
        for i in range(N):
            t = tmin + i*dt
            tlist.append(t)
        #print tlist
        return tlist
        # compute range of t-values with N pts
        
        
def drawBisector(myscreen, bis):
    N = 200
    """
    t= bis.minT()
    tmax = 400
    dt = float(tmax)/float(N)
    """
    ppts = []
    mpts = []
    for t in bis.trange(N):
        ppts.append( bis.Point(t,1)[0] ) # with positive sign
        mpts.append( bis.Point(t,1)[1] ) # with negative sign
        #t= t+dt
    #for p in ppts:
    #    drawVertex(myscreen, p, ovdvtk.green, rad=1)
    for n in range( len(ppts)-1):
        p1=ppts[n]
        p2=ppts[n+1]
        myscreen.addActor( ovdvtk.Line( p1=( p1.x,p1.y,0), p2=(p2.x,p2.y,0), color=ovdvtk.green ) )
    
    #for p in mpts:
    #    drawVertex(myscreen, p, ovdvtk.red, rad=1)
    for n in range( len(mpts)-1):
        p1=mpts[n]
        p2=mpts[n+1]
        myscreen.addActor( ovdvtk.Line( p1=( p1.x,p1.y,0), p2=(p2.x,p2.y,0), color=ovdvtk.red ) )
        
def drawLinePointTest():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  1000 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    c1 = Circle(c=ovd.Point(10,-30), r=0, cw=+1) # first site is zero-radius circle, i.e. point-site
    drawVertex(myscreen, c1.c, ovdvtk.yellow)
    l1 = Line( math.cos(1),   math.sin(1)   , 1 , +1) # second site is a line
    drawLine(myscreen, l1, ovdvtk.yellow )
    c1l1 = CircleLine( c1, l1) # the bisector between the sites
    bill = Bisector( c1l1 )
    drawBisector( myscreen, bill )
    myscreen.render()    
    myscreen.iren.Start()

def drawPointPointTest():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  1000 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    c1 = Circle(c=ovd.Point(-10,20), r=0, cw=+1) # first site
    drawVertex(myscreen, c1.c, ovdvtk.yellow)
    c2 = Circle(c=ovd.Point(20,4),  r=0,cw=+1) # second site
    drawVertex(myscreen, c2.c, ovdvtk.yellow)
    c1l1 = CircleCircle( c1, c2) # bisector
    bill = Bisector( c1l1 )
    drawBisector( myscreen, bill )
    myscreen.render()
    myscreen.iren.Start()

def drawLineLineTest():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  1000 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    
    l1 = Line( math.cos(1),   math.sin(1)   , 1 , +1) # first line-site
    drawLine(myscreen, l1, ovdvtk.yellow )
    
    l2 = Line( math.cos(3),   math.sin(3)   , 1 , -1) # second line-site
    drawLine(myscreen, l2, ovdvtk.orange )
    
    l1l2 = LineLine( l1, l2 ) # bisectors
    l2l1 = LineLine( l2, l1 ) # it should not matter if we call with (l1,l2) or (l2,l1) (?)
    print l1l2
    print l2l1
    b1= Bisector( l1l2 )
    b2= Bisector( l2l1 )
    drawBisector( myscreen, b1 )
    drawBisector( myscreen, b2 )
    myscreen.render()
    myscreen.iren.Start()

def drawPointArcTest():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  1000 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    c1 = Circle(c=ovd.Point(-20,20), r=0, cw=+1) # first site
    drawVertex(myscreen, c1.c, ovdvtk.yellow)
    
    c2 = Circle(c=ovd.Point(30,50), r=20, cw=+1) # second site
    drawVertex(myscreen, c1.c, ovdvtk.orange)
    drawCircle(myscreen, c2.c, c2.r, ovdvtk.orange)
    #myscreen.addActor( ovdvtk.Circle( center=(c.c.x,c.c.y,c.c.z), radius=c.r, color=circleColor ) )
        
    c1c2 = CircleCircle( c1, c2 ) # bisectors
    
    b1= Bisector( c1c2 )
    
    drawBisector( myscreen, b1 )
    
    myscreen.render()
    myscreen.iren.Start()

def drawPointArcTest2():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  1000 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    c1 = Circle(c=ovd.Point(-10,20), r=0, cw=+1) # first site
    drawVertex(myscreen, c1.c, ovdvtk.yellow)
    

    
    c2 = Circle(c=ovd.Point(30,50), r=20, cw=+1) # second site
    drawVertex(myscreen, c1.c, ovdvtk.orange)
    drawCircle(myscreen, c2.c, c2.r, ovdvtk.orange)
    
    c1c2 = CircleCircle( c2, c1 ) # bisectors
    
    b1= Bisector( c1c2 )

    drawBisector( myscreen, b1 )

    myscreen.render()
    myscreen.iren.Start()

def drawLineArcTest():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  1000 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    l1 = Line( math.cos(1),   math.sin(1)   , -50 , +1) # first line-site
    drawLine(myscreen, l1, ovdvtk.yellow )
    
    c2 = Circle(c=ovd.Point(30,50), r=20, cw=+1) # first site
    drawVertex(myscreen, c2.c, ovdvtk.orange)
    drawCircle(myscreen, c2.c, c2.r, ovdvtk.orange)
    #myscreen.addActor( ovdvtk.Circle( center=(c.c.x,c.c.y,c.c.z), radius=c.r, color=circleColor ) )
        
    c1c2 = CircleLine( c2, l1 ) # bisectors
    #CircleLine
    b1= Bisector( c1c2 )
    #b2= Bisector( l2l1 )
    drawBisector( myscreen, b1 )
    #drawBisector( myscreen, b2 )
    myscreen.render()
    myscreen.iren.Start()
    
# s1= l1: line-site
# s2= p2: point-site (and end-point of l1)
# s3: l3: line-site
def drawSeparatorSolver1(alfa = 6.4):
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  100 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #             a              b            c    k       
    l1 = Line( math.cos(1),   math.sin(1)   , 1 , -1) # first line-site
    drawLine(myscreen, l1, ovdvtk.yellow )
    
    # pick a point on the line which represents the end-point
    y = 20
    p2 = ovd.Point( (-l1.b*y-l1.c)/l1.a , y )
    drawVertex(myscreen, p2, ovdvtk.orange , rad=1)
    myscreen.camera.SetPosition(p2.x+0.001, p2.y,  300 ) 
    myscreen.camera.SetFocalPoint(p2.x, p2.y, 0)
    
    # a second line-site
    #alfa = 4.1
    l3 = Line( math.cos(alfa),   math.sin(alfa)   , +50 , -1) # second line-site
    drawLine(myscreen, l3, ovdvtk.orange )
    
    sv = ovd.Point()
    if l1.k == -1:
        sv.x = l1.a
        sv.y = l1.b
    else:
        sv.x = -l1.a
        sv.y = -l1.b
    # draw rays representing positive/negative separator:
    pos_sep =[p2, p2+100*sv ]
    neg_sep =[p2, p2-100*sv ]
    drawEdge(myscreen, pos_sep, edgeColor=ovdvtk.green)
    drawEdge(myscreen, neg_sep, edgeColor=ovdvtk.red)
    
    # 1st degree equation gives t directly:
    tsln = -(l3.a*p2.x+l3.b*p2.y+l3.c) / ( sv.x*l3.a + sv.y*l3.b + l3.k  )
    print tsln
    psln = p2 + tsln*sv
    drawVertex(myscreen, psln, ovdvtk.pink , rad=1)
    drawCircle(myscreen, psln, tsln, ovdvtk.pink )
    drawCircle(myscreen, p2, tsln, ovdvtk.red )
    
    myscreen.render()
    myscreen.iren.Start()

# s1= l1: line-site
# s2= p2: point-site (and end-point of l1)
# s3: p3: point-site
def drawSeparatorSolver2(px=10,py=20):
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  100 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #             a              b            c    k       
    l1 = Line( math.cos(1),   math.sin(1)   , 1 , -1) # first line-site
    drawLine(myscreen, l1, ovdvtk.yellow )
    
    # pick a point on the line which represents the end-point
    y = 20
    p2 = ovd.Point( (-l1.b*y-l1.c)/l1.a , y )
    drawVertex(myscreen, p2, ovdvtk.orange , rad=1)
    myscreen.camera.SetPosition(p2.x+0.001, p2.y,  300 ) 
    myscreen.camera.SetFocalPoint(p2.x, p2.y, 0)
    
    
    # a second point-site
    #alfa = 4.1
    #l3 = Line( math.cos(alfa),   math.sin(alfa)   , +50 , -1) # second line-site
    p3 = ovd.Point(px,py)
    drawVertex(myscreen, p3, ovdvtk.yellow , rad=1)
    #drawLine(myscreen, l3, ovdvtk.orange )
    sv = ovd.Point()
    if l1.k == -1:
        sv.x = l1.a
        sv.y = l1.b
    else:
        sv.x = -l1.a
        sv.y = -l1.b
    # draw rays representing positive/negative separator:
    pos_sep =[p2, p2+100*sv ]
    neg_sep =[p2, p2-100*sv ]
    drawEdge(myscreen, pos_sep, edgeColor=ovdvtk.green)
    drawEdge(myscreen, neg_sep, edgeColor=ovdvtk.red)
    
    # 2nd degree equation gives t :
    #  A t^2 +  B t + C = 0
    dx = p2.x - p3.x
    dy = p2.y - p3.y
    tsln = -(dx*dx+dy*dy) / (2*( dx*sv.x+dy*sv.y  ))

    print tsln
    psln = p2 + tsln*sv
    drawVertex(myscreen, psln, ovdvtk.pink , rad=1)
    drawCircle(myscreen, psln, tsln, ovdvtk.pink )
    drawCircle(myscreen, p3, tsln, ovdvtk.orange )
    drawCircle(myscreen, p2, tsln, ovdvtk.red )
    myscreen.render()
    myscreen.iren.Start()

# solve: a x^2 + b x + c = 0
# copied from c-code in numeric.hpp
def quadratic_roots( a,  b,  c):
    roots = []
    if ((a == 0) and (b == 0)):
        return roots;
    
    if (a == 0):
        roots.append( -c / b )
        return roots
    
    if (b == 0):
        sqr = -c / a
        if (sqr > 0): 
            roots.append( math.sqrt(sqr) )
            roots.append( -roots[0] )
            return roots
        elif (sqr == 0): 
            roots.append( 0 )
            return roots
        else:
            return roots
        
    disc = b*b - 4*a*c # discriminant
    if (disc > 0): 
        q=0
        if (b > 0):
            q = (b + math.sqrt(disc)) / -2
        else:
            q = (b - math.sqrt(disc)) / -2
        roots.append( q / a )
        roots.append( c / q ) 
        return roots
    elif (disc == 0): 
        roots.append( -b / (2*a) )
        return roots

    return roots # no real roots
    
    
def drawBitangents():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  100 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    """
    c1 = ovd.Point(10,20)
    r1=20
    c2 = ovd.Point(6,13)
    r2=23
    """
    
    # external ma-pocket fails with this input:
    c1 = ovd.Point(0, -18)
    r1 = 16.7033
    c2 = ovd.Point(0, -17.2167)
    r2 = 15.9299


    drawCircle(myscreen, c1, r1, ovdvtk.red)
    drawCircle(myscreen, c2, r2, ovdvtk.green)
    
    # when machining c2 the maximum cut-width is 
    # w_max = | c2 - c1 | + r2 - r1
    
    # we are looking for the bi-tangents of the two circles, given by
    # A x + B y + C = 0
    # with the constraint A^2 + B^2 = 1
    # the line-should be a bi-tangent, so:
    # +/- r1 = A c1x + B c1y + C
    # +/- r2 = A c2x + B c2y + C
    # which in matrix form is
    #   Mc=v  
    # ( c1x c1y ) ( A ) = ( +/- r1 - C )
    # ( c2x c2y ) ( B ) = ( +/- r2 - C )
    #
    # now solve for A and B using Cramers rule
    #
    #  A = det(M_1) / det(M) =  det ( +/-r1-C   c1y ) / det(M)  =  mC + n
    #                               ( +/-r2-C   c2y )
    # and
    #  B = det(M_2) / det(M) = det ( c1x  +/-r1-C ) / det(M) = pC + q
    #                              ( c2x  +/-r2-C )
    # this can be inserted into the constraint
    # A^2 + B^2 = 1 
    # (mC + n )^2 + (pC+q)^2
    # m^2 C^2 + 2mnC + n^2 + p^2 C^2 + 2pqC + q^2 =1
    #
    # C^2 ( m^2 + p^2 ) + C 2(mn+pq) + n^2 + q^2 = 1
    #
    # +r1 +r2 case:
    #  det( r1-C c1y ) = (r1-C)c2y - (r2-C)c1y = C (c1y-c2y) + (c2y*r1 - c1y*r2)
    #     ( r2-C c2y )
    #  det( c1x r1-C ) = c1x(r2-C) - c2x(r1-C) = C ( c2x-c1x ) + (c1x*r2 - c2x*r1)
    #     ( c2x r2-C )
    lines = []
    detM = c1.x*c2.y - c2.x*c1.y
    print "detM= ",detM
    m = ( c1.y-c2.y ) / detM
    p = ( c2.x-c1.x ) / detM
    
    n = ( c2.y*r1 - c1.y*r2 ) / detM
    q = ( c1.x*r2 - c2.x*r1 ) / detM
    roots1 = quadratic_roots( m*m+p*p, 2*(m*n+p*q),  n*n+q*q-1)
    print "roots1 ", roots1
    for r in roots1:
        lines.append( root_to_line(r,m,n,p,q) )
        
    """
    n = ( c2.y*(-1)*r1 - c1.y*(-1)*r2 ) / detM
    q = ( c1.x*(-1)*r2 - c2.x*(-1)*r1 ) / detM
    roots2 = quadratic_roots( m*m+p*p, 2*(m*n+p*q),  n*n+q*q-1)
    print "roots2 ", roots2
    for r in roots2:
        lines.append( root_to_line(r,m,n,p,q) )
    """

    for l in lines:
        drawLine( myscreen, l, ovdvtk.yellow )
        
    # now figure out the tangent points
    # this is the point on the tangent closest to the center.
    # Ax + By + C = 0
    # (xc-x^2) + (yc-y^2) = r^2
    #
    # from C, go a distance r along the normal to the line.
    p1 = c1 - r1*ovd.Point(lines[0].a,lines[0].b)
    p2 = c1 - r1*ovd.Point(lines[1].a,lines[1].b)
    drawCircle(myscreen,p1,0.5,ovdvtk.pink)
    drawCircle(myscreen,p2,0.5,ovdvtk.pink)
    
    p3 = c2 - r2*ovd.Point(lines[0].a,lines[0].b)
    p4 = c2 - r2*ovd.Point(lines[1].a,lines[1].b)
    drawCircle(myscreen,p3,0.5,ovdvtk.pink)
    drawCircle(myscreen,p4,0.5,ovdvtk.pink)
    
    #drawLine( myscreen, l2, ovdvtk.orange )
    #print roots
    
    myscreen.render()
    myscreen.iren.Start()
    
def root_to_line(root,m,n,p,q):
    C = root
    A = m*C+n
    B = p*C+q
    print A, " ", B, " ", C ," ", A*A+B*B
    l = Line(A,B,C,1)
    return l

def drawBitangents2():
    myscreen = ovdvtk.VTKScreen()
    myscreen.camera.SetPosition(0.01, 0,  100 ) 
    myscreen.camera.SetFocalPoint(0, 0, 0)
    myscreen.camera.SetClippingRange(-100,3000)
    
    c1 = ovd.Point(10,20)
    r1=20
    c2 = ovd.Point(6,13)
    r2=23
    # external ma-pocket fails with this input:
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
    [bd1,bd2] = bitanget_direction(c1,c2,r1,r2)
    # from C, go a distance r along the normal to the line.
    p1 = c1 + r1*bd1
    p2 = c1 + r1*bd2
    drawCircle(myscreen,p1,0.5,ovdvtk.pink)
    drawCircle(myscreen,p2,0.5,ovdvtk.red)

    p3 = c2 + r2*bd1
    p4 = c2 + r2*bd2
    drawCircle(myscreen,p3,0.5,ovdvtk.lgreen)
    drawCircle(myscreen,p4,0.5,ovdvtk.green)
    
    myscreen.render()
    myscreen.iren.Start()

def bitanget_direction(c1,c2,r1,r2):
    if r1==r2:
        c1c2 = c2-c1
        c1c2 = c1c2 * (1/c1c2.norm())
        bd1 = c1c2.xy_perp()
        bd2 = -1*c1c2.xy_perp()
        return [bd1, bd2]
    #r_small=0
    #r_large=0
    #c_small=c1
    #c_large=c2
    print " c1 ",c1
    print " c2 ",c2
    """
    if r1>r2:
        r_large = r1
        c_large = c1
        r_small = r2
        c_small = c2
    else:
        r_large = r2
        c_large = c2
        r_small = r1
        c_small = c1
    """
    dr = abs(r1-r2) #r_large-r_small
    # bitangent direction
    # is from c_small to a circle with radius r_large-r_small at c_large
    c1c2_dist = (c1-c2).norm()
    bitang_length = math.sqrt( c1c2_dist*c1c2_dist + dr*dr )
    # now calculate the height of the triangle
    area = dr*bitang_length
    height = area / c1c2_dist
    bitang_c1c2 = math.sqrt( bitang_length*bitang_length - height*height )
    cdir =  c2-c1 #(c_large-c_small)
    if r1>r2:
        cdir = c1-c2
    #print cdir
    cdir = cdir* (1/ cdir.norm()) #(c_large-c_small).normalize()
    #print cdir
    bit1 = bitang_c1c2*cdir + height*cdir.xy_perp()
    bit2 = bitang_c1c2*cdir - height*cdir.xy_perp()
    bit1 = bit1 * (1/bit1.norm())
    bit2 = bit2 * (1/bit2.norm())
    if r1>r2:
        return [bit1.xy_perp(), -1*bit2.xy_perp()]
    else:
        return [-1*bit2.xy_perp(), bit1.xy_perp()]
    #    return [bit1.xy_perp(), bit2.xy_perp()]
    #if r1>r2: 
    
if __name__ == "__main__":  
    #drawPointPointTest()# point-point bisector is a LINE
    
    #drawLinePointTest() # point-line bisector is a PARABOLA
    
    #drawLineLineTest() # line-line bisector is LINE(LINE)
    
    #drawSeparatorSolver1(4.3)
    #drawSeparatorSolver2(4.3)
    #drawBitangents2()
    #drawPointArcTest()
    drawPointArcTest2()
    #drawLineArcTest()
