import openvoronoi as ovd
import ovdvtk
import time
import vtk
import math

# draw a line from pt1 to pt2 with the given color
def drawLine(myscreen, pt1, pt2, lineColor):
    myscreen.addActor( ovdvtk.Line(p1=(pt1.x,pt1.y,0),p2=(pt2.x,pt2.y,0),color=lineColor) ) 

# rotate a point in 2D by cos/sin. from emc2 gcodemodule.cc
def rotate(x, y,  c,  s):
    tx = x * c - y * s;
    y = x * s + y * c;
    x = tx;
    return [x,y]

# draw an arc using many line-segments
# start at pt1, and at pt2, center at cen
# radius r, cw is a bool cw/ccw flag
def drawArc(myscreen, pt1, pt2, r, cen,cw,arcColor):
    start = pt1-cen
    end = pt2-cen
    theta1 = math.atan2(start.x,start.y)
    theta2 = math.atan2(end.x,end.y)
    alfa=[] # the list of angles
    da=0.1
    CIRCLE_FUZZ = 1e-9
    # idea from emc2 / cutsim g-code interp G2/G3
    if (cw == False ): 
        while ( (theta2 - theta1) > -CIRCLE_FUZZ): 
            theta2 -= 2*math.pi
    else:
        while( (theta2 - theta1) < CIRCLE_FUZZ): 
            theta2 += 2*math.pi
    
    dtheta = theta2-theta1
    arclength = r*dtheta
    dlength = min(0.01, arclength/10)
    steps = int( float(arclength) / float(dlength))
    rsteps = float(1)/float(steps)
    dc = math.cos(-dtheta*rsteps) # delta-cos  
    ds = math.sin(-dtheta*rsteps) # delta-sin
    
    previous = pt1
    tr = [start.x, start.y]
    for i in range(steps):
        tr = rotate(tr[0], tr[1], dc, ds) #; // rotate center-start vector by a small amount
        x = cen.x + tr[0] 
        y = cen.y + tr[1] 
        current = ovd.Point(x,y)
        myscreen.addActor( ovdvtk.Line(p1=(previous.x,previous.y,0),p2=(current.x,current.y,0),color=arcColor) )
        previous = current 

# draw many offsets
def drawOffsets(myscreen, ofs):
    # draw loops
    nloop = 0
    lineColor = ovdvtk.green
    arcColor = ovdvtk.grass
    for lop in ofs:
        n = 0
        N = len(lop)
        first_point=[]
        previous=[]
        for p in lop:
            # p[0] is the Point
            # p[1] is -1 for lines, and r for arcs
            if n==0: # don't draw anything on the first iteration
                previous=p[0]
            else:
                cw=p[3]
                cen=p[2]
                r=p[1]
                p=p[0]
                if r==-1: # this offset element is a line
                    drawLine(myscreen, previous, p, lineColor)
                else: # this offset element is an arc
                    drawArc(myscreen, previous, p, r,cen,cw, arcColor)
                previous=p
            n=n+1
        print "rendered loop ",nloop, " with ", len(lop), " points"
        nloop = nloop+1

if __name__ == "__main__":  
    #w=2500 # screen resolution for big screens
    #h=1500
    
    #w=1920
    #h=1080
    w=1024
    h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) # a VTK window for drawing 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )   # the OpenVoronoi text, revision, and date
    
    scale=1
    myscreen.render()

    far = 1
    camPos = far
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    # these actions on the vod object control how the VD is drawn using VTK
    vod.drawFarCircle()
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=1
    vod.drawGenerators=1
    vod.offsetEdges = 0 # for debug. a bool flag to set null-edge drawing on/off. use together with setEdgeOffset()
    vd.setEdgeOffset(0.05) # for debug. a non-zero value will draw null-edges as circular arcs
    # null-edges are an internal openvoronoi construction to avoid high-degree vertices in the VD-graph
    # they are not relevant for upstream or downstream algorithms
    
    # input points (vertices/sites)
    p1=ovd.Point(-0.1,-0.2)
    p2=ovd.Point(0.2,0.1)
    p3=ovd.Point(0.4,0.2)
    p4=ovd.Point(0.6,0.6)
    p5=ovd.Point(-0.6,0.3)
    pts = [p1,p2,p3,p4,p5] # a list of all points in the input
    
    #t_after = time.time()
    #print ".done in {0:.3f} s.".format( t_after-t_before )
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts: # add all points before adding line-segments
        id_list.append( vd.addVertexSite( p ) )
        #print m," added vertex", seg_id[0]
        m=m+1

    t_after = time.time()
    times.append( t_after-t_before )
    print "all point sites inserted. "
    print "VD check: ", vd.check()
    
    t_before = time.time()
    # now add line-segments, by using the integer indexes returned by vd.addVertexSite() above
    vd.addLineSite( id_list[0], id_list[1])
    vd.addLineSite( id_list[1], id_list[2])
    vd.addLineSite( id_list[2], id_list[3])
    vd.addLineSite( id_list[3], id_list[4])
    vd.addLineSite( id_list[4], id_list[0])
    vd.check()
    t_after = time.time()
    line_time = t_after-t_before
    if line_time < 1e-3:
        line_time = 1
    times.append( line_time )
    
    # we now have a VD of the input sites
    # we can now run downstream algorithms on the VD such as
    # Offset, MedialAxis, etc.
    
    of = ovd.Offset( vd.getGraph() ) # pass the created graph to the Offset class
    of.str() # text output, for debug
    ofs = of.offset(0.123) # generate offsets at the given distance.
    drawOffsets(myscreen, ofs) # draw the generated offsets

    vod.setVDText2(times)
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
