import openvoronoi as ovd
import ovdvtk
import math
import vtk
import time

pts = []
# calculate Hilbert curve, using recursion
# based on: http://www.andyshelley.co.uk/axishilbert/index.php
#    Copyright 2010   Andy Shelley <andy@andyshelley.co.uk>
#    GPL v2 or later
# self.Hilbert    (0, self.Blocks, 90, self.Side, self.Feed)
def Hilbert ( ca, level, angle, size, p): 
    if level == 0:
        return

    ca = (ca + 360 - angle) % 360
    Hilbert( ca, level - 1, -angle, size,p )
    step = size
    if (ca == 0):
        p[0]=p[0]+step
    elif (ca == 90):
        p[1]=p[1]-step
    elif (ca == 180):
        p[0]=p[0]-step
    else:
        p[1]=p[1]+step
    pts.append( [ p[0], p[1] ] )

    ca = (ca + angle) % 360
    Hilbert ( ca, level - 1, angle, size, p)
    if (ca == 0):
        p[0]=p[0]+step
    elif (ca == 90):
        p[1]=p[1]-step
    elif (ca == 180):
        p[0]=p[0]-step
    else:
        p[1]=p[1]+step
    pts.append( [ p[0], p[1] ] )
    Hilbert ( ca, level - 1, angle, size, p)
    
    ca = (ca + angle) % 360
    if (ca == 0):
        p[0]=p[0]+step 
    elif (ca == 90):
        p[1]=p[1]-step
    elif (ca == 180):
        p[0]=p[0]-step
    else:
        p[1]=p[1]+step
    pts.append( [ p[0], p[1] ] )
    Hilbert ( ca, level - 1, -angle, size, p)
    ca = (ca + 360 - angle) % 360

def insert_curve(vd,curve):
    pts=[]
    for p in curve:
        pts.append( ovd.Point( p[0], p[1] ) )
    id_list = []
    print "inserting ",len(pts)," point-sites...",
    m=0
    t_before = time.time()
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        print " ",m," added vertex ", id_list[ len(id_list) -1 ]
        m=m+1
    t_after = time.time()
    pt_time = t_after-t_before
    print "done."
    j=0
    t_before = time.time()
    print "inserting ",len(id_list)," line-segments..."
    for n in range( len(id_list)-1 ):
        n_nxt = n+1
        #if n==(len(id_list)-1):
        #    n_nxt=0
        print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
        vd.addLineSite( id_list[n], id_list[n_nxt])
        j=j+1
    t_after = time.time()
    seg_time = t_after-t_before
    print "done."
    return [pt_time, seg_time]
    
if __name__ == "__main__":
    
    w=1024
    h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    #ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    # draw a unit-circle
    ca = ovdvtk.Circle(center=(0,0,0) , radius=1, color=(0,1,1), resolution=50 )
    myscreen.addActor(ca)   
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )  
    
    
    scale=1
    far = 1
    camPos = far
    zmult = 3
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    size = 1.1
    
    # bounding-box
    """
    myscreen.addActor( ovdvtk.Line(p1=(-size/2,-size/2,0),p2=(-size/2,size/2,0),color=ovdvtk.orange) )
    myscreen.addActor( ovdvtk.Line(p1=(-size/2,size/2,0),p2=(size/2,size/2,0),color=ovdvtk.orange) )
    myscreen.addActor( ovdvtk.Line(p1=(size/2,size/2,0),p2=(size/2,-size/2,0),color=ovdvtk.orange) )
    myscreen.addActor( ovdvtk.Line(p1=(size/2,-size/2,0),p2=(-size/2,-size/2,0),color=ovdvtk.orange) )
    """
    
    startpt = [-size/2,-size/2]
    level = 2
    step = size / math.pow(2, level)

    Hilbert( 0, level, 90, step, startpt)
    print "Generated hilbert curve, level= ", level
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.03)
    vod.drawFarCircle()

    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=1
    vod.drawGenerators=1
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)
    
    times = insert_curve(vd,pts)
    
    vd.check()
    #print pts
    #drawCurve2(myscreen,pts,ovdvtk.yellow)
    vod.setVDText2(times)
    vod.setAll()
    print "PYTHON All DONE."

    myscreen.render()        
    myscreen.iren.Start()
