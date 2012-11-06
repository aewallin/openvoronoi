import truetypetracer as ttt                 # https://github.com/aewallin/truetype-tracer
import openvoronoi as ovd  # https://github.com/aewallin/openvoronoi
import ovdvtk
import time
import vtk

def drawLoops(myscreen,loops,loopColor):
    # draw the loops
    nloop = 0
    for lop in loops:
        n = 0
        N = len(lop)
        first_point=[]
        previous=[]
        for p in lop:
            if n==0: # don't draw anything on the first iteration
                previous=p 
                first_point = p
            elif n== (N-1): # the last point
                myscreen.addActor( ovdvtk.Line(p1=(previous[0],previous[1],0),p2=(p[0],p[1],0),color=loopColor) ) # the normal line
                # and a line from p to the first point
                myscreen.addActor( ovdvtk.Line(p1=(p[0],p[1],0),p2=(first_point[0],first_point[1],0),color=loopColor) )
            else:
                myscreen.addActor( ovdvtk.Line(p1=(previous[0],previous[1],0),p2=(p[0],p[1],0),color=loopColor) )
                previous=p
            n=n+1
        print "rendered loop ",nloop, " with ", len(lop), " points"
        nloop = nloop+1

def translate(segs,x,y):
    out = []
    for seg in segs:
        seg2 = []
        for p in seg:
            p2 = []
            p2.append(p[0] + x)
            p2.append(p[1] + y)
            seg2.append(p2)
            #seg2.append(seg[3] + y)
        out.append(seg2)
    return out

def modify_segments(segs):
    segs_mod =[]
    for seg in segs:
        first = seg[0]
        last = seg[ len(seg)-1 ]
        assert( first[0]==last[0] and first[1]==last[1] )
        seg.pop()
        segs_mod.append(seg)
    return segs_mod

def draw_ttt(myscreen, text, x,y,scale):
    wr = ttt.SEG_Writer()
    wr.arc = False
    wr.conic = False
    wr.cubic = False
    wr.scale = float(1)/float(scale)
    # "L" has 36 points by default
    wr.conic_biarc_subdivision = 10 # this has no effect?
    wr.conic_line_subdivision = 10 # this increasesn nr of points to 366
    wr.cubic_biarc_subdivision = 10 # no effect?
    wr.cubic_line_subdivision = 10 # no effect?
    wr.setFont(2)
    s3 = ttt.ttt(text,wr) 
    ext = wr.extents
    dx = ext.maxx-ext.minx

    segs = wr.get_segments()
    segs = translate(segs, x, y)
    print "number of polygons: ", len(segs)
    np = 0
    for s in segs:
        print " polygon ",np," has ",len(s)," points"
        np=np+1        
    segs = modify_segments(segs)
    drawLoops(myscreen, segs, ovdvtk.yellow )
    
# this script only draws geometry from ttt
# no voronoi-diagram is created!
if __name__ == "__main__": 
    print "ttt version = ",ttt.version()
    #w=2500
    #h=1500
    
    #w=1920
    #h=1080
    #w=1024
    #h=1024
    w=800
    h=600

    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
        
    scale=1
    far = 1
    camPos = far
    zmult = 3
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    # draw a unit-circle
    ca = ovdvtk.Circle(center=(0,0,0) , radius=1, color=(0,1,1), resolution=50 )
    myscreen.addActor(ca)   
    draw_ttt(myscreen, "A", 0,0,10000)
    #draw_ttt(myscreen, "ABCDEFGHIJKLMNOPQRSTUVWXYZ", -0.5,0,80000)
    #draw_ttt(myscreen, "abcdefghijklmnopqrstuvwxyz", -0.5,-0.1,80000)
    #draw_ttt(myscreen, "1234567890*", -0.5,-0.2,80000)
    #draw_ttt(myscreen, "m", -0.5,-0.2,80000)
    print "PYTHON All DONE."

    myscreen.render()        
    myscreen.iren.Start()
