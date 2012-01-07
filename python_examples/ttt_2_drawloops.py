import ttt
import openvoronoi as ovd
import ovdvtk
import time
import vtk
"""
def drawSegment(myscreen, seg):
    #p1x = seg[0]
    #p1y = seg[1]
    #p2x = seg[2]
    #p2y = seg[3]
    for pt in seg:
    actor = ovdvtk.Line( p1=( p1x,p1y, 0), p2=(p2x,p2y, 0), color=ovdvtk.yellow)
    myscreen.addActor(actor)
"""
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
    
if __name__ == "__main__":  
    #w=2500
    #h=1500
    
    #w=1920
    #h=1080
    w=1024
    h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.revision() )
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #w2if.Modified()
    #lwr.SetFileName("tux1.png")
    
    scale=1
    myscreen.render()
    #random.seed(42)
    far = 1
    camPos = far
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    # draw a unit-circle
    ca = ovdvtk.Circle(center=(0,0,0) , radius=1, color=(0,1,1), resolution=50 )
    myscreen.addActor(ca)   
    
    
    wr = ttt.SEG_Writer()

    # wr.scale = 3
    wr.arc = False
    wr.conic = False
    wr.cubic = False
    wr.scale = float(1)/float(80000)
    s3 = ttt.ttt("ABCDEFGHIJKLMNOPQRSTUVWXYZ",wr) 
    ext = wr.extents
    print ext
    dx = ext[1]-ext[0]
    
    segs = wr.get_segments()
    segs = translate(segs, -0.5, -0.4)

    #print s
    #print segs
    #exit()
    print "number of polygons: ", len(segs)
    np = 0
    for s in segs:
        print " polygon ",np," has ",len(s)," points"
        np=np+1
        
    segs_mod =[]
    for seg in segs:
        first = seg[0]
        last = seg[ len(seg)-1 ]
        assert( first[0]==last[0] and first[1]==last[1] )
        seg.pop()
        segs_mod.append(seg)
        #drawSegment(myscreen, seg)
    segs = segs_mod
    
    drawLoops(myscreen, segs, ovdvtk.yellow )
    
    
    print "PYTHON All DONE."

    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
     
    myscreen.iren.Start()
