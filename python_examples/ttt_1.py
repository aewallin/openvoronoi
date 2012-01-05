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

def insert_polygon(vd, polygon):
    pts=[]
    for p in polygon:
        pts.append( ovd.Point( p[0], p[1] ) )
        
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        print m," added vertex "
        m=m+1
    t_after = time.time()
    times.append( t_after-t_before )
    
    print "polygon is: "
    for idx in id_list:
        print idx," ",
    print "."
    

    
    print "all point sites inserted. ",
    vd.check()
    #vd.debug_on()

    t_before = time.time()

    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        
        vd.addLineSite( id_list[n], id_list[n_nxt])
        """
        if n == n_problem:
            vd.addLineSite( id_list[n], id_list[n_nxt],n_step)
            
            t_after = time.time()
            line_time = t_after-t_before
            times.append( max(line_time,0.1) )
            vod.setVDText2(times)
            vod.setAll()
            myscreen.render()
            myscreen.iren.Start()
            raw_input("Press ENTER to exit")
        """
    t_after = time.time()
    times.append( t_after-t_before )
    
    vd.check()
    return times
    
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
    
    
    wr = ttt.SEG_Writer()

    # wr.scale = 3
    wr.arc = False
    wr.conic = False
    wr.cubic = False
    wr.scale = float(1)/float(4000)
    s3 = ttt.ttt("L",wr)
    #print s3
    segs = wr.get_segments()
    segs = translate(segs, -0.5, -0.4)
    #print s
    print segs
    #exit()
    segs_mod =[]
    for seg in segs:
        first = seg[0]
        last = seg[ len(seg)-1 ]
        assert( first[0]==last[0] and first[1]==last[1] )
        seg.pop()
        segs_mod.append(seg)
        #drawSegment(myscreen, seg)
    segs = segs_mod
    
    #drawLoops(myscreen, segs, ovdvtk.yellow )
    
    vd = ovd.VoronoiDiagram(far,120)
    print vd.version()
    
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)
    
    times = [0,0]
    for polygon in segs:
        tp = insert_polygon(vd, polygon )
        times[0] = times[0] + tp[0]
        times[1] = times[1] + tp[1]
        
    vod.setVDText2(times)
    vod.setAll()
    
    print "PYTHON All DONE."

    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
     
    myscreen.iren.Start()
