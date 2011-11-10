import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os

def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd500_zoomout"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    #lwr.Write()

def randomGenerators(far=1, Nmax=2):
    pradius = (1.0/math.sqrt(2))*far
    plist=[]
    for n in range(Nmax):
        x=-pradius+2*pradius*random.random()
        y=-pradius+2*pradius*random.random()
        plist.append( ovd.Point(x,y) )
    return plist

def intersects(s1,s2):
    p1 = s1[0]
    p2 = s1[1]
    p = p1
    r = p2-p1
    q1 = s2[0]
    q2 = s2[1]
    q = q1
    s = q2-q1
    # t = (q-p) cross (s) / (r cross s)
    # u = (q-p) cross (r) / (r cross s)
    if ( r.cross(s) == 0 ): #parallel lines
        if ( (q-p).cross(r) == 0 ): #collinear
            return 1   
        else:
            return 0 # parallel lines that never intersect

    t = (q-p).cross(s) / (r.cross(s))
    u = (q-p).cross(r) / (r.cross(s))
    if ( (0<=t) and (t<=1) and (0<=u) and (u<=1) ):
        return 1
    return 0

# test if s intersects with any of the segments in seg
def segmentIntersects(segs, s):
    for sg in segs:
        if intersects(sg,s):
            return 1
    
    return 0 # no intersections found
    
if __name__ == "__main__":  
    #print ocl.revision()
    myscreen = ovdvtk.VTKScreen(width=1024, height=720) #(width=1920, height=1080)
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.revision() )
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #w2if.Modified()
    #lwr.SetFileName("tux1.png")
    
    scale=1
    myscreen.render()
    random.seed(42)
    far = 1
    camPos = far
    zmult = 4
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print vd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    vod.textScale = 0.002
    vod.vertexRadius = 0.001
    
    t_before = time.time()
    
    Nmax = 19

    segs = []
    id_list = []
    for n in range(Nmax):
        seg = randomGenerators()
        while segmentIntersects(segs, seg):
            seg = randomGenerators()
        segs.append(seg)
        seg_id=[]
        seg_id.append( vd.addVertexSite( seg[0] ) )
        seg_id.append( vd.addVertexSite( seg[1] ) )
        id_list.append( seg_id )
        #print seg[0].x," , ",seg[1].x

    #exit()
    
    nsegs =18 # Nmax
    n=1
    for s in id_list:
        if n<= nsegs:
            vd.addLineSite(s[0],s[1])
        n=n+1
    s = id_list[18]
    vd.addLineSiteStep( s[1], s[0], 12) 
    t_after = time.time()
    calctime = t_after-t_before
    if Nmax==0:
        Nmax=1
    print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
    
    vod.setAll()
        
    print "PYTHON All DONE."

    myscreen.render()    
    myscreen.iren.Start()
