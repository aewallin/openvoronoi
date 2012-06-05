import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os
import gc



def drawLine(myscreen, p1, p2):
    myscreen.addActor( ovdvtk.Line( p1 = (p1.x,p1.y,0), p2=(p2.x,p2.y,0), color = ovdvtk.yellow ) )

def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd_lin50"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    lwr.Write()

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
    

def drawFrame(Nmax, myscreen, vd, vod, framenr, anim, step):
    #Nmax = 2

    segs = []
    id_list = []
    for n in range(Nmax):
        seg = randomGenerators()
        while segmentIntersects(segs, seg):
            seg = randomGenerators()
        segs.append(seg)
    anim_count=0
    for seg in segs:
        seg_id=[]
        if (anim_count == anim):
            r1 = vd.addVertexSite( seg[0] , step )
            vod.setAll()
            myscreen.render()    
            writeFrame( w2if, lwr, framenr )
            gc.collect()
            return
        else:
            seg_id.append( vd.addVertexSite( seg[0] ) )
            anim_count = anim_count+1
        
        if (anim_count == anim):
            vd.addVertexSite( seg[1] , step )
            vod.setAll()
            myscreen.render()    
            writeFrame( w2if, lwr, framenr )
            gc.collect()
            return
        else:
            seg_id.append( vd.addVertexSite( seg[1] ) )
            anim_count = anim_count+1
            
        #seg_id.append( vd.addVertexSite( seg[1] ) )
        id_list.append( seg_id )

    for s in id_list:
        vd.addLineSite(s[0],s[1])

    vod.setAll()
        
    print "PYTHON All DONE."

    myscreen.render()    
    #myscreen.iren.Start()

def drawLinesegFrame(Nmax, myscreen, vd, vod, framenr, anim, step):
    segs = []
    id_list = []
    
    for n in range(Nmax):
        seg = randomGenerators()
        while segmentIntersects(segs, seg):
            seg = randomGenerators()
        segs.append(seg)
    anim_count=0
    for seg in segs:
        seg_id=[]
        seg_id.append( vd.addVertexSite( seg[0] ) )
        seg_id.append( vd.addVertexSite( seg[1] ) )
        id_list.append( seg_id )

    for s in id_list:
        if (anim_count==anim):
            vd.addLineSite(s[0],s[1],step)
            vod.setAll()
            myscreen.render()    
            writeFrame( w2if, lwr, framenr )
            gc.collect()
            return
        else:
            vd.addLineSite(s[0],s[1])
            anim_count=anim_count+1

        
    print "PYTHON All DONE."

    myscreen.render()    

if __name__ == "__main__":
    Nmax = 50 # number of points to insert
    
    nframe=0 
    w=1920
    h=1080
    #w=1024
    #h=720
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    
    for npt in range(Nmax*2): # animate insertion of point-sites
        for step in [x+1 for x in range(6)]:
            scale=1
            myscreen.render()
            random.seed(42)
            far = 1
            camPos = far
            zmult = 4
            # camPos/float(1000)
            myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
            myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
            myscreen.camera.SetFocalPoint(0, 0, 0)

            vd = ovd.VoronoiDiagram(far,120)

            vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
            vod.drawFarCircle()
            vod.vertexRadius = 0.005
            vod.textScale = 0.01
            vod.drawVertexIndex=0
            drawFrame(Nmax,myscreen, vd, vod, nframe,npt,step)
            vd.reset_vertex_count()
            gc.collect()
            nframe = nframe+1
            # remove all actors
            acts = vod.getActors()
            for a in acts:
                myscreen.removeActor(a)
                
    # now animate insertion of line-segments
    for npt in range(Nmax): 
        for step in [x+1 for x in range(10)]:
            scale=1
            myscreen.render()
            random.seed(42)
            far = 1
            camPos = far
            zmult = 4
            # camPos/float(1000)
            myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
            myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
            myscreen.camera.SetFocalPoint(0, 0, 0)

            vd = ovd.VoronoiDiagram(far,120)

            vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
            vod.drawFarCircle()
            vod.vertexRadius = 0.005
            vod.textScale = 0.02
            vod.drawVertexIndex=0
            drawLinesegFrame(Nmax,myscreen, vd, vod, nframe, npt,step)
            vd.reset_vertex_count()
            gc.collect()
            nframe = nframe+1
            # remove all actors
            acts = vod.getActors()
            for a in acts:
                myscreen.removeActor(a)
