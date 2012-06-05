import openvoronoi as ovd
import ovdvtk

import time
import vtk
import datetime
import math
import random
import os
import sys
import pickle
import gzip

def drawCircle(myscreen, c, r, circlecolor):
    ca = ovdvtk.Circle(center=(c.x,c.y,0) , radius=r, color=circlecolor, resolution=50 )
    myscreen.addActor(ca)
    
if __name__ == "__main__":  
    #w=2500
    #h=1500
    
    w=1920
    h=1080
    #w=1024
    #h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    
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
    zmult = 1.8
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0.22, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()

    
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)
    
    
    linesegs = 1 # switch to turn on/off line-segments
    
    segs = []
    #ovd.Point(1,1)
    eps=0.9
    p1=ovd.Point(-0.1,-0.2)
    p2=ovd.Point(0.2,0.1)
    p3=ovd.Point(0.4,0.2)
    p4=ovd.Point(0.6,0.6)
    p5=ovd.Point(-0.6,0.3)

    pts = [p1,p2,p3,p4,p5]
    
    #t_after = time.time()
    #print ".done in {0:.3f} s.".format( t_after-t_before )
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts:
        
        id_list.append( vd.addVertexSite( p ) )
        #print m," added vertex", seg_id[0]
        m=m+1
   
    t_after = time.time()
    times.append( t_after-t_before )
    
    print "all point sites inserted. "
    vd.check()
    
    t_before = time.time()
    vd.addLineSite( id_list[0], id_list[1])
    vd.addLineSite( id_list[1], id_list[2])
    vd.addLineSite( id_list[2], id_list[3])
    vd.addLineSite( id_list[3], id_list[4])
    vd.addLineSite( id_list[4], id_list[0])
    t_after = time.time()
    times.append( t_after-t_before )
    vd.check()
    
    pi = ovd.PolygonInterior(True)
    vd.filter_graph(pi)
    ma = ovd.MedialAxis()
    vd.filter_graph(ma)
    
    vod.setVDText2(times)
    
    vod.setAll()
    
    mapocket = ovd.MedialAxisPocket(vd.getGraph())
    mapocket.setWidth(0.005)
    
    
    
    mapocket.run()

    mic_components = mapocket.get_mic_components()
    for mic_list in mic_components:
        nframe=0
        
        for n in range( len(mic_list) ):
            mic = mic_list[n]
            if n == 0:
                print "hello", mic[0]," r = ",mic[1]
                drawCircle( myscreen, mic[0], mic[1] , ovdvtk.red )
            else:
                drawCircle( myscreen, mic[0], mic[1] , ovdvtk.green )
            w2if.Modified()
            lwr.SetFileName("frames/%06d.png" % ( nframe ) )
            #lwr.Write()
            time.sleep(0.1)
            myscreen.render()
        print "mic done."
    
    print "PYTHON All DONE."

    myscreen.render()   

     
    myscreen.iren.Start()
