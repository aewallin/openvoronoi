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

if __name__ == "__main__":  
    #w=2500
    #h=1500
    
    #w=1920
    #h=1080
    w=1024
    h=1024
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
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()

    
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=1
    vod.drawGenerators=1
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
    #exit()
    
    #print "   ",2*Nmax," point-sites sites took {0:.3f}".format(times[0])," seconds, {0:.2f}".format( 1e6*float( times[0] )/(float(2*Nmax)*float(math.log10(2*Nmax))) ) ,"us/n*log(n)"
    print "all point sites inserted. ",
    vd.check()
    
    #nsegs = Nmax
    #nsegs = 5 #Nmax
    #n=1
    t_before = time.time()
    
    #vd.debug_on()
    vd.addLineSite( id_list[0], id_list[1])
    
    # 4 augment vertex set
    # 5 process null-faces
    # 6 create faces and pseudo-edges
    # 7 add new vertices
    # 8,9,10,11 separators
    # 11, all separators added
    # 12, add edges
    # 13, remove delete-set
    # 13, reset status, delete split-verts
    
    vd.check()
    
    #vd.debug_on()
    vd.addLineSite( id_list[1], id_list[2])
    vd.check()

        
    vd.addLineSite( id_list[2], id_list[3])
    vd.check()
    
    #vd.debug_on()
    
    vd.addLineSite( id_list[3], id_list[4])
    vd.check()
    
    vd.addLineSite( id_list[4], id_list[0])
    vd.check()
    
    t_after = time.time()
    line_time = t_after-t_before
    if line_time < 1e-3:
        line_time = 1
    times.append( line_time )
    
    #s = id_list[nsegs]
    #vd.debug_on()
    #vd.addLineSite( s[0], s[1], 10) 
    #seg = id_list[nsegs]
    #vd.addLineSite(seg[0],seg[1],10)
    # 1 identify start/endvert
    # 2 add line-segment edges/sites to graph
    # 3 identify seed-vertex
    # 4 create delete-tree
    # 5 create new vertices
    # 6 add startpoint pos separator
    # 7 add startoiubt neg separator
    # 8 add end-point pos separator
    # 9 add end-point neg separator
    # 10 add new edges
    # 11 delete delete-tree edges
    # 12 reset status
            
    vod.setVDText2(times)
    
    #err = vd.getStat()
    #print err 
    """
    print "got errorstats for ",len(err)," points"
    if len(err)>1:
        minerr = min(err)
        maxerr = max(err)
        print "min error= ",minerr
        print "max error= ",maxerr
    
    print "num vertices: ",vd.numVertices() 
    print "num SPLIT vertices: ",vd.numSplitVertices() 
    """
    #calctime = t_after-t_before
    
    #pi = ovd.PolygonInterior(False)
    #vd.filter_graph(pi)
    #pi = ovd.IslandFilter()
    #vd.filter_graph(pi)
    
    vod.setAll()
        
    print "PYTHON All DONE."

    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
     
    myscreen.iren.Start()
