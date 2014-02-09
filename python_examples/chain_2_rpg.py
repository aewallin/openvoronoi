import openvoronoi as ovd
import ovdvtk
import ovdgenerators as gens
import randompolygon as rpg  # random polygon generator see https://github.com/aewallin/CGAL_RPG

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
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)
    linesegs = 1 # switch to turn on/off line-segments
    
    Npts = 500
    poly = rpg.rpg(Npts, 0)
    # N=3 working seeds: 25-43 45
    
    
    n_problem=Npts # go all the way to the end
    #n_problem=2
    n_step=10
    
    pts=[]
    for p in poly:
        pts.append( ovd.Point( p[0], p[1] ) )
        
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        print m," added vertex "
        m=m+1
    print "polygon is: "
    for idx in id_list:
        print idx,"-",
    print "."
    
    t_after = time.time()
    times.append( t_after-t_before )
    
    print "all point sites inserted. ",
    vd.check()
    #vd.debug_on()

    t_before = time.time()
    if linesegs==1:
        for n in range(len(id_list)):
            n_nxt = n+1
            if n==(len(id_list)-1):
                n_nxt=0
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
            else:
                vd.addLineSite( id_list[n], id_list[n_nxt])
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
    
    err = vd.getStat()
    #print err 
    print "got errorstats for ",len(err)," points"
    if len(err)>1:
        minerr = min(err)
        maxerr = max(err)
        print "min error= ",minerr
        print "max error= ",maxerr
    
    print "num vertices: ",vd.numVertices() 
    print "num SPLIT vertices: ",vd.numSplitVertices() 
        
    calctime = t_after-t_before
    
    vod.setAll()
        
    print "PYTHON All DONE."

    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
     
    myscreen.iren.Start()
