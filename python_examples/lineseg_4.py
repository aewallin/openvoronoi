import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os
import sys

import ovdgenerators as gens

def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd500_zoomout"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    #lwr.Write()


    
if __name__ == "__main__":  
    #print ocl.revision()
    #w=1920
    #h=1080
    w=1024
    h=720
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
    vod.vertexRadius = 0.000031
    vod.drawVertices=0
    vod.drawVertexIndex=0
    
    Nmax = 400
    linesegs = 0
    print "waiting for ",Nmax," segments..",
    sys.stdout.flush()
    segs = gens.randomSegments(far,Nmax)
    print ".done."

    id_list = []
    m=0
    t_before = time.time()
    for seg in segs:
        seg_id=[]
        seg_id.append( vd.addVertexSite( seg[0] ) )
        print m," added vertex", seg_id[0]
        m=m+1
        seg_id.append( vd.addVertexSite( seg[1] ) )
        print m," added vertex", seg_id[1]
        m=m+1
        id_list.append( seg_id )
        #print seg[0].x," , ",seg[1].x

    #exit()
    
    nsegs = Nmax
    #nsegs = 14 #Nmax
    n=1
    for s in id_list:
        if n<= nsegs and linesegs==1:
            vd.addLineSite(s[0],s[1])
            print n," added line-segment"
        n=n+1
    t_after = time.time()
    
    
    #s = id_list[nsegs]
    #vd.addLineSite( s[0], s[1], 5) 
    
    err = vd.getStat()
    #print err 
    print "got errorstats for ",len(err)," points"
    if len(err)>1:
        minerr = min(err)
        maxerr = max(err)
        print "min error= ",minerr
        print "max error= ",maxerr
    
    print "num vertices: ",vd.numVertices() # Nmax=200 gives 1856(187)
    print "num SPLIT vertices: ",vd.numSplitVertices() 
    # nmax= 20 gives 175(13)
    
    # 4 delete-tree
    # 5 create new vertices
    # 6 add startpoint separator
    # 7 add endpoint separator
    
    calctime = t_after-t_before
    if Nmax==0:
        Nmax=1
    print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
    
    vod.setAll()
        
    print "PYTHON All DONE."

    myscreen.render()    
    myscreen.iren.Start()
