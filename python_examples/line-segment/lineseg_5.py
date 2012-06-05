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
import ovdgenerators as gens

def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd500_zoomout"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    #lwr.Write()

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
    vod.drawVertices=1
    vod.drawVertexIndex=1
    vod.drawGenerators=1
    
    
    linesegs = 1 # switch to turn on/off line-segments
    
    segs = []
    s1 = [ovd.Point(-0.5,-0.7) , ovd.Point(-0.5,+0.7) ]
    s2 = [ovd.Point(+0.5,-0.7) , ovd.Point(+0.5,+0.7) ]
    s3 = [ovd.Point(+0.4,-0.7) , ovd.Point(+0.4,+0.7) ]
    s4 = [ovd.Point(-0.3,-0.7) , ovd.Point(-0.3,+0.7) ]
    s5 = [ovd.Point(+0.3,-0.7) , ovd.Point(+0.3,+0.7) ]
    
    #s3 = [ovd.Point(+0.4,+0.0) , ovd.Point(-0.4,+0.0) ]
    segs.append(s1)
    segs.append(s2)    
    segs.append(s3)
    segs.append(s4)
    segs.append(s5)
    
    #t_after = time.time()
    #print ".done in {0:.3f} s.".format( t_after-t_before )
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for seg in segs:
        seg_id=[]
        seg_id.append( vd.addVertexSite( seg[0] ) )
        #print m," added vertex", seg_id[0]
        m=m+1
        seg_id.append( vd.addVertexSite( seg[1] ) )
        #print m," added vertex", seg_id[1]
        m=m+1
        id_list.append( seg_id )
        #print seg[0].x," , ",seg[1].x
    t_after = time.time()
    times.append( t_after-t_before )
    #exit()
    
    #print "   ",2*Nmax," point-sites sites took {0:.3f}".format(times[0])," seconds, {0:.2f}".format( 1e6*float( times[0] )/(float(2*Nmax)*float(math.log10(2*Nmax))) ) ,"us/n*log(n)"
    print "all point sites inserted. ",
    vd.check()
    
    #nsegs = Nmax
    nsegs = 5 #Nmax
    n=1
    t_before = time.time()
    for s in id_list:
        if n<= nsegs and linesegs==1:
            vd.addLineSite(s[0],s[1])
            print n," added line-segment"
        n=n+1
    t_after = time.time()
    line_time = t_after-t_before
    if line_time < 1e-3:
        line_time = 1
    times.append( line_time )
    
    vd.check()
    
    #s = id_list[nsegs]
    #vd.debug_on()
    #vd.addLineSite( s[0], s[1], 10) 
    #seg = id_list[nsegs]
    #vd.addLineSite(seg[0],seg[1],10)
    # 4 delete-tree
    # 5 create new vertices
    # 6 add startpoint separator
    # 7 add endpoint separator
    # 8 add new edges
    # 9 delete delete-tree edges
    # 10 reset status
            
    #vod.setVDText2(times)
    
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
    

    
    calctime = t_after-t_before
    #if Nmax==0:
    #    Nmax=1
    #print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
    
    vod.setAll()
        
    print "PYTHON All DONE."

    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
     
    myscreen.iren.Start()
