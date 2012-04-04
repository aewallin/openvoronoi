import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os

def drawLine(myscreen, p1, p2):
    myscreen.addActor( ovdvtk.Line( p1 = (p1.x,p1.y,0), p2=(p2.x,p2.y,0), color = ovdvtk.yellow ) )

def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd500_zoomout"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    #lwr.Write()

if __name__ == "__main__":  
    #print ocl.revision()
    myscreen = ovdvtk.VTKScreen(width=1024, height=720) #(width=1920, height=1080)
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #w2if.Modified()
    #lwr.SetFileName("tux1.png")
    
    scale=1
    myscreen.render()
    random.seed(2)
    far = 1
    camPos = far
    zmult = 4
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    vd.check()
    print "created."
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=1
    vod.drawVertexIndex=1
    vod.drawGenerators=1
    
    vod.offsetEdges = 1
    vd.setEdgeOffset(0.05)
    
    Nmax = 6
    
    #plist = randomGenerators(far, Nmax)
    #plist = regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    
    #plist = randomGenerators(far, Nmax) 
    plist = []
    plist.append( ovd.Point(0.1,0.1) )
    plist.append( ovd.Point(-0.1,0.1) )
    #plist.append( ovd.Point(-0.15, -0.15) )
    #+ regularGridGenerators(far, Nmax) + circleGenerators(far, Nmax)

    #plist = [ovd.Point(0,0)]
    print plist
    times=[]
    t_before = time.time() 
    n=0
    id_list=[]
    #vd.debug_on()
    for p in plist: 
        print n," adding ",p
        id_list.append( vd.addVertexSite( p ) )
        n=n+1
    t_after = time.time()
    calctime = t_after-t_before
    times.append(calctime)
    
    id1 = id_list[0]
    id2 = id_list[1]
    #id3 = id_list[2]
    #id4 = id_list[3]
    #print "add segment ",id1, " to ", id2
    vd.debug_on()
    c1 = ovd.Point(0,0)
    vd.addArcSite( id1, id2 , c1, True , 13 )
    
    #vd.addLineSite( id3, id4 )
    t_after = time.time()
    calctime = t_after-t_before
    times.append(calctime)
    if Nmax==0:
        Nmax=1
    print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
    
    vod.setVDText2(times)
    
    vod.setAll()
    myscreen.render()
            

        
    print "PYTHON All DONE."

    myscreen.render()    
    myscreen.iren.Start()
