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
    #print "created."
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    
    vod.textScale = 0.0005
    vod.vertexRadius = 0.00031
    vod.drawVertices=1
    vod.drawVertexIndex=1
    vod.drawGenerators=1
    
    vod.offsetEdges = 1
    vd.setEdgeOffset(0.0053)

    plist = []
    d=0.02
    plist.append( ovd.Point( -0.1+d, 0.1) ) # 0
    plist.append( ovd.Point( 0.1-d,  0.1) ) # 1
    plist.append( ovd.Point( 0.1, 0.1-d) )  # 2
    plist.append( ovd.Point( 0.1,-0.1+d) )  # 3
    plist.append( ovd.Point( 0.1-d,-0.1) )  # 4
    plist.append( ovd.Point( -0.1+d,-0.1) ) # 5
    plist.append( ovd.Point( -0.1,-0.1+d) ) # 6
    plist.append( ovd.Point( -0.1, 0.1-d) ) # 7
    
    #plist.append( ovd.Point(-0.03,-0.03) )
    #plist.append( ovd.Point(-0.15, -0.15) )
    #+ regularGridGenerators(far, Nmax) + circleGenerators(far, Nmax)

    #plist = [ovd.Point(0,0)]
    #print plist
    times=[]
    t_before = time.time() 
    n=0
    id_list=[]
    #vd.debug_on()
    for p in plist: 
        print n," adding PointSite ",p
        id_list.append( vd.addVertexSite( p ) )
        n=n+1
    t_after = time.time()
    calctime = t_after-t_before
    times.append(calctime)
    
    
    vd.addLineSite(id_list[0],id_list[1])
    vd.addLineSite(id_list[2],id_list[3])
    vd.addLineSite(id_list[4],id_list[5])
    vd.addLineSite(id_list[6],id_list[7])
    #vd.addLineSite(id4,id3)
    #vd.addLineSite(id4,id2)
    
    #print "add segment ",id1, " to ", id2
    vd.debug_on()
    c1 = ovd.Point(0.1-d,0.1-d)

    vd.addArcSite( id_list[1], id_list[2] , c1, False  ) 
    #vd.addArcSite( id_list[2], id_list[1] , c1, True ) 

    #c2 = ovd.Point(0.1-d,-0.1-d)
    #vd.addArcSite( id_list[3], id_list[4] , c2, False  ) 
    
    #vd.check()
    
    
    #vd.addArcSite( id2, id1 , c1, False , 14) # left-to-right arc
    
    # 1 lookup vertex descriptors
    # 2 create Sites
    # 3 find seed vertex
    # 4 find delete-tree
    # 5 process null-faces
    # 6 add pseudo-edges
    # 7 NEW vertices
    # 8 start pos separator
    # 9 start neg separator
    # 10 endp pos separator
    # 11 endp neg separator
    # 12 add non-separator edges
    # 13 remove delete-tree
    # 14 reset status
    
    #vd.addLineSite( id3, id4 )
    t_after = time.time()
    calctime = t_after-t_before
    times.append(calctime)
    #if Nmax==0:
    #    Nmax=1
    #print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
    
    vod.setVDText2(times)
    
    vod.setAll()
    myscreen.render()
            

        
    print "PYTHON All DONE."

    myscreen.render()    
    myscreen.iren.Start()
