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
    filename = current_dir + "/frames/vd_lineseg"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    lwr.Write()

def randomGenerators(far, Nmax):
    pradius = (1.0/math.sqrt(2))*far
    plist=[]
    for n in range(Nmax):
        x=-pradius+2*pradius*random.random()
        y=-pradius+2*pradius*random.random()
        plist.append( ovd.Point(x,y) )
    return plist
    
if __name__ == "__main__":  
    myscreen = ovdvtk.VTKScreen(width=1024, height=720) #(width=1920, height=1080)
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    
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
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    
    #vod.clearance_disk=1
    vod.vertexRadius = 0.005
    vod.textScale = 0.02
    Nmax = 60
    
    plist = randomGenerators(far, Nmax)
    
    t_before = time.time() 
    n=0
    id_list=[]
    npts = Nmax-1
    for p in plist: 
        print n," adding ",p
        id_list.append( vd.addVertexSite( p ) )
        n=n+1
    
    nstep = 10
    vd.addLineSite( id_list[0], id_list[30], nstep ) 
    
    t_after = time.time()
    calctime = t_after-t_before
    if Nmax==0:
        Nmax=1
    print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
    
    vod.setAll()
    myscreen.render()
    writeFrame( w2if, lwr, nstep )        

        
    print "PYTHON All DONE."

    myscreen.render()    
    myscreen.iren.Start()
