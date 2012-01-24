import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os


# this might be useful:
# http://code.activestate.com/recipes/409413-a-python-based-descriptive-statistical-analysis-to/

def drawLine(myscreen, p1, p2):
    myscreen.addActor( ovdvtk.Line( p1 = (p1.x,p1.y,0), p2=(p2.x,p2.y,0), color = ovdvtk.yellow ) )

def writeLargeFrame( myscreen, w2if, lwr, n , zoom=1):
    renderlarge = vtk.vtkRenderLargeImage()
    renderlarge.SetInput( myscreen.ren )
    renderlarge.SetMagnification(zoom)
    writer = vtk.vtkPNGWriter()
    writer.SetFileName("large_frame.png")
    writer.SetInputConnection( renderlarge.GetOutputPort() )
    writer.Write()
    print "Wrote large frame!"

def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd500_zoomout"+ ('%05d' % n)+".png"
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
    
def histogram(L):
    d = {}
    for x in L:
        if x in d:
            d[x] += 1
        else:
            d[x] = 1
    return d

if __name__ == "__main__":  

    width=1920
    height=1080
    pixmult=1
    width=pixmult*1024
    height=pixmult*1024
    myscreen = ovdvtk.VTKScreen(width, height)
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
    zmult = 4
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
    
    Nmax = 100
    plist = randomGenerators(far, Nmax)    
    times=[]
    t_before = time.time() 
    n=0
    id_list=[]
    for p in plist: 
        #print n," adding ",p
        id_list.append( vd.addVertexSite( p ) )
        n=n+1
    t_after = time.time()
    calctime = t_after-t_before
    times.append(calctime)
    if Nmax==0:
        Nmax=1
    #print " VD done in ", calctime," s, ", 1e6*calctime/float(Nmax)*(math.log(Nmax)/math.log(2.0))," us per n*log2(n)"
    stat = vd.getFaceStats()
    data=[]
    for s in stat:
        data.append( s[2] )
        #print s
    hist = histogram(data)
    print hist
    times.append(0)
    vod.setVDText2(times)
    vod.setAll()
    myscreen.render()
            
    #writeFrame(  w2if, lwr, 2 )
    #writeLargeFrame( myscreen, w2if, lwr, 2 , zoom=20)
    print "PYTHON All DONE."

    myscreen.render()    
    myscreen.iren.Start()
