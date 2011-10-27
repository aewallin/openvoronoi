import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os


def drawFarCircle(myscreen, r, circleColor):
    myscreen.addActor( ovdvtk.Circle( center=(0,0,0), radius=r, color=circleColor ) )

def drawLine(myscreen, p1, p2):
    myscreen.addActor( ovdvtk.Line( p1 = (p1.x,p1.y,0), p2=(p2.x,p2.y,0), color = ovdvtk.yellow ) )

def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd500_zoomout"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    #lwr.Write()

def regularGridGenerators(far, Nmax):
    # REGULAR GRID
    rows = int(math.sqrt(Nmax))
    print "rows= ",rows
    gpos=[-0.7*far ,  1.4*far/float(rows-1) ]  # start, stride
    plist = []
    for n in range(rows):
        for m in range(rows):
            x=gpos[0]+gpos[1]*n
            y=gpos[0]+gpos[1]*m
            # rotation
            #alfa = 0
            #xt=x
            #yt=y
            #x = xt*math.cos(alfa)-yt*math.sin(alfa)
            #y = xt*math.sin(alfa)+yt*math.cos(alfa)
            plist.append( ovd.Point(x,y) )
    random.shuffle(plist)
    return plist

def randomGenerators(far, Nmax):
    pradius = (1.0/math.sqrt(2))*far
    plist=[]
    for n in range(Nmax):
        x=-pradius+2*pradius*random.random()
        y=-pradius+2*pradius*random.random()
        plist.append( ovd.Point(x,y) )
    return plist
    
def circleGenerators(far, Nmax):
    # POINTS ON A CIRCLE
    #"""
    #cpos=[50,50]
    #npts = 100
    dalfa= float(2*math.pi)/float(Nmax-1)
    #dgamma= 10*2*math.pi/npts
    #alfa=0
    #ofs=10
    plist=[]
    radius=0.81234*float(far)
    for n in range(Nmax):
        x=float(radius)*math.cos(float(n)*float(dalfa))
        y=float(radius)*math.sin(float(n)*float(dalfa))
        plist.append( ovd.Point(x,y) )
    #random.shuffle(plist)
    return plist
    
    
if __name__ == "__main__":  
    #print ocl.revision()
    myscreen = ovdvtk.VTKScreen()
    ovdvtk.drawOCLtext(myscreen)
    
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
    zmult = 5
    myscreen.camera.SetPosition(camPos/float(1000), camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    vd = ovd.VoronoiDiagram(far,120)
    print vd.version()
    vod = ovdvtk.VD(myscreen,vd,float(scale))
    drawFarCircle(myscreen, vd.getFarRadius(), ovdvtk.orange)
    
    Nmax = 4
    
    plist = randomGenerators(far, Nmax)
    #plist = regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    
    #plist = randomGenerators(far, Nmax) 
    #plist = []
    #plist.append( ovd.Point(0.0,0.1) )
    #plist.append( ovd.Point(0,0.9) )
    #plist.append( ovd.Point(-0.15, -0.15) )
    #+ regularGridGenerators(far, Nmax) + circleGenerators(far, Nmax)

    #plist = [ovd.Point(0,0)]
    
    t_before = time.time() 
    n=0
    id_list=[]
    for p in plist: 
        print n," adding ",p
        id_list.append( vd.addVertexSite( p ) )
        n=n+1
    id1 = id_list[0]
    id2 = id_list[1]
    id3 = id_list[2]
    id4 = id_list[3]
    print "add segment ",id1, " to ", id2
    vd.addLineSite( id1, id2 )
    #drawLine(myscreen, plist[0], plist[1] )
    
    #vd.addLineSite3( id3, id4 )
    #drawLine(myscreen, plist[2], plist[3] )
    t_after = time.time()
    calctime = t_after-t_before
    if Nmax==0:
        Nmax=1
    print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
    
    vod.setAll(vd)
    myscreen.render()
            

        
    print "PYTHON All DONE."

    myscreen.render()    
    myscreen.iren.Start()
