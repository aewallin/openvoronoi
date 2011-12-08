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

def sortedDict(adict):
    keys = adict.keys()
    keys.sort()
    l = []
    for key in keys:
        d=[]
        d.append(key)
        d.append(adict[key])
        l.append(d)
    return l

if __name__ == "__main__":  
    #print ocl.revision()
    #myscreen = ovdvtk.VTKScreen(width=1024, height=720) #(width=1920, height=1080)
    #ovdvtk.drawOCLtext(myscreen)
    
    #w2if = vtk.vtkWindowToImageFilter()
    #w2if.SetInput(myscreen.renWin)
    #lwr = vtk.vtkPNGWriter()
    #lwr.SetInput( w2if.GetOutput() )
    #w2if.Modified()
    #lwr.SetFileName("tux1.png")
    
    scale=1
    #myscreen.render()
    random.seed(42)
    far = 1
    #camPos = far
    zmult = 4
    # camPos/float(1000)
    #myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    #myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    #myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print vd.version()
    
    # for vtk visualization
    #vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    #vod.drawFarCircle()
    Nmax = 40000
    plist = randomGenerators(far, Nmax)    
    t_before = time.time() 
    n=0
    #id_list=[]
    for p in plist: 
        #print n," adding ",p
        vd.addVertexSite( p )
        n=n+1
    t_after = time.time()
    calctime = t_after-t_before
    #if Nmax==0:
    #    Nmax=1
    print " VD done in ", calctime," s, ", calctime/(Nmax*(math.log(Nmax)/math.log(2)))," s per n*log2(n())"
    stat = vd.getFaceStats()
    data=[]
    for s in stat:
        data.append( s[2] )
        #print s
    hist = histogram(data)
    hist2 = sortedDict(hist)
    print hist2

    #vod.setAll()
    #myscreen.render()
            

        
    print "PYTHON All DONE."

    #myscreen.render()    
    #myscreen.iren.Start()
