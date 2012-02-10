import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os

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

    scale=1
    random.seed(42)
    far = 1
    zmult = 4
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    Nmax = 40000
    print "calculating VD for ",Nmax," point-sites..."
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
    print " VD done in ", calctime," s, ", 1e6*calctime/(Nmax*(math.log(Nmax)/math.log(2)))," us per n*log2(n())"
    stat = vd.getFaceStats()
    data=[]
    for s in stat:
        data.append( s[2] )
    hist = histogram(data)
    hist2 = sortedDict(hist)
    print "distribution of faces with n edges in poisson vd with",Nmax," point-sites:"
    print hist2

            

        
    print "PYTHON All DONE."

    #myscreen.render()    
    #myscreen.iren.Start()
