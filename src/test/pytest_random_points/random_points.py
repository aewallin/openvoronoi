import openvoronoi as ovd
import time
import math
import random
import sys

def randomGenerators(far, Nmax):
    pradius = (1.0/math.sqrt(2))*far
    plist=[]
    for n in range(Nmax):
        x=-pradius+2*pradius*random.random()
        y=-pradius+2*pradius*random.random()
        plist.append( ovd.Point(x,y) )
    return plist


# this test computes the VD for N random points, with seed s
#
# call with random_points.py s N

if __name__ == "__main__": 
    s = int(sys.argv[1])
    n_pts = int(sys.argv[2])
    random.seed(s)
    Nmax = n_pts
    
    print "random_points.py N=",Nmax," seed=",s
    vd = ovd.VoronoiDiagram(1,120)
    
    plist = randomGenerators(1, Nmax)    
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
    print " VD done in ", calctime," s, ", 1e6*calctime/(float(Nmax)*(math.log(Nmax)/math.log(2.0)))," us per n*log2(n)"
    c = vd.check()
    print " VD check: ", c
    if c:
        exit(0)
    else:
        exit(-1)
