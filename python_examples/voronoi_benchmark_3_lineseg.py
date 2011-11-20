import openvoronoi as ovd
import time
import math
import random
import csv
import sys

import ovdgenerators as gens

def timeVoronoiPoints(Nmax):
    far = 1
    vd = ovd.VoronoiDiagram(far, int( math.floor( math.sqrt(2)*math.sqrt(Nmax) ) ) )
    print "waiting for ",Nmax," random points..",
    sys.stdout.flush()
    t_before = time.time() 
    plist = gens.randomGenerators(far, Nmax)
    t_after = time.time()
    print ".done in ",(t_after - t_before)," seconds"
    #plist = regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    t_before = time.time() 
    for p in plist: 
        vd.addVertexSite( p )
    t_after = time.time()
    return [(t_after - t_before)]

def timeVoronoiSegs(Nmax):
    far = 1
    vd = ovd.VoronoiDiagram(far, int( math.floor( math.sqrt(2)*math.sqrt(Nmax) ) ) )
    print "waiting for ",Nmax," random segments..",
    sys.stdout.flush()
    t_before = time.time() 
    segs = gens.randomSegments(far,Nmax)
    t_after = time.time()
    print ".done in {0:.3f} seconds".format( (t_after - t_before) ) 
    id_list = []
    t_before = time.time() 
    for s in segs:
        segid=[]
        segid.append( vd.addVertexSite( s[0] ) )
        segid.append( vd.addVertexSite( s[1] ) )
        id_list.append( segid )
    t_after = time.time()
    point_site_time = (t_after - t_before)
    
    # insert linesegs
    t_before = time.time() 
    for s in id_list:
        vd.addLineSite( s[0], s[1] )
    t_after = time.time()
    line_site_time = (t_after - t_before)
    return [point_site_time , line_site_time]
    
if __name__ == "__main__":  
    far = 1
    Nmax_exp = 20
    # 10 -> 32 linesites
    # 14 -> 128
    # 18 -> 512
    # 20 -> 1024
    # 24 -> 4096  lineseg benchmark crashes with segfault at 1448 generators!
    # 28 -> 16384
    # 32 -> 65536
    # 33 -> 92681
    # 34 -> 131072
    # 35 -> 185363
    # 36 -> 262144
    
    random.seed(1)
    exp_list = [0.5*x for x in range(5,Nmax_exp+1)]
    Nmax_list=[]
    n=5
    for e in exp_list:
        Nmax_list.append( [ n, int( math.floor( (math.pow(2,e) ) ) ) ] )
        n=n+1
        
    #print Nmax_list
    #exit()
    csvWriter = csv.writer(open('results_rand_opt.csv', 'wb'), delimiter=',' )
    for case in Nmax_list:
        n=case[0]
        Nmax=case[1]
        #times = timeVoronoiPoints(Nmax)
        times = timeVoronoiSegs(Nmax)
        print n," voronoi-diagram for ",Nmax," sites took {0:.3f}".format(sum(times)) ," seconds, {0:.0f}".format( 1e6*float( sum(times) )/(float(Nmax)*float(math.log10(Nmax))) ) ,"us/n*log(n)"
        row = []
        row.append(Nmax)
        for t in times:
            row.append(t)
        csvWriter.writerow( row )
        

        
    print "PYTHON All DONE."


