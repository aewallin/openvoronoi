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
    plist = gens.randomGenerators(far, Nmax)
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
    print "waiting for ",Nmax," segments..",
    sys.stdout.flush()
    segs = gens.randomSegments(far,Nmax)
    print ".done."
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
    # 10 -> 22 linesites
    # 15 -> 128
    # 20 -> 724
    # 25 -> 4096  lineseg benchmark crashes with segfault at 1448 generators!
    # 31 -> 32768
    # 36 -> 185363
    random.seed(1)
    exp_list = [0.5*x for x in range(5,Nmax_exp)]
    Nmax_list=[]
    for e in exp_list:
        Nmax_list.append( int( math.floor( (math.pow(2,e) ) ) ) )

    print Nmax_list
    #exit()
    csvWriter = csv.writer(open('results_rand_opt.csv', 'wb'), delimiter=',' )
    for Nmax in Nmax_list:
        #times = timeVoronoiPoints(Nmax)
        times = timeVoronoiSegs(Nmax)
        print Nmax," gens took ", sum(times) ," seconds, ", 1e6*float( sum(times) )/(float(Nmax)*float(math.log10(Nmax)))," us/n*log(n)"
        row = []
        row.append(Nmax)
        for t in times:
            row.append(t)
        csvWriter.writerow( row )
        

        
    print "PYTHON All DONE."


