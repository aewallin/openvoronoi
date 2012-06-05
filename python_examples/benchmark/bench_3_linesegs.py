import openvoronoi as ovd
import time
import math
import random
import csv
import sys
import gc

import ovdgenerators
    
def timeVoronoiPoints(Nmax):
    far = 1
    vd = ovd.VoronoiDiagram(far, int( math.floor( math.sqrt(2)*math.sqrt(Nmax) ) ) )
    print "waiting for {0} random points..".format(Nmax)
    sys.stdout.flush()
    t_before = time.time() 
    plist = ovdgenerators.randomGenerators(far, Nmax)
    t_after = time.time()
    print ".done in ",(t_after - t_before)," seconds"
    #plist = regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    t_before = time.time() 
    for p in plist: 
        vd.addVertexSite( p )
    t_after = time.time()
    return [(t_after - t_before)]

def timeVoronoiSegs(Nmax, segtype=1):
    far = 1
    vd = ovd.VoronoiDiagram(far, int( math.floor( math.sqrt(2)*math.sqrt(Nmax) ) ) )
    vd.reset_vertex_count()
    print "waiting for ",Nmax," random segments..",
    sys.stdout.flush()
    t_before = time.time() 
    segs=[]
    if segtype==1:
        segs = ovdgenerators.randomSegments(far,Nmax,seed=1)
    elif segtype==2:
        segs = ovdgenerators.randomSegments2(far,Nmax,seed=1)
    t_after = time.time()
    print ".done in {0:.3f} seconds".format( (t_after - t_before) ) 
    #print " first seg is ", segs[0]
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

def log2(x):
    return math.log10(x)/math.log10(2)
    
if __name__ == "__main__":  
    print ovd.version()
    far = 1
    Nmax_exp_start = 13
    Nmax_exp_end = 20
    # 10 -> 32 linesites
    # 14 -> 128
    # 18 -> 512
    # 20 -> 1024        (this takes 19s to generate the dataset!)
    # 24 -> 4096        ( this and below takes much longer!!)
    # 28 -> 16384
    # 32 -> 65536
    # 33 -> 92681
    # 34 -> 131072      
    # 35 -> 185363
    # 36 -> 262144
    
    random.seed(1)
    exp_list = [0.5*x for x in range(Nmax_exp_start,Nmax_exp_end+1)]
    Nmax_list=[]
    n=5
    for e in exp_list:
        Nmax_list.append( [ n, int( math.floor( (math.pow(2,e) ) ) ) ] )
        n=n+1
    print "Benchmarking for : "    
    print Nmax_list
    #exit()
    csvWriter = csv.writer(open('results_rand_opt.csv', 'wb'), delimiter=',' )
    for case in Nmax_list:
        n=case[0]
        Nmax=case[1]
        #times = timeVoronoiPoints(Nmax)
        times = timeVoronoiSegs(Nmax,segtype=1)
        print n," voronoi-diagram for ",Nmax," sites took {0:.3f}".format(sum(times)) ,
        print " seconds, {0:.2f}".format( 1e6*float( sum(times) )/(float(Nmax)*float(log2(Nmax))) ) ,"us/n*log2(n)"

        
        if len(times)==2:
            print "   ",2*Nmax," point-sites sites took {0:.3f}".format(times[0])," seconds, {0:.2f}".format( 1e6*float( times[0] )/(float(2*Nmax)*float(log2(2*Nmax))) ) ,"us/n*log2(n)"
            print "   ",Nmax," line-sites sites took {0:.3f}".format(times[1])," seconds, {0:.2f}".format( 1e6*float( times[1] )/(float(Nmax)*float(log2(Nmax))) ) ,"us/n*log2(n)"
        row = []
        row.append(Nmax)
        for t in times:
            row.append(t)
        csvWriter.writerow( row )
        #gc.clean()

        
    print "PYTHON All DONE."


