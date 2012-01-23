import openvoronoi as ovd

import time
import math
import sys
import pickle
import gzip

if __name__ == "__main__":  
    Nmax  = int(sys.argv[1])
    vd = ovd.VoronoiDiagram(1,120)
    print "random_linesegments.py:  ",Nmax," random segments."
    sys.stdout.flush()
    filename = "../src/test/data/randomsegments_{0}.pickle.gz".format(Nmax) #  load pre-computed segments 
    # (produced with lineseg_dataset_generator.py)
    f = gzip.open(filename, 'rb')
    pstring = f.read()
    segs = pickle.loads( pstring )
    f.close()

    id_list = []
    for seg in segs:
        seg_id=[]
        seg_id.append( vd.addVertexSite( seg[0] ) )
        seg_id.append( vd.addVertexSite( seg[1] ) )
        id_list.append( seg_id )
    print "all point-sites inserted."
    assert( vd.check() )
    nsegs = Nmax
    for s in id_list:
        vd.addLineSite(s[0],s[1])
    print "all line-sites inserted."
    c = vd.check()
    print " VD check: ", c
    if c:
        exit(0)
    else:
        exit(-1)
    

