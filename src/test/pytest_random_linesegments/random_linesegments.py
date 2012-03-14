import openvoronoi as ovd

import time
import math
import sys
import pickle
import gzip
import os

def test(Nmax):
    vd = ovd.VoronoiDiagram(1,120)
    print "random_linesegments.py:  ",Nmax," random segments."
    print "in directory: ", os.getcwd()
    sys.stdout.flush()
    f1 = "../../src/test/data/randomsegments_{0}.pickle.gz".format(Nmax)
    f2 = "../data/randomsegments_{0}.pickle.gz".format(Nmax)
    filename = ""
    if os.path.exists( f1 ):
        filename = f1 # "../src/test/data/randomsegments_{0}.pickle.gz".format(Nmax) #  load pre-computed segments 
    elif os.path.exists( f2 ):
        filename = f2
    # (produced with lineseg_dataset_generator.py)
    print "trying to open data file ",filename
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
    print sys.version
    print ovd.version()
    if c:
        print "Test passed."
        #exit(0)
    else:
        print "Test failed."
        #exit(-1)
    print "end of test()"
    
if __name__ == "__main__":  
    Nmax = 128 
    if len(sys.argv) == 2:
        Nmax  = int(sys.argv[1])
    test(Nmax)
    print "end of script"

