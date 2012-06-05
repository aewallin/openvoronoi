import openvoronoi as ovd
import time
import math
import sys
import pickle
import gzip
import ovdgenerators as gens

if __name__ == "__main__":  
    print ovd.version()
    
    far = 1
    Nmax = int(math.pow(2,10)) # number of line-segments to generate
    
    print "waiting for ",Nmax," segments..",
    sys.stdout.flush()
    t_before = time.time()
    segs = gens.randomSegments(far,Nmax)
    t_after = time.time()
    print ".done in {0:.3f} s.".format( t_after-t_before )
    
    filename = "randomsegments_{0}.pickle.gz".format(Nmax)
    pstring = pickle.dumps( segs,  2 ) # 2 is binary format
    f = gzip.open(filename, 'wb')
    f.write(pstring)
    f.close()
    print "PYTHON All DONE."
