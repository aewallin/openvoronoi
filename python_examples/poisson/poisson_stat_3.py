import openvoronoi as ovd
import time
import math
import random
import pickle
import gzip

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
    
def getVoronoiStats(Nmax,seed=1):
    random.seed(seed)
    far = 1
    vd = ovd.VoronoiDiagram(1,int( math.floor( math.sqrt(2)*math.sqrt(Nmax) ) ))
    
    plist = randomGenerators(1, Nmax)    
    t_before = time.time() 

    for p in plist: 
        vd.addVertexSite( p )

    t_after = time.time()
    calctime = t_after-t_before
    print " VD done in ", calctime," s, ", 1e6*calctime/(Nmax*(math.log(Nmax)/math.log(2)))," us per n*log2(n())"
    stat = vd.getFaceStats()
    data=[]
    for s in stat:
        data.append( s[2] )

    hist = histogram(data)
    hist2 = sortedDict(hist)
    return hist

def writeResults(seed,Nmax,data):
    filename = "poisson/N{0}_S{1}.pickle.gz".format(Nmax,seed)
    pstring = pickle.dumps( data,  2 ) # 2 is binary format
    f = gzip.open(filename, 'wb')
    f.write(pstring)
    f.close()

if __name__ == "__main__":  
    print ovd.version()
    Nmax=10000
    max_seed = 10
    for s in range(0,max_seed):
        d = getVoronoiStats(Nmax,s)
        print "seed= ",s," Nmax=",Nmax
        writeResults(s,Nmax,d)
    print "PYTHON All DONE."
