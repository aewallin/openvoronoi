import openvoronoi as ovd
import time
import math
import random
import csv
import ovdgenerators

"""
def randomGenerators(far, Nmax):
    # random points
    #random.seed(seed)
    pradius = (1.0/math.sqrt(2))*far
    plist=[]
    for n in range(Nmax):
        x=-pradius+2*pradius*random.random()
        y=-pradius+2*pradius*random.random()
        plist.append( ovd.Point(x,y) )
    return plist
"""
def timeVoronoi(Nmax):
    far = 1
    vd = ovd.VoronoiDiagram(far, int( math.floor( math.sqrt(2)*math.sqrt(Nmax) ) ) )
    
    plist = ovdgenerators.randomGenerators(far, Nmax)
    #plist = regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    t_before = time.time() 
    for p in plist: 
        vd.addVertexSite( p )
    t_after = time.time()
    return (t_after - t_before)

def timeVoronoi_batch(Nmax):
    far = 1
    vd = ovd.VoronoiDiagram(far, int( math.floor( math.sqrt(16)*math.sqrt(Nmax) ) ) )
    
    #vd = ovd.VoronoiDiagram(far, 1200 ) # use a fixed number of bins
    plist = gens.randomGenerators(far, Nmax)
    #plist = gens.regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    for p in plist: 
        vd.pushVertexSite( p )
    t_before = time.time() 
    vd.run()
    t_after = time.time()
    return (t_after - t_before)

if __name__ == "__main__":  
    far = 1
    write_csv_output = 0 # turn on for output to CSV file
    Nmax_exp = 30
    exp_list = [0.5*x for x in range(5,Nmax_exp)]
    
    Nmax_list=[]
    for e in exp_list:
        Nmax_list.append( int( math.floor( (math.pow(2,e) ) ) ) )
    print ovd.version()
    print "Benchmarking for Nmax in: "
    print Nmax_list
    
    if write_csv_output:
        csvWriter = csv.writer(open('results_rand_189.csv', 'wb'), delimiter=',' )
    
    for Nmax in Nmax_list:
        t = timeVoronoi(Nmax)
        print Nmax," points took %.3f seconds = %.3f us/n*log2(n)" %  ( t,  1e6*float(t)/(float(Nmax)*float(math.log10(Nmax)/math.log10(2))) ) 
        if write_csv_output:
            csvWriter.writerow( [ Nmax, t ] )
        

        
    print "PYTHON All DONE."


