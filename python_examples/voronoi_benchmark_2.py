import openvoronoi as ovd
import time
import math
import random
import csv

def regularGridGenerators(far, Nmax):
    # REGULAR GRID
    rows = int(math.sqrt(Nmax))
    #print "rows= ",rows
    gpos=[-0.7*far ,  1.4*far/float(rows-1) ]  # start, stride
    plist = []
    for n in range(rows):
        for m in range(rows):
            x=gpos[0]+gpos[1]*n
            y=gpos[0]+gpos[1]*m
            # rotation
            """
            alfa = 0
            xt=x
            yt=y
            x = xt*math.cos(alfa)-yt*math.sin(alfa)
            y = xt*math.sin(alfa)+yt*math.cos(alfa)
            """
            plist.append( ovd.Point(x,y) )
    random.seed()
    random.shuffle(plist)
    return plist

def randomGenerators(far, Nmax):
    random.seed()
    pradius = (1.0/math.sqrt(2))*far
    plist=[]
    for n in range(Nmax):
        x=-pradius+2*pradius*random.random()
        y=-pradius+2*pradius*random.random()
        plist.append( ovd.Point(x,y) )
    return plist
    
def circleGenerators(far, Nmax, shuffle=1, seed=0):
    # POINTS ON A CIRCLE
    dalfa= float(2*math.pi)/float(Nmax-1)
    plist=[]
    radius=0.81234*float(far)
    for n in range(Nmax):
        x=float(radius)*math.cos(float(n)*float(dalfa))
        y=float(radius)*math.sin(float(n)*float(dalfa))
        plist.append( ovd.Point(x,y) )
    if shuffle:
        if seed:
            random.seed(seed)
        else:
            random.seed()
        random.shuffle(plist)
    return plist
    
def timeVoronoi(Nmax):
    far = 1
    vd = ovd.VoronoiDiagram(far, int( math.floor( math.sqrt(2)*math.sqrt(Nmax) ) ) )
    plist = randomGenerators(far, Nmax)
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
    #plist = randomGenerators(far, Nmax)
    plist = regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    for p in plist: 
        vd.pushVertexSite( p )
    t_before = time.time() 
    vd.run()
    t_after = time.time()
    return (t_after - t_before)

if __name__ == "__main__":  
    far = 1

    Nmax_exp = 40
    exp_list = [0.5*x for x in range(5,Nmax_exp)]
    Nmax_list=[]
    for e in exp_list:
        Nmax_list.append( int( math.floor( (math.pow(2,e) ) ) ) )

    print Nmax_list
    #exit()
    csvWriter = csv.writer(open('results_rand_opt.csv', 'wb'), delimiter=',' )
    for Nmax in Nmax_list:
        t = timeVoronoi_batch(Nmax)
        print Nmax," gens took ", t ," seconds, ", 1e6*float(t)/(float(Nmax)*float(math.log10(Nmax)))," us/n*log(n)"
        csvWriter.writerow( [ Nmax, t ] )
        

        
    print "PYTHON All DONE."


