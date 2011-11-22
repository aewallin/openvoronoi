import openvoronoi as ovd
import math
import random

# this file has function for generating test-cases for OpenVoronoi

def randomSegments(far=1,Nmax=1,seed=1):
    random.seed(seed)
    segs = []
    for n in range(Nmax):
        seg = randomGenerators(far,2)
        while segmentIntersects(segs, seg):
            seg = randomGenerators(far,2)
        segs.append(seg)
    return segs

def randomSegments2(far=1,Nmax=1,seed=1):
    random.seed(seed)
    segs = []
    for n in range(Nmax):
        seg=randomSegment2(far,Nmax)
        while segmentIntersects(segs, seg):
            seg = randomSegment2(far,Nmax)
        segs.append(seg)
    return segs

def randomSegment2(far,Nmax):
    s1 = randomPoint(far)
    s2 = s1+randomPoint(float(far)/math.sqrt(float(Nmax)))
    while ( s2.norm() > far ):
        s1 = randomPoint(far)
        s2 = s1+randomPoint(float(far)/math.sqrt(float(Nmax)))

    return [s1,s2]
    
def intersects(s1,s2):
    """ test if line-segment s1 intersects with line-segment s2 """
    p1 = s1[0]
    p2 = s1[1]
    p = p1
    r = p2-p1
    q1 = s2[0]
    q2 = s2[1]
    q = q1
    s = q2-q1
    # t = (q-p) cross (s) / (r cross s)
    # u = (q-p) cross (r) / (r cross s)
    if ( r.cross(s) == 0 ): #parallel lines
        if ( (q-p).cross(r) == 0 ): #collinear
            return 1   
        else:
            return 0 # parallel lines that never intersect

    t = (q-p).cross(s) / (r.cross(s))
    u = (q-p).cross(r) / (r.cross(s))
    if ( (0<=t) and (t<=1) and (0<=u) and (u<=1) ):
        return 1
    return 0

# test if s intersects with any of the segments in seg
def segmentIntersects(segs, s):
    for sg in segs:
        if intersects(sg,s):
            return 1
    
    return 0 # no intersections found

def randomPoint(far):
    # random points
    #random.seed(seed)
    pradius = (1.0/math.sqrt(2))*far
    x=-pradius+2*pradius*random.random()
    y=-pradius+2*pradius*random.random()
    #plist.append( ovd.Point(x,y) )
    return ovd.Point(x,y)

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

def regularGridGenerators(far, Nmax, shuffle=0):
    # REGULAR GRID
    rows = int(math.sqrt(Nmax))
    #print "rows= ",rows
    gpos=[-0.7*far ,  1.4*far/float(rows-1) ]  # start, stride
    plist = []
    for n in range(rows):
        for m in range(rows):
            x=gpos[0]+gpos[1]*n
            y=gpos[0]+gpos[1]*m
            plist.append( ovd.Point(x,y) )
    if shuffle==1:
        random.seed()
        random.shuffle(plist)
    return plist
