import openvoronoi as ovd    # https://github.com/aewallin/openvoronoi
import ttt                   # https://github.com/aewallin/truetype-tracer
import time

scale = 140
feed = 200
clearance_height=10
plunge_height = 2
plunge_feed = 100

def line_to(p,z):
    print "G1 X% 8.6f Y% 8.6f Z% 8.6f F%.0f" % (scale*p.x, scale*p.y, scale*z, feed)

# no arc-output in medial axis!
"""
def arc_to( p, r, cen, cw ):
    if (cw):
        print "G2 X% 8.5f Y% 8.5f R% 8.5f" % (scale*p.x, scale*p.y, scale*r)
    else:
        print "G3 X% 8.5f Y% 8.5f R% 8.5f" % (scale*p.x, scale*p.y, scale*r)
"""

def rapid_to(p):
    print "G0 X% 8.4f Y% 8.4f " % (scale*p.x, scale*p.y)
    


def pen_up():
    print "G0Z% 8.4f " % (clearance_height)



def pen_down():
    print "G0Z% 8.4f" % (plunge_height)
    print "G1Z0F% 8.0f" % (plunge_feed)
    
    
def plunge(p,z):
    print "G1 X% 8.4f Y% 8.4f Z% 8.4f" % (scale*p.x, scale*p.y, scale*z)
    

# naive medial-axis printer
# take all edges from the filtered vd
# don't print LINESITE or NULLEDGE edges
# print all other edges in the order they are output with (pen-up, rapid, pen-down) between every move
# NO depth information (i.e. we don't move in Z-axis)
def printMedial(vd):
    #edges = vd.getVoronoiEdges()

    maw = ovd.MedialAxisWalk(  vd.getGraph() )
    toolpath = maw.walk()
    #print w
    #print "( got %d edges )" % (len(edges))

    for chain in toolpath:
        # edges returned in this format:
        #            edge_data.append( point_list );   point_list = [ [p, clerance-disk] , ... ]
        #            edge_data.append( g[edge].type );
        #            edge_data.append( g[v1].status ); // source status
        #            edge_data.append( g[v2].status ); // target status

        n = 0
        for move in chain:
            #print move
            for point in move:
                #print point
                #exit()
                if n==0: # don't draw anything on the first iteration
                    p = point[0]
                    z = point[1]
                    pen_up();
                    rapid_to( p );
                    pen_down()
                    plunge( p, -z ) # now we are at the correct height, at the startpoint of the first move
                else:
                    p = point[0]
                    z = point[1]
                    line_to( p, -z )
                n=n+1

    return


# this function inserts point-sites into vd
def insert_polygon_points(vd, polygon):
    pts=[]
    for p in polygon:
        pts.append( ovd.Point( p[0], p[1] ) )
    id_list = []
    #print "inserting ",len(pts)," point-sites:"
    m=0
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        #print " ",m," added vertex ", id_list[ len(id_list) -1 ]
        m=m+1    
    return id_list

# this function inserts line-segments into vd
def insert_polygon_segments(vd,id_list):
    j=0
    #print "inserting ",len(id_list)," line-segments:"
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        #print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
        vd.addLineSite( id_list[n], id_list[n_nxt])
        j=j+1

# this function takes all segments from ttt and inserts them into vd
def insert_many_polygons(vd,segs):
    polygon_ids =[]
    t_before = time.time()
    for poly in segs:
        poly_id = insert_polygon_points(vd,poly)
        polygon_ids.append(poly_id)
    t_after = time.time()
    pt_time = t_after-t_before
    
    t_before = time.time()
    for ids in polygon_ids:
        insert_polygon_segments(vd,ids)
    
    t_after = time.time()
    seg_time = t_after-t_before
    
    return [pt_time, seg_time]

# this translates segments from ttt
def translate(segs,x,y):
    out = []
    for seg in segs:
        seg2 = []
        for p in seg:
            p2 = []
            p2.append(p[0] + x)
            p2.append(p[1] + y)
            seg2.append(p2)
        out.append(seg2)
    return out
    
# modify by deleting last point (since it is identical to the first point)
def modify_segments(segs):
    segs_mod =[]
    for seg in segs:
        first = seg[0]
        last = seg[ len(seg)-1 ]
        assert( first[0]==last[0] and first[1]==last[1] )
        seg.pop()
        seg.reverse() # to get interior or exterior offsets. Try commenting out this and see what happens.
        segs_mod.append(seg)
    return segs_mod
    
# get segments from ttt
def ttt_segments(text,scale):
    wr = ttt.SEG_Writer()
    wr.arc = False   # approximate arcs with lines
    wr.conic = False # approximate conic with arc/line
    wr.cubic = False # approximate cubic with arc/line
    
    wr.scale = float(1)/float(scale)
    ttt.ttt(text,wr) 
    segs = wr.get_segments()
    return segs

def preamble():
    print "G21 F% 8.0f" % (feed) # G20 F6 for inch
    print "G64 P0.001"
    pen_up()
    print "G0 X0 Y0"
    
def postamble():
    pen_up()
    print "M2"
    
if __name__ == "__main__":  
    vd = ovd.VoronoiDiagram(1,120) # parameters: (r,bins)  
    # float r  = radius within which all geometry is located. it is best to use 1 (unit-circle) for now.
    # int bins = number of bins for grid-search (affects performance, should not affect correctness)
    
    print "( TTT++",ttt.version()," )"
    print "( OpenVoronoi",vd.version()," )"
    
    # get segments from ttt. NOTE: must set scale so all geometry fits within unit-circle!
    segs = ttt_segments(  "EMC2", 15000) # (text, scale) all coordinates are divided by scale
    segs = translate(segs, -0.06, 0.05)
    segs = modify_segments(segs)
    
    times = insert_many_polygons(vd,segs) # insert segments into vd
    print "( all sites inserted in %.3f seconds )" % ( sum(times))
    print "( VD check: ", vd.check(), " )"
    
    ovd.PolygonInterior( vd.getGraph() ) # filter so that only polygon interior remains
    ovd.MedialAxis( vd.getGraph() ) # filter so that only medial axis remains
    
    preamble()
    
    printMedial( vd )
    
    postamble()
