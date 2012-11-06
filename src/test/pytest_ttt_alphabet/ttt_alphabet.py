import truetypetracer as ttt
import openvoronoi as ovd
import time
import sys

def translate(segs,x,y):
    out = []
    for seg in segs:
        seg2 = []
        for p in seg:
            p2 = []
            p2.append(p[0] + x)
            p2.append(p[1] + y)
            seg2.append(p2)
            #seg2.append(seg[3] + y)
        out.append(seg2)
    return out


def insert_polygon_points(vd, polygon):
    pts=[]
    for p in polygon:
        pts.append( ovd.Point( p[0], p[1] ) )
    id_list = []
    print "inserting ",len(pts)," point-sites:"
    m=0
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        print " ",m," added vertex ", id_list[ len(id_list) -1 ]
        m=m+1    
    return id_list

def insert_polygon_segments(vd,id_list):
    j=0
    print "inserting ",len(id_list)," line-segments:"
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
        vd.addLineSite( id_list[n], id_list[n_nxt])
        j=j+1

def modify_segments(segs):
    segs_mod =[]
    for seg in segs:
        first = seg[0]
        last = seg[ len(seg)-1 ]
        assert( first[0]==last[0] and first[1]==last[1] )
        seg.pop()
        seg.reverse()
        segs_mod.append(seg)
        #drawSegment(myscreen, seg)
    return segs_mod
    
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
    
def ttt_segments(text,scale,conic_subdiv):
    wr = ttt.SEG_Writer()

    # wr.scale = 3
    wr.arc = False
    wr.conic = False
    wr.cubic = False
    wr.conic_biarc_subdivision = 10 # this has no effect?
    wr.conic_line_subdivision = conic_subdiv # this increases nr of points 
    wr.cubic_biarc_subdivision = 10 # no effect?
    wr.cubic_line_subdivision = 10 # no effect?
    wr.setFont(3)
    wr.scale = float(1)/float(scale)
    s3 = ttt.ttt(text,wr) 
    segs = wr.get_segments()
    return segs
    
    
if __name__ == "__main__":  

    conic_subdiv = 200
    if len(sys.argv) == 2:
        conic_subdiv  = int(sys.argv[1])
    
    
    scale = 25000
    segs = ttt_segments(  "ABCDEFGHIJKLM", scale, conic_subdiv)
    segs2 = ttt_segments( "NOPQRSTUVWXYZ", scale, conic_subdiv)
    segs3 = ttt_segments( "abcdefghijklm", scale, conic_subdiv)
    #segs3 = ttt_segments( "m", 6400)
    segs4 = ttt_segments( "nopqrstuvwxyz", scale, conic_subdiv) # NOPQRSTUVWXYZ", 64000)
    segs5 = ttt_segments( "0123456789+-*/", scale, conic_subdiv)
    dx =  float(50000)/float(scale)
    xt=-0.3
    segs = translate(segs, xt*dx, 0.05*dx)
    segs = modify_segments(segs)
    segs2 = translate(segs2, xt*dx, -0.05*dx)
    segs2 = modify_segments(segs2)
    segs3 = translate(segs3, xt*dx, -0.15*dx)
    segs3 = modify_segments(segs3)
    segs4 = translate(segs4, xt*dx, -0.22*dx)
    segs4 = modify_segments(segs4)
    segs5 = translate(segs5, xt*dx, -0.32*dx)
    segs5 = modify_segments(segs5)
    
    vd = ovd.VoronoiDiagram(1,120)
        
    all_segs=segs+segs2 +segs3 +segs4+segs5
    insert_many_polygons(vd,all_segs)

    c = vd.check()
    print " VD check: ", c
    if c:
        exit(0)
    else:
        exit(-1)
