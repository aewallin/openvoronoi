import truetypetracer as ttt                # https://github.com/aewallin/truetype-tracer
import openvoronoi as ovd # https://github.com/aewallin/openvoronoi
import sys
import time
import math

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

def rotate(segs,angle):
	theta = 2*math.pi*angle/float(360)
	out = []
	for seg in segs:
		seg2 = []
		for p in seg:
			p2 = []
			p2.append(p[0]*math.cos(theta)-p[1]*math.sin(theta) )
			p2.append(p[0]*math.sin(theta)+p[1]*math.cos(theta))
			seg2.append(p2)
		out.append(seg2)
	return out
    
def center(segs, exts, tscale):
	minx = tscale*exts.minx
	maxx = tscale*exts.maxx
	miny = tscale*exts.miny
	maxy = tscale*exts.maxy
	meanx = minx+0.5*(maxx-minx)
	meany = miny+0.5*(maxy-miny)
	out = translate(segs, -meanx, -meany)
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
        #print " ",m," added vertex ", id_list[ len(id_list) -1 ]
        m=m+1   
    #print vd.numFaces()," faces after all points inserted"
    return id_list

def insert_polygon_segments(vd,id_list):
    print "inserting ",len(id_list)," line-segments:"
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        vd.addLineSite( id_list[n], id_list[n_nxt])

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
    
def ttt_segments(text,scale):
    wr = ttt.SEG_Writer()

    # wr.scale = 3
    wr.arc = False
    wr.conic = False
    wr.cubic = False
    wr.scale = float(1)/float(scale)
    # "L" has 36 points by default
    wr.conic_biarc_subdivision = 10 # this has no effect?
    wr.conic_line_subdivision = 50 # =10 increasesn nr of points to 366, = 5 gives 729 pts
    wr.cubic_biarc_subdivision = 10 # no effect?
    wr.cubic_line_subdivision = 10 # no effect?
    s3 = ttt.ttt(text,wr) 
    segs = wr.get_segments()
    exts = wr.extents
    return (segs, exts)
    
    
if __name__ == "__main__":
    n=1
    rot=-20.333
    if len(sys.argv) == 2:
        n = int(sys.argv[1])
    elif len(sys.argv) == 3:
        n = int(sys.argv[1])
        rot = int(sys.argv[2])

    alphabet = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    char = alphabet[n]
    print "ttt_single_glyph.py: inserting glyph ", char

    scale = 3000
    #segs = ttt_segments(  char, scale)
    (segs, exts) = ttt_segments(  char , scale) # 25000

    segs = center(segs, exts, 1/float(scale) )
    segs = rotate(segs, rot)
    segs = modify_segments(segs)
    
    #segs = translate(segs, -0.5, -0.5)
    #segs = modify_segments(segs)
    vd = ovd.VoronoiDiagram(1,120)
    
    insert_many_polygons(vd,segs)
    
    # optional output to svg file
    filename = "svg_test.svg"
    #ovd.vd2svg(filename,vd)
    #print "wrote to file %s" % filename
    
    c = vd.check()
    print " VD check: ", c
    if c:
        exit(0)
    else:
        exit(-1)
