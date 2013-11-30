import truetypetracer as ttt                # https://github.com/aewallin/truetype-tracer
import openvoronoi as ovd # https://github.com/aewallin/openvoronoi
#import ovdvtk 

import time
#import vtk
import random
import string
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
    print vd.numFaces()," faces after all points inserted"
    return id_list

def insert_polygon_segments(vd,id_list):
    j=0
    jmax=9999999 # for debugging, set jmax to the problematic case to stop algorithm in the middle
    print "inserting ",len(id_list)," line-segments:"
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        
        if (j<jmax):
            #vd.debug_on()
            print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
            
            if id_list[n] == 115869: # 51456: 115869
                vd.debug_on()
                vd.addLineSite( id_list[n], id_list[n_nxt],6)
                vod.setVDText2([1,1])
                vod.setAll()
                #verts=[92555, 51680,92624,52559,51474,92620,52805]
                #for v in verts:
                    #print "drawing ",v
                    #print vod
                    #print dir(vod)
                #    vod.drawVertexIdx(v)
                print "PYTHON All DONE."
                myscreen.render()   
                myscreen.iren.Start()
            else:
                #pass
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
    
def ttt_segments(text,scale):
    wr = ttt.SEG_Writer()

    # wr.scale = 3
    wr.arc = False
    wr.conic = False
    wr.cubic = False
    wr.scale = float(1)/float(scale)
    # "L" has 36 points by default
    wr.conic_biarc_subdivision = 10 # this has no effect?
    wr.conic_line_subdivision = 200 # =10 increasesn nr of points to 366, = 5 gives 729 pts
    wr.cubic_biarc_subdivision = 10 # no effect?
    wr.cubic_line_subdivision = 10 # no effect?
    s3 = ttt.ttt(text,wr) 
    segs = wr.get_segments()
    ext = wr.extents
    return [ext, segs]

def scale_segs(segs, current_length, desired_length):
    out=[]
    scale = float(desired_length) / float(current_length)
    for seg in segs:
        seg2 = []
        for p in seg:
            p2 = []
            p2.append(p[0] * scale)
            p2.append(p[1] * scale)
            seg2.append(p2)
            #seg2.append(seg[3] + y)
        out.append(seg2)
    return [out,scale]
    
def get_random_row(row_length):
    # construct some random strings
    chars = ""
    for n in range(row_length):
        c = random.choice(string.ascii_lowercase) # http://stackoverflow.com/questions/2823316/generate-a-random-letter-in-python
        chars+=c
        chars+=" "
    return chars

def get_scaled_segs( chars, length):
    # generate segs with scale 1
    ret = ttt_segments(  chars , 1)
    extents = ret[0]
    segs = ret[1]
    # translate so lower left corner is at (0,0)
    segs = translate(segs, -extents.minx, -extents.miny )
    # scale to desired length
    current_length = extents.maxx-extents.minx
    current_height = extents.maxy-extents.miny
    [segs,scale] = scale_segs(segs, current_length, length)

    # remove duplicate points
    segs = modify_segments(segs)
    return [segs, extents,scale]
    
if __name__ == "__main__": 
    print ttt.version()
    conic_subdiv = 200
    seed = 42
    if len(sys.argv) == 2:
        seed  = int(sys.argv[1]) # if seed was specified on command-line, take it
    
    random.seed( seed )
    row_length = 15
    n_rows = 10
    
    length = 1
    dx = -0.5
    start_y = -0.5
    current_y = start_y
    segs=[]
    for n in range(n_rows):
        chars = get_random_row(row_length)
        [rowsegs, extents, scale] = get_scaled_segs( chars, length)
        rowsegs_t = translate(rowsegs, dx, current_y )
        current_y = current_y + 1.1*(extents.maxy-extents.miny)*scale
        segs+=rowsegs_t
    
    vd = ovd.VoronoiDiagram(1,120)
    print ovd.version()
    
    times = insert_many_polygons(vd,segs)
    assert( vd.check() )
    print "PYTHON All DONE."
