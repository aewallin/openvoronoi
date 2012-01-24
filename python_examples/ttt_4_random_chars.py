import ttt                # https://github.com/aewallin/truetype-tracer
import openvoronoi as ovd # https://github.com/aewallin/openvoronoi
import ovdvtk 

import time
import vtk
import random
import string

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
            if id_list[n] == 133630:
                vd.debug_on()
                vd.addLineSite( id_list[n], id_list[n_nxt],5)
                vod.setVDText2([1,1])
                vod.setAll()
                print "PYTHON All DONE."
                myscreen.render()   
                myscreen.iren.Start()
            else:
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
    wr.conic_line_subdivision = 50 # =10 increasesn nr of points to 366, = 5 gives 729 pts
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
    return out
    
def get_random_row(row_length):
    # construct some random strings
    chars = ""
    for n in range(row_length):
        c = random.choice(string.ascii_lowercase) # http://stackoverflow.com/questions/2823316/generate-a-random-letter-in-python
        chars+=c
    return chars

def get_scaled_translated_segs( chars, length, dx, n_row):
    # generate segs with scale 1
    ret = ttt_segments(  chars , 1)
    extents = ret[0]
    segs = ret[1]
    # translate so lower left corner is at (0,0)
    segs = translate(segs, -extents[0], -extents[2] )
    # scale to desired length
    current_length = extents[1]-extents[0]
    current_height = extents[3]-extents[2]
    segs = scale_segs(segs, current_length, length)
    # translate to final position
    start_y=-0.5
    dy = 1.5*float(length)/float(current_length)*current_height
    print " row to y= ",start_y+n_row*dy
    segs = translate(segs, dx, start_y+n_row*dy )
    # remove duplicate points
    segs = modify_segments(segs)
    return segs
    
if __name__ == "__main__":  
    #w=2500
    #h=1500
    
    #w=1920
    #h=1080
    w=1024
    h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    
    scale=1

    far = 1
    camPos = far
    zmult = 3
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    random.seed(42)
    row_length = 12
    n_rows = 3
    # works (seed=42) 10/3 11/3
    # segfaults: 12/3
    # length = 10 fits ca 4 rows
    
    length = 1
    dx = -0.5
    #dy = 0.2
    #starty = -0.5

    segs=[]
    for n in range(n_rows):
        chars = get_random_row(row_length)
        rowsegs = get_scaled_translated_segs( chars, length, dx, n)
        segs+=rowsegs
        print chars
    exit()
    #print chars
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    vod.textScale = 0.0002
    vod.vertexRadius = 0.0011
    vod.drawVertices=0
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vod.drawNullEdges = 1
    vd.setEdgeOffset(0.001)
    
    times = insert_many_polygons(vd,segs)
    vd.check()
    vod.setVDText2(times)
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
