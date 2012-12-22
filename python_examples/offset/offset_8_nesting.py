import openvoronoi as ovd
import ovdvtk
import time
import vtk
import math

import offset2vtk
import truetypetracer as ttt   # https://github.com/aewallin/truetype-tracer

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

def ttt_segments(text,scale):
    wr = ttt.SEG_Writer()
    wr.arc = False
    wr.conic = False
    wr.cubic = False

    wr.conic_biarc_subdivision = 10 # this has no effect?
    wr.conic_line_subdivision = 50 # this increases nr of points 
    wr.cubic_biarc_subdivision = 10 # no effect?
    wr.cubic_line_subdivision = 10 # no effect?
    wr.setFont(3)
    
    wr.scale = float(1)/float(scale)
    ttt.ttt(text,wr) 
    segs = wr.get_segments()
    return segs

def modify_segments(segs):
    segs_mod =[]
    for seg in segs:
        first = seg[0]
        last = seg[ len(seg)-1 ]
        assert( first[0]==last[0] and first[1]==last[1] )
        seg.pop()
        seg.reverse() # to get interior or exterior offsets
        segs_mod.append(seg)
        #drawSegment(myscreen, seg)
    return segs_mod

if __name__ == "__main__":  
    #w=2500 # screen resolution for big screens
    #h=1500
    
    #w=1920
    #h=1080
    w=1024
    h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) # a VTK window for drawing 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )   # the OpenVoronoi text, revision, and date
    
    scale=1
    myscreen.render()

    far = 1
    camPos = far
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    # these actions on the vod object control how the VD is drawn using VTK
    vod.drawFarCircle()
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0 # for debug. a bool flag to set null-edge drawing on/off. use together with setEdgeOffset()
    vd.setEdgeOffset(0.05) # for debug. a non-zero value will draw null-edges as circular arcs
    # null-edges are an internal openvoronoi construction to avoid high-degree vertices in the VD-graph
    # they are not relevant for upstream or downstream algorithms
    
    print "all sites inserted. "
    print "VD check: ", vd.check()
    
    # segments from ttt
    segs = ttt_segments(  "H", 30000)
    segs = translate(segs, -0.06, -0.05)
    segs = modify_segments(segs)
    
    # a box around the character
    box = 0.1
    p0 = [-box,-box]
    p1 = [box,-box]
    p2 = [box,box]
    p3 = [-box,box]
    seg0 = [p3,p2,p1,p0]
    segs.append(seg0)
    
    times = insert_many_polygons(vd,segs)
    
    pi = ovd.PolygonInterior(  False )
    vd.filter_graph(pi)
    
    
    of = ovd.Offset( vd.getGraph() ) # pass the created graph to the Offset class
    of.str() # text output, for debug
    dists =[ 0.004*x for x in range(1,15) ]
    ofs_loops=[]
    ofsl = []
    for d in dists:
        d_offsets = of.offset(d)
        if len(d_offsets) != 0:
            ofs_loops.extend( d_offsets ) # compute offset at d, and add to loop-list
            ofsl.extend( of.offset_loop_list(d) )
    #print ofsl
    
    sorter = ovd.OffsetSorter(vd.getGraph())
    for loop in ofsl:
        sorter.add_loop( loop )
    
    sorter.sort_loops()
    # this generates a graph test.dot
    # to generate a png image run
    #  dot -Tpng test.dot > test.png
    ofs_loops2 = sorter.get_loops()
    print "number of loops= ",len(ofs_loops2)
    
    offset2vtk.drawOffsets(myscreen, ofs_loops2) # draw the generated offsets
    
    """
    for loop in ofs_loops:
        first_vert=loop[0]
        print "loop at dist=", first_vert[2], " with ",len(loop)," vertices:"
        for v in loop[1:]:
            print " face ",v[4]
    """
    vod.setVDText2(times)
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
