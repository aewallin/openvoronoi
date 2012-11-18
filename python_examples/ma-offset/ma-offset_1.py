import openvoronoi as ovd
import ovdvtk
import truetypetracer as ttt   # https://github.com/aewallin/truetype-tracer

import time
import vtk
import math

import offset2vtk

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
    #print ocl.revision()
    #w=2500
    #h=1500
    
    #w=1920
    #h=1080
    #w=1024
    #h=1024
    w=800
    h=600
    
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )   
    
    scale=1
    myscreen.render()

    far = 1
    camPos = far
    zmult = 0.5
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()

    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)
    
    
    # segments from ttt
    segs = ttt_segments(  "L", 40000)
    #segs = translate(segs, -0.06, 0.05)
    segs = modify_segments(segs)
    
    times = insert_many_polygons(vd,segs)
    print "all sites inserted. "
    print "VD check: ", vd.check()
    
    pi = ovd.PolygonInterior(  True )
    vd.filter_graph(pi)
    
    of = ovd.Offset( vd.getGraph() ) # pass the created graph to the Offset class
    ofs_dist=0.0001
    ofs = []
    for n in range(50):
        ofsx = of.offset(ofs_dist) # generate offsets at the given distance.
        ofs.extend(ofsx)
        ofs_dist=ofs_dist+0.0001
        
    offset2vtk.drawOffsets2(myscreen, ofs) # draw the generated offsets
    
    
    #ma = ovd.MedialAxis(1)
    #vd.filter_graph(ma)
    
    vod.setVDText2(times)
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
