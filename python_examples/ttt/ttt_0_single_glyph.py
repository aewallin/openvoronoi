import truetypetracer as ttt                # https://github.com/aewallin/truetype-tracer
import openvoronoi as ovd # https://github.com/aewallin/openvoronoi
import ovdvtk 

import math
import time
import vtk

def drawLoops(myscreen,loops,loopColor):
    # draw the loops
    nloop = 0
    for lop in loops:
        n = 0
        N = len(lop)
        first_point=[]
        previous=[]
        for p in lop:
            if n==0: # don't draw anything on the first iteration
                previous=p 
                first_point = p
            elif n== (N-1): # the last point
                myscreen.addActor( ovdvtk.Line(p1=(previous[0],previous[1],0),p2=(p[0],p[1],0),color=loopColor) ) # the normal line
                # and a line from p to the first point
                myscreen.addActor( ovdvtk.Line(p1=(p[0],p[1],0),p2=(first_point[0],first_point[1],0),color=loopColor) )
            else:
                myscreen.addActor( ovdvtk.Line(p1=(previous[0],previous[1],0),p2=(p[0],p[1],0),color=loopColor) )
                previous=p
            n=n+1
        print "rendered loop ",nloop, " with ", len(lop), " points"
        nloop = nloop+1

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
        print " ",m," added vertex ", id_list[ len(id_list) -1 ], " ({0}, {1})".format(p.x,p.y) 
        m=m+1   
    print vd.numFaces()," faces after all points inserted"
    return id_list

def insert_polygon_segments(vd,id_list):
    j=0
    #jmax=9999999 # for debugging, set jmax to the problematic case to stop algorithm in the middle
    print "inserting ",len(id_list)," line-segments:"
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        
        #if (j<jmax):
            #vd.debug_on()
        #print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
        if 0: #id_list[n]==304: #j == 0:
            print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
            vd.debug_on()
            vd.addLineSite( id_list[n], id_list[n_nxt])
            vod.setVDText2([1,1])
            #ovd.PolygonInterior( vd.getGraph() , True )
            #ovd.MedialAxis( vd.getGraph() )
            vod.drawIncidentVertexIds()
            for v in [1294, 1329]:
                vod.drawVertexIdx(v)
            for v in vd.getFaceVertices(3858):
                vod.drawVertexIdx(v)
                
            vod.setAll()
            print "PYTHON All DONE."
            myscreen.render()   
            myscreen.iren.Start()
        else:
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
    
def ttt_segments(text,scale, subdivision=100):
    wr = ttt.SEG_Writer()

    # wr.scale = 3
    #wr.arc = False
    #wr.conic = False
    #wr.cubic = False
    wr.scale = scale #float(1)/float(scale)
    # "L" has 36 points by default
    #wr.conic_biarc_subdivision = 10 # this has no effect?
    wr.conic_line_subdivision = subdivision  # =10 increasesn nr of points to 366, = 5 gives 729 pts
    
    # "C"  subd  points
    #       50   248
    #       40   311
    #       30   416
    #       20   627
    
    #wr.cubic_biarc_subdivision = 10 # no effect?
    #wr.cubic_line_subdivision = 10 # no effect?
    #wr.setFont(0)
    s3 = ttt.ttt(text,wr)
    exts = wr.extents
    #print exts
    segs = wr.get_segments()
    return (segs,exts)
    

if __name__ == "__main__":  
    print ovd.version() + " " + ovd.build_type()
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

    tscale = 2e-4
    (segs, exts) = ttt_segments(  "P", tscale) # 25000

    
    segs = center(segs, exts, tscale)
    segs = rotate(segs, 130)
    segs = modify_segments(segs)
    vd = ovd.VoronoiDiagram(far,120)
    
    
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    vod.textScale = 0.00002
    vod.vertexRadius = 0.0011
    vod.drawVertices=0
    vod.drawVertexIndex=1
    vod.drawGenerators=0
    vod.offsetEdges = 1
    vod.drawNullEdges = 1
    vd.setEdgeOffset(0.0001)
    
    all_segs=segs #+segs2 #+segs3 +segs4+segs5
    #all_segs=segs
    #all_segs=segs3 #+segs4
    #all_segs = segs6
    #insert_many_polygons(vd,all_segs)
    times = insert_many_polygons(vd,all_segs)
    #vd.check()
    vod.setVDText2(times)
    
    pi = ovd.PolygonInterior(  True )
    #vd.filter_graph(pi)
    ma = ovd.MedialAxis()
    #vd.filter_graph(ma)
    
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
