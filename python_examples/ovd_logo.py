import truetypetracer as ttt # https://github.com/aewallin/truetype-tracer
import openvoronoi as ovd    # https://github.com/aewallin/openvoronoi
import ovdvtk 

import time
import vtk
import math

# rotate by cos/sin. from emc2 gcodemodule.cc
def rotate(x, y,  c,  s):
    tx = x * c - y * s;
    y = x * s + y * c;
    x = tx;
    return [x,y]

# return a list of points corresponding to an arc
def arc_pts(  pt1, pt2, r, cen,cw): # (start, end, radius, center, cw )
    # draw arc as many line-segments
    start = pt1-cen
    end = pt2-cen
    theta1 = math.atan2(start.x,start.y)
    theta2 = math.atan2(end.x,end.y)
    alfa=[] # the list of angles
    da=0.1
    CIRCLE_FUZZ = 1e-9
    # idea from emc2 / cutsim g-code interp G2/G3
    if (cw == False ): 
        while ( (theta2 - theta1) > -CIRCLE_FUZZ): 
            theta2 -= 2*math.pi
    else:
        while( (theta2 - theta1) < CIRCLE_FUZZ): 
            theta2 += 2*math.pi
    
    dtheta = theta2-theta1
    arclength = r*dtheta
    dlength = min(0.001, arclength/10)
    steps = int( float(arclength) / float(dlength))
    rsteps = float(1)/float(steps)
    dc = math.cos(-dtheta*rsteps) # delta-cos  
    ds = math.sin(-dtheta*rsteps) # delta-sin
    
    previous = pt1
    tr = [start.x, start.y]
    pts=[]
    for i in range(steps):
        #f = (i+1) * rsteps #; // varies from 1/rsteps..1 (?)
        #theta = theta1 + i* dtheta
        tr = rotate(tr[0], tr[1], dc, ds) #; // rotate center-start vector by a small amount
        x = cen.x + tr[0] 
        y = cen.y + tr[1] 
        
        current = ovd.Point(x,y)
        #myscreen.addActor( ovdvtk.Line(p1=(previous.x,previous.y,0),p2=(current.x,current.y,0),color=arcColor) )
        pts.extend([previous, current])
        previous = current 
    return pts

def drawOffsets2(myscreen, ofs):
    # draw loops
    nloop = 0
    lineColor = ovdvtk.lgreen
    arcColor = ovdvtk.green #grass
    ofs_points=[]
    for lop in ofs:
        points=[]
        n = 0
        N = len(lop)
        first_point=[]
        previous=[]
        for p in lop:
            # p[0] is the Point
            # p[1] is -1 for lines, and r for arcs
            if n==0: # don't draw anything on the first iteration
                previous=p[0]
                #first_point = p[0]
            else:
                cw=p[3]
                cen=p[2]
                r=p[1]
                p=p[0]
                if r==-1: # r=-1 means line-segment
                    points.extend( [previous,p] ) #drawLine(myscreen, previous, p, lineColor)
                else: # otherwise we have an arc
                    points.extend( arc_pts( previous, p, r,cen,cw) )

                previous=p
            n=n+1
        ofs_points.append(points)
        #print "rendered loop ",nloop, " with ", len(lop), " points"
        nloop = nloop+1
        
    # now draw each loop with polydata
    oPoints = vtk.vtkPoints()
    lineCells=vtk.vtkCellArray()
    #self.colorLUT = vtk.vtkLookupTable()
    print len(ofs_points)," loops to render:"
    idx = 0
    last_idx = 0
        
    for of in ofs_points:
        epts  = of 
        segs=[]
        first = 1
        print " loop with ", len(epts)," points"
        for p in epts:
            oPoints.InsertNextPoint( p.x, p.y, 0)
            if first==0:
                seg = [last_idx,idx]
                segs.append(seg)
            first = 0
            last_idx = idx
            idx = idx + 1
            
        # create line and cells
        for seg in segs:
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, seg[0])
            line.GetPointIds().SetId(1, seg[1])
            #print " indexes: ", seg[0]," to ",seg[1]
            lineCells.InsertNextCell(line)
    
    linePolyData = vtk.vtkPolyData()
    linePolyData.SetPoints(oPoints)
    linePolyData.SetLines(lineCells)
    linePolyData.Modified() 
    linePolyData.Update()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(linePolyData)
    edge_actor = vtk.vtkActor()
    edge_actor.SetMapper(mapper)
    edge_actor.GetProperty().SetColor( ovdvtk.lgreen)
    myscreen.addActor( edge_actor )

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
        #print " ",m," added vertex ", id_list[ len(id_list) -1 ]
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
            #print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
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
    wr.conic_line_subdivision = 20 # =10 increasesn nr of points to 366, = 5 gives 729 pts
    wr.cubic_biarc_subdivision = 10 # no effect?
    wr.cubic_line_subdivision = 10 # no effect?
    s3 = ttt.ttt(text,wr) 
    segs = wr.get_segments()
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
    
    tscale = 20000
    segs = ttt_segments(  "Open", tscale)
    segs = translate(segs, -0.46, 0)
    segs2 = ttt_segments( "Voronoi", tscale)
    segs2 = translate(segs2, -0.25, -0.16)
    segs = modify_segments(segs)
    segs2 = modify_segments(segs2)

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
    
    all_segs=segs +segs2 #+segs3 +segs4+segs5
    #all_segs=segs
    #all_segs=segs3 #+segs4
    #all_segs = segs6
    #insert_many_polygons(vd,all_segs)
    times = insert_many_polygons(vd,all_segs)
    #vd.check()
    vod.setVDText2(times)
    
    pi_filt = ovd.PolygonInterior(  True )
    vd.filter_graph(pi_filt)
    ma_filt = ovd.MedialAxis(  )
    vd.filter_graph(ma_filt)
    vod.setAll() # this draws the medial-axis interior to the letters
    
    vd.filterReset()
    # another vod to draw outside text:
    vod2 = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    #vod.drawFarCircle()
    vod2.textScale = 0.0002
    vod2.vertexRadius = 0.0011
    vod2.drawVertices=0
    vod2.drawVertexIndex=0
    vod2.drawGenerators=0
    vod2.offsetEdges = 0
    vod2.drawNullEdges = 1
    #vd.setEdgeOffset(0.001)
    pi_filt = ovd.PolygonInterior( False )
    is_filt = ovd.IslandFilter(  )
    vd.filter_graph(pi_filt)
    vd.filter_graph(is_filt)
    vod2.setAll() 
    vd.filterReset()
    pi_filt = ovd.PolygonInterior(  False )
    vd.filter_graph(pi_filt)
    
    # some offsets
    of = ovd.Offset( vd.getGraph() ) # pass the created graph to the Offset class
    #of.str()
    ofs_list=[]
    t_before = time.time()
    for t in [0.02*x for x in range(1,4)]:
        ofs = of.offset(t)
        ofs_list.append(ofs)
    
    t_after = time.time()
    oftime = t_after-t_before
    for ofs in ofs_list:
        drawOffsets2(myscreen, ofs)
        
    
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
