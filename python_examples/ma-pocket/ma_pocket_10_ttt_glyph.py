import openvoronoi as ovd
import ovdvtk
import truetypetracer as ttt

import time
import vtk
import math
import sys
import datetime
import ngc_writer

import ma_pocket_helper as maxh


def drawPoint( myscreen, c, pcolor , rad = 0.02):
    ca = ovdvtk.Sphere(center=(c.x,c.y,0) , radius=rad, color=pcolor)
    myscreen.addActor(ca)

# rotate by cos/sin. from emc2 gcodemodule.cc
def rotate(x, y,  c,  s):
    tx = x * c - y * s;
    y = x * s + y * c;
    x = tx;
    return [x,y]
    
def drawArc(myscreen, pt1, pt2, r, cen,cw,arcColor):
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
    dlength = min(0.01, arclength/10)
    steps = int( float(arclength) / float(dlength))
    if steps==0: steps=1
    rsteps = float(1)/float(steps)
    dc = math.cos(-dtheta*rsteps) # delta-cos  
    ds = math.sin(-dtheta*rsteps) # delta-sin
    
    previous = pt1
    tr = [start.x, start.y]
    for i in range(steps):
        tr = rotate(tr[0], tr[1], dc, ds) #; // rotate center-start vector by a small amount
        x = cen.x + tr[0] 
        y = cen.y + tr[1] 
        current = ovd.Point(x,y)
        myscreen.addActor( ovdvtk.Line(p1=(previous.x,previous.y,0),p2=(current.x,current.y,0),color=arcColor) )
        #ngc_writer.xy_line_to(current.x,current.y)
        previous = current 

def rapid_to_next(myscreen, prv_tang, nxt_tang, c1, r1, c2, r2, prv, nxt):
    # rapid from prev, to nxt
    # while staying inside c1(r1) and c2(r)

    rad_default = 0.003
    rad = min( rad_default, 0.9*r1 , 0.9*r2)
    
    prv_tang.normalize()
    nxt_tang.normalize()

    prv_normal = -1*prv_tang.xy_perp()
    nxt_normal = nxt_tang.xy_perp()
    
    cen1 = prv + rad*prv_normal # + rad1*prv_tang
    cen2 = nxt - rad* nxt_normal #rapid_tang # + rad1*prv_tang

    rapid_tang = cen2-cen1
    rapid_tang.normalize()

    trg1 = cen1 + rad* rapid_tang.xy_perp() #prv_tang
    src2 = cen2 + rad* rapid_tang.xy_perp() 

    drawArc(myscreen, prv, trg1, rad, cen1, True, ovdvtk.blue) # lead-out arc
    ngc_writer.xy_arc_to(trg1.x, trg1.y, rad, cen1.x, cen1.y, True)
    
    ovdvtk.drawLine(myscreen, trg1, src2, ovdvtk.cyan) # rapid
    ngc_writer.xy_rapid_to(src2.x,src2.y)
    
    drawArc(myscreen, src2, nxt, rad, cen2, True, ovdvtk.lblue) # lead-in arc
    ngc_writer.xy_arc_to(nxt.x, nxt.y, rad, cen2.x, cen2.y, True)
    

def rapid_to_new_branch(myscreen, prv_tang, nxt_tang, c1, r1, c2, r2, prv, nxt):
    # rapid from prev, to nxt
    # while staying inside c1(r1) and c2(r)
    rad_default = 0.003
    rad1 = min( rad_default, 0.9*r1 ) # wrong? we get the new-branch r1 here, while we would want the old-branch r1
    rad2 = min( rad_default, 0.9*r2 )    
    prv_tang.normalize()
    nxt_tang.normalize()

    prv_normal = -1*prv_tang.xy_perp()
    nxt_normal = nxt_tang.xy_perp()
    
    cen1 = prv + rad1*prv_normal # + rad1*prv_tang
    
    cen2 = nxt - rad2* nxt_normal #rapid_tang # + rad1*prv_tang

    rapid_tang = cen2-cen1
    rapid_tang.normalize()
    
    trg1 = cen1 + rad1* prv_tang 
    src2 = cen2 - rad2* nxt_tang 
    
    drawArc(myscreen, prv, trg1, rad1, cen1, True, ovdvtk.orange) # lead-out arc
    ngc_writer.pen_up()
    ovdvtk.drawLine(myscreen, trg1, src2, ovdvtk.magenta) # rapid
    ngc_writer.xy_rapid_to(src2.x,src2.y)
    ngc_writer.pen_down()
    drawArc(myscreen, src2, nxt, rad2, cen2, True, ovdvtk.mag2) # lead-in arc

def final_lead_out(myscreen, prv_tang, nxt_tang, c1, r1, c2, r2, prv, nxt):
    rad_default = 0.003
    rad1 = min( rad_default, 0.9*r1 ) # wrong? we get the new-branch r1 here, while we would want the old-branch r1
    prv_tang.normalize()
    prv_normal = -1*prv_tang.xy_perp()
    cen1 = prv + rad1*prv_normal # + rad1*prv_tang
    trg1 = cen1 + rad1* prv_tang 
    drawArc(myscreen, prv, trg1, rad1, cen1, True, ovdvtk.red) # lead-out arc


# clear out initial MIC of a component
def spiral_clear(myscreen, out_tangent, in_tangent, c1, r1, c2, r2, out1, in1):
    # MIC: center= c1, radius=r1
    print "( spiral clear )"
    #ngc_writer.pen_up()
    # end spiral at in1
    # archimedean spiral
    # r = a + b theta
    in1_dir = in1-c1
    in1_theta = math.atan2(in1_dir.y,in1_dir.x)
    drawPoint( myscreen, c1, ovdvtk.red, 0.002 ) # red dot at center of MIC
    #drawPoint( myscreen, in1, ovdvtk.blue, 0.006 )
    # width = 2*pi*b
    # => b = width/(2*pi)
    b=0.01/(2*math.pi)
    # r = a + b in1_theta = r_max
    # =>
    # a = r_max-b*in1_theta
    a = r1 - b*in1_theta
    
    # figure out the start-angle
    theta_min = in1_theta
    theta_max = in1_theta
    dtheta = 0.1
    min_r = 0.001
    while True:
        r = a+b*theta_min
        if r < min_r:
            break
        else:
            theta_min = theta_min - dtheta
    #print "start_theta = ", theta_min
    
    Npts = (theta_max - theta_min)/dtheta
    Npts = int(Npts)
    #print "spiral has ",Npts," points"
    p = ovd.Point(c1)
    ngc_writer.xy_rapid_to(p.x,p.y)
    ngc_writer.pen_down()
    
    theta_end = 0
    for n in range(Npts+1):
        theta = theta_min + n*dtheta
        r = a + b*theta
        theta = theta - 2* abs(in1_theta - math.pi/2 )
        trg = c1 + r*ovd.Point( -math.cos(theta), math.sin(theta) )
        ovdvtk.drawLine(myscreen, p,trg,ovdvtk.pink)
        ngc_writer.xy_line_to(trg.x,trg.y)
        p = trg
        theta_end = theta
    
    
    # add a complete circle after the spiral.
    print "( spiral-clear: final circle )"
    Npts = (2*math.pi)/dtheta
    Npts = int(Npts)
    for n in range(Npts+2):
        theta = theta_end + (n+1)*dtheta
        #theta = theta_min + n*dtheta
        r = r1 #a + b*theta
        #theta = theta - 2* abs(in1_theta - math.pi/2 )
        trg = c1 + r*ovd.Point( -math.cos(theta), math.sin(theta) )
        ovdvtk.drawLine(myscreen, p,trg,ovdvtk.pink)
        ngc_writer.xy_line_to(trg.x,trg.y)
        
        #if n != Npts+1:
        #    drawPoint(myscreen, trg, ovdvtk.orange)
        #else:
        #    drawPoint(myscreen, trg, ovdvtk.orange,0.004)
        p = trg
        #if n == Npts-2:
        #    break
    print "( spiral clear done. )"

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

# return a list of points corresponding to an arc
# don't return the initial points, we already have that!
def arc_pts2(  pt1, pt2, r, cen,cw): # (start, end, radius, center, cw )
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
        pts.append(current)
        previous = current 
    return pts

# faster drawing of offsets using vtkPolyData
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
                cw=p[3]  # cw/ccw flag
                cen=p[2] # center
                r=p[1]   # radius
                p=p[0]   # target point
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
    #print len(ofs_points)," loops to render:"
    idx = 0
    last_idx = 0
        
    for of in ofs_points:
        epts  = of 
        segs=[]
        first = 1
        #print " loop with ", len(epts)," points"
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


def insert_polygon_points2(vd, polygon):
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
    #print vd.numFaces()," faces after all points inserted"
    return id_list

def insert_polygon_segments2(vd,id_list):
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        vd.addLineSite( id_list[n], id_list[n_nxt])

# give ofsets ofs
# insert points and line-segments in the vd
def insert_offset_loop(vd,ofs):
    polygon_ids =[]
    # create segs from ofs
    segs = []
    previous = ovd.Point()
    for ofloop in ofs:
        loop = []
        first = True
        for of in ofloop:
            #print of
            if first:
                #loop.append( of[0] )
                previous = of[0]
                first = False
            else:
                cw=of[3]  # cw/ccw flag
                cen=of[2] # center
                r=of[1]   # radius
                p=of[0]   # target point
                if r==-1: # r=-1 means line-segment
                    loop.append(p) #points.extend( [previous,p] ) #drawLine(myscreen, previous, p, lineColor)
                else: # otherwise we have an arc
                    loop.extend( arc_pts2( previous, p, r,cen,cw) )
                    #points.extend( arc_pts( previous, p, r,cen,cw) )
                previous = p
                #loop.append(p)
    
        segs.append(loop)
        
    #print segs
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

def insert_polygon_points(vd, pts):
    id_list = []
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
    return id_list

def insert_polygon_segments(vd,id_list):
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        vd.addLineSite( id_list[n], id_list[n_nxt])

# a simple class with a write method
class WritableObject:
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)

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
        poly_id = insert_polygon_points2(vd,poly)
        polygon_ids.append(poly_id)
    t_after = time.time()
    pt_time = t_after-t_before
    
    t_before = time.time()
    for ids in polygon_ids:
        insert_polygon_segments2(vd,ids)
    
    t_after = time.time()
    seg_time = t_after-t_before
    
    return [pt_time, seg_time]

# does point c lie left or right of the a->b line?
def is_right(a,b,c):
    det = (a.x-c.x)*(b.y-c.y) - (b.x-c.x)*(a.y-c.y)
    return det>0
    
if __name__ == "__main__":  
    #w=2500
    #h=1500
    w=1920
    h=1080
    #w=1024
    #h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )

    scale=1
    myscreen.render()
    #random.seed(42)
    far = 1
    camPos = far
    zmult = 1.8
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0.22, 0)
    
    # redirect stdout to file
    # example with redirection of sys.stdout
    foo = WritableObject()                   # a writable object
    sys.stdout = foo                         # redirection

    print "( Medial-Axis pocketing. )"
    print "( %s )" % (datetime.datetime.now())
    print "( OpenVoronoi %s  )" % (ovd.version())
    print "( TOOL/MILL,10,0,50 ) " # for cutsim
    print "( COLOR,0,255,255 ) " # for cutsim
    print "( STOCK/BLOCK,700.0000,400.0000,10.0000,350.0000,160.0000,5.0000 ) " # for cutsim
    
    
    toolRadius = 0.008
    cutWidth = toolRadius
    
    [segs, extents, scale] = get_scaled_segs( "P", 0.3) # geometry from ttt
    dx = -0.3
    dy = 0
    segs = translate(segs, dx, dy )
    
    vd = ovd.VoronoiDiagram(1,120)
    times = insert_many_polygons(vd,segs) # construct vd
    vd.check()
    
    print "( VD constructed in   %.3f s.  )" % (sum(times))
    
    pi = ovd.PolygonInterior(True) # for interior max/offsetting
    vd.filter_graph(pi)
    
    # an interior offset, just for visualization
    of = ovd.Offset( vd.getGraph() ) 
    t_before = time.time()
    ofs = of.offset(toolRadius) 
    t_after = time.time()
    print "( OFFSET in %.3f ms.  )" % (1e3*(t_after-t_before))
    
    
    
    maxh.drawOffsets2(myscreen, ofs)
    myscreen.render()
    #myscreen.iren.Start()
    
    # visualization of the VD
    vod3 = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod3.textScale = 0.0002
    vod3.vertexRadius = 0.0031
    vod3.drawVertices=0
    vod3.drawVertexIndex=1
    vod3.drawGenerators=0
    vod3.offsetEdges = 0
    vod3.setVDText2(times)
    
    # ma filter
    ma = ovd.MedialAxis(1.0)
    vd.filter_graph(ma)
    
    vod3.setAll()
    myscreen.render()   
    #myscreen.iren.Start()
    
    mapocket = ovd.MedialAxisPocket(vd.getGraph())
    mapocket.setCutterRadius( toolRadius )
    mapocket.setCutWidth( cutWidth ) # 0.3*toolRadius
    mapocket.debug(True)
    t_before = time.time()
    mapocket.run()
    
    mic_components = mapocket.get_mic_components()
    mic_list = mic_components[0] # we assume only one component here
    
    t_after = time.time()
    print "( MA-pocket done in %.3f s. Got %d MICs )" % ((t_after-t_before),len(mic_list) )
    
    maxmic = mic_list[0] # first MIC
    previous_center = maxmic[0]
    previous_radius = maxmic[1]
    
    # draw the initial MIC. to be cleared with a spiral-path
    #ovdvtk.drawCircle( myscreen, maxmic[0], maxmic[1] , ovdvtk.red )
    
    myscreen.render()
    #myscreen.iren.Start()
    ngc_writer.scale = 10/0.03
    ngc_writer.preamble()
    ngc_writer.feed = 2000
    print "T1 M6" # change to tool 1
    
    # the rest of the MICs are then cleared
    nframe=0
    first = True
    previous_out1 = ovd.Point()
    out_tangent = ovd.Point()
    in_tangent = ovd.Point()

    for n in range(1,len(mic_list)): # first MIC is already done
        mic = mic_list[n] 
        
        if 0: #nframe == 40: # interrupt, for debug
            break
            
        cen2 = mic[0] # new center
        r2 = mic[1]   # new radius
        out1 = mic[2] # bitangent point
        in1 = mic[3]
        out2 = mic[4]
        in2 = mic[5]
        in_tangent = in2-in1
        previous_center = mic[6] # old center
        previous_radius = mic[7] # old radius
        new_branch = mic[8]           # true/false indicates if we are starting on new branch
        prev_branch_center = mic[9]   # old branch MIC center
        prev_branch_radius = mic[10]  # old branch MIC radius

        # Make a positioning move so we are ready for cutting
        # the current MIC
        if first:
            # spiral-clear the start-MIC. The spiral should end at in1
            spiral_clear(myscreen, out_tangent, in_tangent, previous_center, previous_radius, cen2, r2, previous_out1, in1)
            #print "No rapid-move on first-iteration."
            first = False
        else:
            if new_branch:
                # new branch re-position move: pen_up(), rapid(), pen_down()
                rapid_to_new_branch(myscreen, out_tangent, in_tangent, prev_branch_center, prev_branch_radius , cen2, r2, previous_out1, in1)
            else:
                # normal arc-rapid-arc to next MIC
                rapid_to_next(myscreen, out_tangent, in_tangent, previous_center, previous_radius, cen2, r2, previous_out1, in1)
        
        # actual cutting of current MIC
        ovdvtk.drawLine(myscreen, in1, in2, ovdvtk.green)         # in bi-tangent
        ngc_writer.xy_line_to(in2.x,in2.y) # this bitangent is wrong just after spiral_clear!!
        #digits = 4
        
        # arc-cut
        # in2: start point
        # out2: end point
        # r2: radius
        # cen2: center
        # True: CW arc
        drawArc(myscreen, in2, out2, r2, cen2, True, ovdvtk.green)  # arc-cut
        #if (r2 > math.sqrt( math.pow(in2.x-out2.x,2) + math.pow(in2.y-out2.y,2) ) ):
        if is_right(in2,out2,cen2):
            # more than a half-circle, handle with two separate arc-moves
            print "(warning circle ambiguous)"
            """
            startvector = in2-cen2
            endvector = out2-cen2
            theta1 = math.atan2(startvector.x,startvector.y)
            theta2 = math.atan2(endvector.x,endvector.y)
            CIRCLE_FUZZ = 1e-9
            
            # FIXME: this would allow CCW cutting also
            #if (cw == False ): 
            #    while ( (theta2 - theta1) > -CIRCLE_FUZZ): 
            #        theta2 -= 2*math.pi
            #else:
            while( (theta2 - theta1) < CIRCLE_FUZZ): 
                theta2 += 2*math.pi
            dtheta = theta2-theta1
            theta_half = dtheta/2.0
            # create a point rotated from start by theta_half
            midvector = rotate(startvector.x, startvector.y, math.cos(-theta_half), math.sin(-theta_half))
            midpt = ovd.Point()
            midpt.x = cen2.x + midvector[0]
            midpt.y = cen2.y + midvector[1]
            #output two circles
            print "( ambiguous arc to %.4f, %.4f )" % ( ngc_writer.scale*out2.x, ngc_writer.scale*out2.y)
            print "( arcing via midpt %.4f, %.4f )" % ( ngc_writer.scale*midpt.x, ngc_writer.scale*midpt.y)
            print "( two circle output )"
            assert( not is_right(in2,midpt,cen2) ) # first arc OK
            assert( not is_right(midpt,out2,cen2) ) # second arc OK
            """
            #ngc_writer.xy_arc_to( midpt.x, midpt.y, r2, cen2.x, cen2.y, True )
            ngc_writer.xy_arc_to( out2.x, out2.y, -r2, cen2.x, cen2.y, True )
        else:
            # a normal case where one G2/3 move is sufficient
            ngc_writer.xy_arc_to( out2.x, out2.y, r2, cen2.x, cen2.y, True )

        ovdvtk.drawLine(myscreen, out2, out1, ovdvtk.green)  # out bi-tangent
        ngc_writer.xy_line_to(out1.x,out1.y)
        
        previous_out1 = out1     # this point is used as the start-point for the rapid on the next iteration
        out_tangent = out1-out2  # out-tangent direction
        
        if n == len(mic_list)-1: # end of operation. do a final lead-out arc.
            final_lead_out(myscreen, out_tangent, in_tangent, previous_center, previous_radius, cen2, r2, previous_out1, in1)
            #print "Final lead-out arc"

        nframe = nframe+1
        #myscreen.render()
    # end MIC-loop

    ngc_writer.postamble()
    sys.stdout = sys.__stdout__              # remember to reset sys.stdout!
    
    # write foo out to a file
    f = open('output.ngc', 'w')
    for item in foo.content:
        if item != '\n':
            print>>f, item
    f.close()
    
    print "python done."
    myscreen.render()   
    myscreen.iren.Start()
