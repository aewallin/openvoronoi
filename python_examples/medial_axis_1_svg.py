import openvoronoi as ovd
import ovdvtk
import ngc_writer

import time
import vtk
import math

from SvgReader import SvgReader

def insert_polygon_points(vd, polygon):
    pts=[]
    for p in polygon:
        pts.append( ovd.Point( p[0], p[1] ) )
    id_list = []
    print("inserting ",len(pts)," point-sites:")
    m=0
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        #print " ",m," added vertex ", id_list[ len(id_list) -1 ]
        m=m+1
    return id_list

def insert_polygon_segments(vd,id_list):
    j=0
    print("inserting ",len(id_list)," line-segments:")
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
        #print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]
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


class NgcWriter:
    clearance_height= 20
    feed_height = 10
    feed = 200
    plunge_feed = 100
    metric = True
    scale = 1

    def __init__(self, filename="ouput.ngc"):
        self.filename = filename
        self.out = open(filename, 'w')

    def write(self, cmd):
        self.out.write("{}\n".format(cmd))
    def pen_up(self):
        self.write("G0 Z{}".format(self.feed_height))
    def xy_rapid_to(self, x,y):
        self.write("G0 X{} Y{}".format(x,y))
    def pen_down(self, z):
        self.write("G0 Z{}".format(z))
    def line_to(self, x,y,z):
        self.write("G1 X{} Y{} Z{} F{}".format(x,y,z, self.feed))

def printMedial(vd):
    maw = ovd.MedialAxisWalk(  vd.getGraph() )
    toolpath = maw.walk()
    ngc =NgcWriter()
    for chain in toolpath:
        n = 0
        for move in chain:
            for point in move:
                if n==0: # don't draw anything on the first iteration
                    p = point[0]
                    zdepth = scale*point[1]
                    ngc.pen_up();
                    ngc.xy_rapid_to( scale*p.x, scale*p.y );
                    ngc.pen_down( z= -zdepth )
                else:
                    p = point[0]
                    z = point[1]
                    ngc.line_to( scale*p.x, scale*p.y, scale*(-z) )
                n=n+1
    return

if __name__ == "__main__":
    w=800
    h=800

    myscreen = ovdvtk.VTKScreen(width=w, height=h)
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    scale=1
    myscreen.render()

    import sys
    svg = "../samples/rectangle.svg" if len(sys.argv) < 2 else sys.argv[1]
    svgr = SvgReader(svg, error_threshold = .6)
    svgr.parse()
    svgr.centerPolys()
    print(svgr.stats())

    far = svgr.radius * 1.2
    camPos = far
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)

    vd = ovd.VoronoiDiagram(far,120)
    print(ovd.version())

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


    times = insert_many_polygons(vd, svgr.polys)
    print("all sites inserted. ")
    print("VD check: ", vd.check())

    pi = ovd.PolygonInterior(  True )
    vd.filter_graph(pi)
    ma = ovd.MedialAxis()
    vd.filter_graph(ma)
    printMedial(vd)

    vod.setVDText2(times)
    vod.setAll()
    print("PYTHON All DONE.")
    myscreen.render()
    myscreen.iren.Start()
