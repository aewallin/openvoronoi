import openvoronoi as ovd
import ovdvtk

from svgpathtools import svg2paths2


import time
import vtk
import math

import offset2vtk

class SvgReader:

    def __init__(self, filename):
        self.paths, self.attributes, svg_attributes = svg2paths2(filename)
        self.polys = []
        self.miny = self.minx = float("+inf")
        self.maxy = self.maxx = float("-inf")
        #print("svg_attributes[idx]:{}".format(svg_attributes))

    def path_to_poly(self, path):
        """Return a list of points that forms a closed polygon

        TODO:
        If the path is not a line the steps parameter will divide the
        segment into \p steps steps
        """
        pts = []
        for segment in path:
            p = segment.point(0)
            self.addPoint(pts, p.real, p.imag)
        if not path.isclosed():
            print("Path is not closed!:")
            p = path[0].point(0)
            self.addPoint(pts, p.real, p.imag)
        return pts

    def addPoint(self, ls, x, y):
        self.miny = min(y, self.miny)
        self.minx = min(x, self.minx)
        self.maxy = max(y, self.maxy)
        self.maxx = max(x, self.maxx)
        ls.append((x,y))

    def parse(self):
        """Parse all paths to linear polygons

        Determine AABB, center and store transform
        """
        for path in self.paths:
            self.polys.append(self.path_to_poly(path))

        center = (self.maxx + self.minx/2, (self.maxy + self.miny)/2)
        self.radius = math.sqrt((self.minx-center[0])*(self.minx-center[0]) + (self.miny-center[1])*(self.miny-center[1]))
        print("miny:{}".format(self.miny))
        print("minx:{}".format(self.minx))
        print("maxy:{}".format(self.maxy))
        print("maxx:{}".format(self.maxx))
        print("center:{}".format(center))
        print("radius:{}".format(self.radius))


        adjusted = []
        for poly in self.polys:
            ls = []
            for p in poly:
                adjustedPoint = (p[0] - center[0],  p[1] - center[1])
                ls.append(adjustedPoint)
            adjusted.append(ls)
        self.polys = adjusted

    def __getitem__(self, index):
        return self.polys[index]

def insert_polygon_points(vd, poly):
    pts=[]
    for p in poly:
        print("Inserting:{}".format(p))
        pts.append( ovd.Point( *p ))
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

def insert_many_polygons(vd, svgr):
    polygon_ids =[]
    t_before = time.time()
    for poly in svgr:
        poly_id = insert_polygon_points(vd, poly)
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
    w=1024
    h=1024
    import sys
    svg = "../../samples/rectangle.svg" if len(sys.argv) < 2 else sys.argv[1]

    svgr = SvgReader(svg)
    svgr.parse()


    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )   
    
    scale=1
    myscreen.render()

    far = svgr.radius * 1.1
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

    times = insert_many_polygons(vd,svgr)
    print "all sites inserted. "
    print "VD check: ", vd.check()
    
    pi = ovd.PolygonInterior(  True )
    vd.filter_graph(pi)
    
    of = ovd.Offset( vd.getGraph() ) # pass the created graph to the Offset class
    step = svgr.radius *.01
    ofs_dist= step;
    ofs = []
    for n in range(50):
        ofsx = of.offset(ofs_dist) # generate offsets at the given distance.
        ofs.extend(ofsx)
        ofs_dist = ofs_dist + step
        
    offset2vtk.drawOffsets2(myscreen, ofs) # draw the generated offsets
    
    
    #ma = ovd.MedialAxis(1)
    #vd.filter_graph(ma)
    
    vod.setVDText2(times)
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
