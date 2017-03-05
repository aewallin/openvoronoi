import openvoronoi as ovd
import ovdvtk

from svgpathtools import svg2paths2


import time
import vtk
import math

import offset2vtk

def sub(p1, p2):
    return (p1[0] - p2[0], p1[0] - p2[0])
def add(p1, p2):
    return (p1[0] + p2[0], p1[0] + p2[0])
def mul(p, f):
    return (p[0] *f, p[1]*f)
def mag(p):
    return math.sqrt(p[0]*p[0]+p[1]*p[1])
def normalize(v):
    return v.mul(1.0/mag(v))
def tup(p):
    return (p.real, p.imag)
def dist(p1, p2):
    return mag(sub(p1, p2))

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
            pts.extend(self.adaptiveAdd(segment, 0, 1))
        if not path.isclosed():
            print("Path is not closed!:")
            p = path[0].point(0)
            self.addPoint(pts, p.real, p.imag)
        return pts

    def adaptiveAdd(self, segment, t0, t1):
        """Adaptively sample segment


        """
        p = tup(segment.point(t0))
        p2 = tup(segment.point(t1))
        v = sub(p2, p)
        d = mag(v)
        if d ==0:
            print("no dist")
            return [p]
        t1_2 = (t0+t1)/2.0
        half = add(p,mul(mul(v, 1.0/d), d*.5))
        half_p = tup(segment.point(t1_2))
        diff = dist(half, half_p)
        error = diff/d
        if d < .01  or error < .01:
            self.registerPoint(*p)
            return [p]
        else:
            #print("Error {}".format(int(100*error)))
            self.registerPoint(*half_p)
            first = self.adaptiveAdd(segment, t0, t1_2)
            res = []
            if (first):
                res = first
            second = self.adaptiveAdd(segment, t1_2, t1)
            if second:
                res.extend(second)
            return res


    def registerPoint(self, x,y):
        self.miny = min(y, self.miny)
        self.minx = min(x, self.minx)
        self.maxy = max(y, self.maxy)
        self.maxx = max(x, self.maxx)

    def addPoint(self, ls, x, y):
        registerPoint(x,y)
        ls.append((x,y))

    def parse(self):
        """Parse all paths to linear polygons

        Determine AABB, center and store transform
        """
        for path in self.paths:
            for continuous_subpath in path.continuous_subpaths():
                poly = self.path_to_poly(continuous_subpath)
                if poly:
                    self.polys.append(poly)

        center = (self.maxx + self.minx/2, (self.maxy + self.miny)/2)
        w = self.maxx - self.minx
        h = self.maxy - self.miny
        self.radius = math.sqrt(w*w + h*h)
        print("miny:{}".format(self.miny))
        print("minx:{}".format(self.minx))
        print("maxy:{}".format(self.maxy))
        print("maxx:{}".format(self.maxx))
        print("center:{}".format(center))
        print("radius:{}".format(self.radius))
        self.offsetx = center[0]
        self.offsety = center[1]



        adjusted = []
        for poly in self.polys:
            ls = []
            for p in poly:
                ls.append(self.offsetPoint(p))
   #             ls.reverse()
            adjusted.append(ls)
        self.polys = adjusted
    def offsetPoint(self, p):
        return (p[0] - self.offsetx,  p[1] - self.offsety)

    def bbox(self):
        return [
            self.offsetPoint( (self.minx-2, self.maxy +2)),
            self.offsetPoint( (self.maxx+2, self.maxy +2)),
            self.offsetPoint( (self.maxx+2, self.miny-2)),
            self.offsetPoint( (self.minx-2, self.miny-2)),
            ]


    def __getitem__(self, index):
        return self.polys[index]

def insert_polygon_points(vd, poly):
    pts=[]
    for p in poly:
        offsetx = center[0]
        offsety = center[1]
        for poly in self.polys:
            ls = []
            for p in poly:
                adjustedPoint = (p[0] - offsetx,  p[1] - offsety)
                ls.append(adjustedPoint)
   #             ls.reverse()
            adjusted.append(ls)
        self.polys = adjusted

    def __getitem__(self, index):
        return self.polys[index]

def insert_polygon_points(vd, poly):
    pts=[]
    for p in poly:
        #print("Inserting:{}".format(p))
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

def insert_bb(vd, svgr):
    polygon_ids =[]
    poly_id = insert_polygon_points(vd, svgr.bbox())
    insert_polygon_segments(vd,poly_id)

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

    far = svgr.radius * 1.5
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
    ofs_dist= step
    ofs = []
    for n in range(15):
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
