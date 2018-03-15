import openvoronoi as ovd
import ovdvtk

import time
import vtk
import math

import offset2vtk

import sys
import os
sys.path.append( os.path.join(os.path.dirname(__file__), "../"))
from SvgReader import SvgReader



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
    svgr.centerPolys()
    print(svgr.stats())

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
