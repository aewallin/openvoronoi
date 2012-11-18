import openvoronoi as ovd
import ovdvtk
import time
import vtk
import math

import offset2vtk

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
    vod.drawVertexIndex=1
    vod.drawGenerators=1
    vod.offsetEdges = 0 # for debug. a bool flag to set null-edge drawing on/off. use together with setEdgeOffset()
    vd.setEdgeOffset(0.05) # for debug. a non-zero value will draw null-edges as circular arcs
    # null-edges are an internal openvoronoi construction to avoid high-degree vertices in the VD-graph
    # they are not relevant for upstream or downstream algorithms
    
    # input points (vertices/sites)
    p1=ovd.Point(-0.1,-0.2)
    p2=ovd.Point(0.2,0.1)
    p3=ovd.Point(0.4,0.2)
    p4=ovd.Point(0.6,0.6)
    p5=ovd.Point(-0.6,0.3)
    pts = [p1,p2,p3,p4,p5] # a list of all points in the input
    
    #t_after = time.time()
    #print ".done in {0:.3f} s.".format( t_after-t_before )
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts: # add all points before adding line-segments
        id_list.append( vd.addVertexSite( p ) )
        #print m," added vertex", seg_id[0]
        m=m+1

    t_after = time.time()
    times.append( t_after-t_before )
    print "all point sites inserted. "
    print "VD check: ", vd.check()
    
    t_before = time.time()
    # now add line-segments, by using the integer indexes returned by vd.addVertexSite() above
    vd.addLineSite( id_list[0], id_list[1])
    vd.addLineSite( id_list[1], id_list[2])
    vd.addLineSite( id_list[2], id_list[3])
    vd.addLineSite( id_list[3], id_list[4])
    vd.addLineSite( id_list[4], id_list[0])
    vd.check()
    t_after = time.time()
    line_time = t_after-t_before
    if line_time < 1e-3:
        line_time = 1
    times.append( line_time )
    
    pi = ovd.PolygonInterior(  True )
    vd.filter_graph(pi)
    
    of = ovd.Offset( vd.getGraph() ) # pass the created graph to the Offset class
    of.str() # text output, for debug
    dists =[0.1,0.2, 0.21, 0.05, 0.14]
    ofs_loops=[]
    ofsl = []
    for d in dists:
        ofs_loops.extend( of.offset(d) )
        ofsl.extend( of.offset_loop_list(d) )
    print ofsl
    
    sorter = ovd.OffsetSorter(vd.getGraph())
    for loop in ofsl:
        sorter.add_loop( loop )
    
    sorter.sort_loops()
    
    offset2vtk.drawOffsets(myscreen, ofs_loops) # draw the generated offsets
    print "number of loops= ",len(ofs_loops)
    for loop in ofs_loops:
        first_vert=loop[0]
        print "loop at dist=", first_vert[2], " with ",len(loop)," vertices:"
        for v in loop[1:]:
            print " face ",v[4]

    vod.setVDText2(times)
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
