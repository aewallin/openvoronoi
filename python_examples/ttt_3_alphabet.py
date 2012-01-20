import ttt
import openvoronoi as ovd
import ovdvtk
import time
import vtk

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
    #print "inserting ",len(pts)," point-sites:"
    m=0
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        #print " ",m," added vertex ", id_list[ len(id_list) -1 ]
        m=m+1    
    return id_list

def insert_polygon_segments(vd,id_list):
    j=0
    #print "inserting ",len(id_list)," line-segments:"
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0
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
    s3 = ttt.ttt(text,wr) 
    segs = wr.get_segments()
    return segs
    
    
if __name__ == "__main__":  
    #w=2500
    #h=1500
    
    w=1600
    h=1024
    #w=1024
    #h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.revision() )
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #w2if.Modified()
    #lwr.SetFileName("tux1.png")
    
    scale=1
    far = 1
    camPos = far
    zmult = 3
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    

    segs = ttt_segments(  "ABCDEFGHIJKLM", 64000)
    segs2 = ttt_segments( "NOPQRSTUVWXYZ", 64000)
    segs3 = ttt_segments( "abcdefghijklm", 64000)
    #segs3 = ttt_segments( "m", 6400)
    segs4 = ttt_segments( "nopqrstuvwxyz", 64000) # NOPQRSTUVWXYZ", 64000)
    segs5 = ttt_segments( "0123456789+-*/", 64000)
    #segs = ttt_segments(  "A", 64000)
    #segs2 = ttt_segments( "B", 64000)
    #segs2=[]
    segs = translate(segs, -0.6, 0.05)
    segs = modify_segments(segs)
    
    segs2 = translate(segs2, -0.6, -0.05)
    segs2 = modify_segments(segs2)
    
    segs3 = translate(segs3, -0.6, -0.15)
    segs3 = modify_segments(segs3)
    
    segs4 = translate(segs4, -0.6, -0.22)
    segs4 = modify_segments(segs4)
    
    segs5 = translate(segs5, -0.6, -0.32)
    segs5 = modify_segments(segs5)
    
    vd = ovd.VoronoiDiagram(far,120)
    print vd.version()
    
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    vod.textScale = 0.000002
    vod.vertexRadius = 0.0011
    vod.drawVertices=0
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vod.drawNullEdges = 0
    vd.setEdgeOffset(0.00001)
    
    all_segs=segs+segs2 +segs3 +segs4+segs5
    #all_segs=segs
    #all_segs=segs3 #+segs4
    #all_segs = segs3
    times = insert_many_polygons(vd,all_segs)
    vd.check()
    
    ovd.PolygonInterior( vd.getGraph() )
    ovd.MedialAxis( vd.getGraph() )
    
    vod.setVDText2(times)
    vod.setAll()
    
    print "PYTHON All DONE."

    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
     
    myscreen.iren.Start()
