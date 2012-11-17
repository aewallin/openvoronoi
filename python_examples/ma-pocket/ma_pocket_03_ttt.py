import openvoronoi as ovd
import ovdvtk

import time
import vtk
import datetime
import math
import random
import os
import sys
import pickle
import gzip
import truetypetracer as ttt

def drawCircle(myscreen, c, r, circlecolor):
    ca = ovdvtk.Circle(center=(c.x,c.y,0) , radius=r, color=circlecolor, resolution=50 )
    myscreen.addActor(ca)


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
        print " ",m," added vertex ", id_list[ len(id_list) -1 ]
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
            print " ",j,"inserting segement ",id_list[n]," - ",id_list[n_nxt]

            if 0: # id_list[n] == 22871: #102187: # 102187/7 #115869: # 51456: 115869
                vd.debug_on()
                vd.addLineSite( id_list[n], id_list[n_nxt],7)
                vod.setVDText2([1,1])
                vod.setAll()
                vod.drawErrorVertices()
                #verts=[92555, 51680,92624,52559,51474,92620,52805]
                #for v in verts:
                    #print "drawing ",v
                    #print vod
                    #print dir(vod)
                #    vod.drawVertexIdx(v)
                print "PYTHON All DONE."
                myscreen.render()   
                myscreen.iren.Start()
            else:
                #pass
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
    wr.conic_line_subdivision = 200 # =10 increasesn nr of points to 366, = 5 gives 729 pts
    wr.cubic_biarc_subdivision = 10 # no effect?
    wr.cubic_line_subdivision = 10 # no effect?
    s3 = ttt.ttt(text,wr) 
    segs = wr.get_segments()
    ext = wr.extents
    return [ext, segs]

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
    

if __name__ == "__main__":  
    #w=2500
    #h=1500
    
    w=1920
    h=1080
    #w=1024
    #h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #w2if.Modified()
    #lwr.SetFileName("tux1.png")
    
    scale=1
    myscreen.render()
    random.seed(42)
    far = 1
    camPos = far
    zmult = 1.8
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0.22, 0)
    
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
    
    [segs, extents, scale] = get_scaled_segs( "M", 0.3)
    dx = -0.3
    dy = 0
    segs = translate(segs, dx, dy )
    
    times = insert_many_polygons(vd,segs)
    vd.check()
    vod.setVDText2(times)
    vod.setAll()
    
    pi = ovd.PolygonInterior(True)
    vd.filter_graph(pi)
    ma = ovd.MedialAxis()
    vd.filter_graph(ma)
    
    vod.setVDText2(times)
    
    vod.setAll()
    
    #myscreen.render()        
    #myscreen.iren.Start()

    
    mapocket = ovd.MedialAxisPocket(vd.getGraph())
    mapocket.debug(True)
    mapocket.setWidth(0.005)
    
    mapocket.run()
    
    #drawCircle( myscreen, maxmic[0], maxmic[1] , ovdvtk.red )
    mic_components = mapocket.get_mic_components()
    for mic_list in mic_components:
        
        nframe=0
        
        for n in range( len(mic_list) ):
            mic = mic_list[n]
            if n == 0:
                print "maxmic at ", mic[0]," r = ",mic[1]
                drawCircle( myscreen, mic[0], mic[1] , ovdvtk.red )
            else:
                drawCircle( myscreen, mic[0], mic[1] , ovdvtk.green )
            w2if.Modified()
            lwr.SetFileName("frames/%06d.png" % ( nframe ) )
            #lwr.Write()
            time.sleep(0.1)
            myscreen.render()
        print "mic done."
        
    """
    nframe=0
    while True:
        mic = mapocket.nxtMic()
        if nframe > 0:
            break
        if len(mic) == 2:
            drawCircle( myscreen, mic[0], mic[1] , ovdvtk.green )
            #w2if.Modified()
            #lwr.SetFileName("frames/%06d.png" % ( nframe ) )
            #lwr.Write()
        else:
            break # end operation when we don't get a valid MIC
        nframe = nframe+1
        
    print "mic done."
    """
    

        
    print "PYTHON All DONE."

    myscreen.render()   

     
    myscreen.iren.Start()
