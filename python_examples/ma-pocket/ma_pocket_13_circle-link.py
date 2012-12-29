import openvoronoi as ovd
import ovdvtk
import ma_pocket_helper as maxh

import time
import vtk

if __name__ == "__main__":  
    #w=2500
    #h=1500
    
    #w=1920
    #h=1080
    w=1024
    h=1024
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
    vod.drawFarCircle()
    vod.setVDText2([1, 1])
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=1
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)
    
    linesegs = 1 # switch to turn on/off line-segments
    
    segs = []
    p1=ovd.Point(-0.1,-0.2)
    p2=ovd.Point(0.2,0.1)
    p3=ovd.Point(0.4,0.2)
    p4=ovd.Point(0.6,0.6)
    p5=ovd.Point(-0.6,0.3)
    pts = [p1,p2,p3,p4,p5]
    
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )
        m=m+1
    t_after = time.time()
    times.append( t_after-t_before )
    
    print "all Point sites inserted. "
    vd.check()
    
    t_before = time.time()    
    vd.addLineSite( id_list[0], id_list[1])
    vd.addLineSite( id_list[1], id_list[2])
    vd.addLineSite( id_list[2], id_list[3])
    vd.addLineSite( id_list[3], id_list[4])
    vd.addLineSite( id_list[4], id_list[0])
    vd.check()
    print "all Line sites inserted. "
    
    toolRadius = 0.01
    t_after = time.time()
    line_time = t_after-t_before
    if line_time < 1e-3:
        line_time = 1
    times.append( line_time )

    pi = ovd.PolygonInterior(True)
    vd.filter_graph(pi)
    
    # an interior offset, just for visualization
    of = ovd.Offset( vd.getGraph() ) 
    t_before = time.time()
    ofs = of.offset(toolRadius) 
    t_after = time.time()
    print "( OFFSET in %.3f ms.  )" % (1e3*(t_after-t_before))
    maxh.drawOffsets2(myscreen, ofs)
    myscreen.render()
    
    
    ma = ovd.MedialAxis()
    vd.filter_graph(ma)
    
    mapocket = ovd.MedialAxisPocket(vd.getGraph())
    mapocket.setCutWidth(0.01)
    mapocket.setCutterRadius(toolRadius)
    mapocket.debug(True)
    mapocket.run()
    mic_components = mapocket.get_mic_components()
    print "python processing of MICs:"
    for mic_list in mic_components:
        for n in range( len(mic_list) ):
            mic = mic_list[n]
            mic_center = mic[0]
            mic_radius = mic[1]
            assert( mic_radius > 0)
            if n == 0:
                #print "First MIC = ", mic_center," r = ",mic_radius
                ovdvtk.drawCircle( myscreen, mic_center, mic_radius , ovdvtk.red )
            else:
                #print "MIC = ", mic_center," r = ",mic_radius
                ovdvtk.drawCircle( myscreen, mic_center, mic_radius , ovdvtk.green )
    print "maxpocket done."
    
    vod.setAll()
    print "PYTHON All DONE."
    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
    myscreen.iren.Start()
