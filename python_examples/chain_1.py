import openvoronoi as ovd
import ovdvtk  # helper library for visualization using vtk

import time
import vtk
import datetime
import math
import random
import os
import sys
import pickle
import gzip

if __name__ == "__main__":
    # size of viewport in pixels
    # w=2500
    # h=1500

    # w=1920
    # h=1080
    w = 1024
    h = 800

    myscreen = ovdvtk.VTKScreen(width=w, height=h)
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version())

    w2if = vtk.vtkWindowToImageFilter()  # for screenshots
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInputConnection(w2if.GetOutputPort())
    # w2if.Modified()
    # lwr.SetFileName("tux1.png")

    scale = 1
    myscreen.render()
    random.seed(42)
    far = 1
    camPos = far
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos / float(1000), zmult * camPos)
    myscreen.camera.SetClippingRange(-(zmult + 1) * camPos, (zmult + 1) * camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)

    vd = ovd.VoronoiDiagram(far, 120)
    print ovd.version(), ovd.build_type()

    # for vtk visualization
    vod = ovdvtk.VD(myscreen, vd, float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()

    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices = 0
    vod.drawVertexIndex = 1
    vod.drawGenerators = 0
    vod.offsetEdges = 1
    vd.setEdgeOffset(0.05)

    linesegs = 1  # switch to turn on/off line-segments

    segs = []
    # ovd.Point(1,1)
    eps = 0.9
    p1 = ovd.Point(-0.1, -0.2)
    p2 = ovd.Point(0.2, 0.1)
    p3 = ovd.Point(0.4, 0.2)
    p4 = ovd.Point(0.6, 0.6)
    p5 = ovd.Point(-0.6, 0.3)

    pts = [p1, p2, p3, p4, p5]

    # t_after = time.time()
    # print ".done in {0:.3f} s.".format( t_after-t_before )
    times = []
    id_list = []
    m = 0
    t_before = time.time()
    print "inserting %d VertexSites one by one: " % len(pts)
    for p in pts:
        id_list.append(vd.addVertexSite(p))
        print "  %02d added vertex %3d at ( %1.3f, %1.3f )" % (m, id_list[m], p.x, p.y)
        m = m + 1

    t_after = time.time()
    times.append(t_after - t_before)
    print "all VertexSites inserted."

    vd.check()

    t_before = time.time()

    # vd.debug_on()
    print "inserting %d LineSites one by one: " % (len(id_list))

    for n in range(len(id_list)):
        print "  %02d source - target = %02d - %02d " % (n, id_list[n - 1], id_list[n])
        vd.addLineSite(id_list[n - 1], id_list[n])
    print "all LineSites inserted."
    vd.check()

    t_after = time.time()
    line_time = t_after - t_before
    if line_time < 1e-3:
        line_time = 1
    times.append(line_time)

    vod.setVDText2(times)

    err = vd.getStat()

    print "getStat() got errorstats for ", len(err), " points"
    if len(err) > 1:
        minerr = min(err)
        maxerr = max(err)
        print "  min error= ", minerr
        print "  max error= ", maxerr

    print "  num vertices: ", vd.numVertices()
    print "  num SPLIT vertices: ", vd.numSplitVertices()

    calctime = t_after - t_before

    vod.setAll()

    print "PYTHON All DONE."

    myscreen.render()
    # w2if.Modified()
    # lwr.SetFileName("{0}.png".format(Nmax))
    # lwr.Write() # write screenshot to file

    myscreen.iren.Start()
