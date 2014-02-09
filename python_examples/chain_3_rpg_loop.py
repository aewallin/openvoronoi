import openvoronoi as ovd
import ovdvtk  # comes with openvoronoi
import time
import vtk
import datetime
import math
import random
import os
import sys
import pickle
import gzip
import ovdgenerators as gens # comes with openvoronoi
import randompolygon as rpg # random polygon generator see https://github.com/aewallin/CGAL_RPG

def draw_vd(vd,times):
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
    #times=[]
    #times.append( 1 )
    #times.append( 1 )
    
    vod.setVDText2(times)
    
    vod.setAll()

    myscreen.render()   
    #w2if.Modified()
    #lwr.SetFileName("{0}.png".format(Nmax))
    #lwr.Write()
     
    myscreen.iren.Start()

def rpg_vd(Npts, seed, debug):
    far = 1
    
    vd = ovd.VoronoiDiagram(far,120)
    vd.reset_vertex_count()    
    poly = rpg.rpg(Npts, seed)

    pts=[]
    for p in poly:
        ocl_pt = ovd.Point( p[0], p[1] )
        pts.append(  ocl_pt )
        print ocl_pt
        
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts:
        #print " adding vertex ",m
        id_list.append( vd.addVertexSite( p ) )
        m=m+1
    """
    print "polygon is: "
    for idx in id_list:
        print idx," ",
    print "."
    """
    t_after = time.time()
    times.append( t_after-t_before )
    
    #print " pts inserted in ", times[0], " s"
    #print " vd-check: ",vd.check()
    if (debug):
        vd.debug_on()

    t_before = time.time()
    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0 # point 0 is the endpoint of the last segment
        #print " adding line-site ", id_list[n]," - ", id_list[n_nxt]
        vd.addLineSite( id_list[n], id_list[n_nxt])
    t_after = time.time()
    times.append( t_after-t_before )


    print " segs inserted in ", times[1], " s"
    is_valid = vd.check()
    print " vd-check: ", is_valid
    
    return [is_valid, vd, times]

def loop_run(Npts, max_seed, debug=False, debug_seed=-1):
    #Npts = 3
    #max_seed = 1000
    seed_range = range(max_seed)
    for seed in seed_range:
        debug2=debug
        if (seed==debug_seed):
            print "debug seed!"
            debug2 = True
        result = rpg_vd(Npts,seed,debug2)
        print "N=",Npts," s=",seed, " ok?=",result
        assert( result[0] == True )

def single_run(Npts, seed, debug=False):
    result = rpg_vd(Npts,seed,debug)
    print "N=",Npts," s=",seed, " ok?=",result
    assert( result[0] == True )
    return result
    

    
if __name__ == "__main__":  
    #loop_run(50,300)
    
    
    r = single_run(50,int(37))
    vd = r[1]
    pi = ovd.PolygonInterior(True)
    vd.filter_graph(pi)
    
    draw_vd(vd,r[2])
        
