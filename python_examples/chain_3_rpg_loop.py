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
import ovdgenerators as gens
import rpg

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
        id_list.append( vd.addVertexSite( p ) )
        #print m," added vertex "
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
            
        vd.addLineSite( id_list[n], id_list[n_nxt])
    t_after = time.time()
    times.append( t_after-t_before )


    #print " segs inserted in ", times[1], " s"
    is_valid = vd.check()
    #print " vd-check: ", is_valid
    
    return is_valid

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
        assert( result == True )

def single_run(Npts, seed, debug=False):
    result = rpg_vd(Npts,seed,debug)
    print "N=",Npts," s=",seed, " ok?=",result
    assert( result == True )

if __name__ == "__main__":  
    loop_run(3,413,False,302)
    
    # s =337 problem!
    #single_run(3,336,1)
    #single_run(3,311)
