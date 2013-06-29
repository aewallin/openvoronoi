import openvoronoi as ovd
import math
import sys

import randompolygon as rpg # random polygon generator see https://github.com/aewallin/CGAL_RPG


def rpg_vd(Npts, seed, debug):
    far = 1
    vd = ovd.VoronoiDiagram(1,120)
    vd.reset_vertex_count()    
    poly = rpg.rpg(Npts, seed)

    pts=[]
    for p in poly:
        ocl_pt = ovd.Point( p[0], p[1] )
        pts.append(  ocl_pt )
        
    id_list = []
    for p in pts:
        id_list.append( vd.addVertexSite( p ) )


    for n in range(len(id_list)):
        n_nxt = n+1
        if n==(len(id_list)-1):
            n_nxt=0 # point 0 is the endpoint of the last segment
        vd.addLineSite( id_list[n], id_list[n_nxt])
    print " all segs inserted."
    is_valid = vd.check()
    print " vd-check: ", is_valid
    
    return is_valid

def loop_run(Npts, max_seed, debug=False, debug_seed=-1):
    seed_range = range(max_seed)
    for seed in seed_range:
        debug2=debug
        if (seed==debug_seed):
            print "debug seed!"
            debug2 = True
        result = rpg_vd(Npts,seed,debug2)
        print "N=",Npts," s=",seed, " ok?=",result
        assert( result == True )
        if ( not result ):
            exit(-1)

def single_run(Npts, seed, debug=False):
    result = rpg_vd(Npts,seed,debug)
    print "N=",Npts," s=",seed, " ok?=",result
    assert( result == True )
    return result

if __name__ == "__main__":  
    n_pts  = int(sys.argv[1])
    seed   = int(sys.argv[2])
    loop_run(n_pts,seed)
    exit(0)
    
    #r = single_run(50,int(37))
    #vd = r[1]
    #ovd.PolygonInterior(vd.getGraph(), True)
    #draw_vd(vd,r[2])
        
