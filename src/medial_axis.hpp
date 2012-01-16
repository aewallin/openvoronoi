/*  
 *  Copyright 2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  
 *  This file is part of OpenVoronoi.
 *
 *  OpenVoronoi is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  OpenVoronoi is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with OpenVoronoi.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <string>
#include <iostream>

#include <boost/python.hpp>

#include "graph.hpp"
#include "site.hpp"

namespace ovd
{

struct medial_filter {
    medial_filter(HEGraph& gi) : g(gi) { }
    bool operator()(const HEEdge& e) const {
        if (g[e].type == LINESITE || g[e].type == NULLEDGE) 
            return true;
        if (g[e].type == SEPARATOR)
            return false;
            
        if (both_endpoints_positive(e)) // these are interior edges which we keep.
            return true;
        
        // this leaves us with edges where one end connects to the polygon (dist==0)
        // and the other end does not.
        // figure out the angle between the adjacent line-segments and decide based on the angle.
        if (segments_parallel(e))
            return false;

        return true; // otherwise we keep the edge
    }
    bool both_endpoints_positive(HEEdge e) const {
        HEVertex src = g.source(e);
        HEVertex trg = g.target(e);
        return (g[src].dist()>0) && (g[trg].dist()>0);
    }
    bool segments_parallel( HEEdge e ) const {

        HEVertex endp1 = find_endpoint(e);
        HEVertex endp2 = find_endpoint( g[e].twin );
        // find the segments
        HEEdge e1 = find_segment(endp1);
        HEEdge e2 = find_segment(endp2);
        e2 = g[e2].twin; // this makes the edges oriented in the same direction 
        
        double dotprod = edge_dotprod(e1,e2);
        return fabs(dotprod)>0.8;
    }
    
    double edge_dotprod(HEEdge e1, HEEdge e2) const {
        HEVertex src1 = g.source(e1);
        HEVertex trg1 = g.target(e1);
        HEVertex src2 = g.source(e2);
        HEVertex trg2 = g.target(e2);
        Point sp1 = g[src1].position;
        Point tp1 = g[trg1].position;
        Point sp2 = g[src2].position;
        Point tp2 = g[trg2].position;
        
        Point dir1 = tp1-sp1;
        Point dir2 = tp2-sp2;
        dir1.normalize();
        dir2.normalize();
        return dir1.dot(dir2);
    }
    
    
    HEEdge find_segment(HEVertex v) const {
        BOOST_FOREACH(HEEdge e, g.out_edges(v)) {
            if ( g[e].type == LINESITE )
                return e;
        }
        assert(0);
        exit(-1);
        return HEEdge();
    }
    
    HEVertex find_endpoint(HEEdge e) const {
        HEEdge next = g[e].next;
        HEEdge prev = g.previous_edge(e);
        HEVertex endp;
        if ( g[next].type == NULLEDGE ) {
            endp = g.target(next);
            assert( g[endp].type == ENDPOINT );
            
        } else if ( g[prev].type == NULLEDGE ) {
            endp = g.source(prev);
            assert( g[endp].type == ENDPOINT );
        } else {
            assert(0);
            exit(-1);
        }
        return endp;
    }
    
private:
    HEGraph& g;
};

/// \brief From a voronoi-diagram, generate offset curve(s).
class MedialAxis {
public:
    MedialAxis(HEGraph& gi): g(gi) {
        medial_filter f(g);
        medial_filter flt(g);
        g.filter_graph(flt);
    }
private:
    MedialAxis(); // don't use.
    HEGraph& g; // original graph
    

};


} // end namespace

// end file polygon_interior.hpp
