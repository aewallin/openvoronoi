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

#include "graph.hpp"
#include "common/numeric.hpp"
#include "site.hpp"
#include "filter.hpp"

namespace ovd
{

/// \brief Filter for retaining the medial-axis of a voronoi diagram
///
/// (approximate) medial-axis Filter.
/// marks the valid-property true for edges belonging to the medial axis
/// and false for other edges.
struct medial_axis_filter : public Filter {
    /// \param thr dot-product threshold used to decide whether the segments
    /// that connect to a given Edge are nearly parallel
    medial_axis_filter( double thr=0.8) :  _dot_product_threshold(thr) { }
    
    /// predicate that decides if an edge is to be included or not.
    bool operator()(const HEEdge& e) const {
        if ( (*g)[e].type == LINESITE || (*g)[e].type == NULLEDGE) 
            return true; // we keep linesites and nulledges
        if ( (*g)[e].type == SEPARATOR)
            return false; // separators are allways removed
            
        if (both_endpoints_positive(e)) // these are interior edges which we keep.
            return true;
        
        // this leaves us with edges where one end connects to the polygon (dist==0)
        // and the other end does not.
        // figure out the angle between the adjacent line-segments and decide based on the angle.
        if (segments_parallel(e))
            return false;

        return true; // otherwise we keep the edge
    }
private:
    /// return true if this is an internal edge, i.e. both endpoints have a nonzero clearance-disk radius 
    bool both_endpoints_positive(HEEdge e) const {
        HEVertex src = g->source(e);
        HEVertex trg = g->target(e);
        return ((*g)[src].dist()>0) && ((*g)[trg].dist()>0);
    }
    /// return true if the segments that connect to the given Edge are nearly parallel
    bool segments_parallel( HEEdge e ) const {
        HEVertex endp1 = find_endpoint(e);
        HEVertex endp2 = find_endpoint( (*g)[e].twin );
        // find the segments
        HEEdge e1 = find_segment(endp1);
        HEEdge e2 = find_segment(endp2);
        e2 = (*g)[e2].twin; // this makes the edges oriented in the same direction 
        double dotprod = edge_dotprod(e1,e2);
        return dotprod >_dot_product_threshold;
    }
    
    /// \brief calculate the dot-product between unit vectors aligned along edges e1->e2
    ///
    /// since e1 and e2 are both line-sites the direction is easy to find
    /// FIXME: more code needed here for tangent calculation if we have arc-sites
    double edge_dotprod(HEEdge e1, HEEdge e2) const {
        HEVertex src1 = g->source(e1);
        HEVertex trg1 = g->target(e1);
        HEVertex src2 = g->source(e2);
        HEVertex trg2 = g->target(e2);
        Point sp1 = (*g)[src1].position;
        Point tp1 = (*g)[trg1].position;
        Point sp2 = (*g)[src2].position;
        Point tp2 = (*g)[trg2].position;
        
        Point dir1 = tp1-sp1;
        Point dir2 = tp2-sp2;
        dir1.normalize();
        dir2.normalize();
        return dir1.dot(dir2);
    }
    
    /// find the LineSite edge that connects to \a v
    HEEdge find_segment(HEVertex v) const {
        BOOST_FOREACH(HEEdge e, g->out_edges(v)) {
            if ( (*g)[e].type == LINESITE )
                return e;
        }
        assert(0);
        exit(-1);
        return HEEdge();
    }
    
    /// find an ::ENDPOINT vertex that connects to Edge e through a ::NULLEDGE at either the source or target of e.
    HEVertex find_endpoint(HEEdge e) const {
        HEEdge next = (*g)[e].next;
        HEEdge prev = g->previous_edge(e);
        HEVertex endp;
        if ( (*g)[next].type == NULLEDGE ) {
            endp = g->target(next);
            assert( (*g)[endp].type == ENDPOINT );
            
        } else if ( (*g)[prev].type == NULLEDGE ) {
            endp = g->source(prev);
            assert( (*g)[endp].type == ENDPOINT );
        } else {
            assert(0);
            exit(-1);
        }
        return endp;
    }

    double _dot_product_threshold; ///< a dot-product threshold in [0,1] for filtering out edges between nearly parallel LineSite segments
};


} // end namespace

// end file 
