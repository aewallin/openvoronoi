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

#include "medial_axis_walk.hpp"

namespace ovd
{


/// \brief add the given edge to the current list of edges.
///
/// for line-edges we add only two endpoints
/// for parabolic edges we add many points
void MedialAxisWalk::append_edge(MedialChain& chain, HEEdge edge)  {
    MedialPointList point_list; // the endpoints of each edge
    HEVertex v1 = g.source( edge );
    HEVertex v2 = g.target( edge );
    // these edge-types are drawn as a single line from source to target.
    if ( (g[edge].type == LINELINE)  || (g[edge].type == PARA_LINELINE) ) {
        MedialPoint pt1( g[v1].position, g[v1].dist() );
        MedialPoint pt2( g[v2].position, g[v2].dist() );
        point_list.push_back(pt1);
        point_list.push_back(pt2);
    } else if ( (g[edge].type == PARABOLA) || (g[edge].type == LINE) ) { // these edge-types are drawn as polylines with _edge_points number of points
        double t_src = g[v1].dist();
        double t_trg = g[v2].dist();
        double t_min = std::min(t_src,t_trg);
        double t_max = std::max(t_src,t_trg);
        
        for (int n=0;n< _edge_points;n++) {
            double t(0);
            if (t_src<=t_trg) // increasing t-value
                t = t_min + ((t_max-t_min)/numeric::sq(_edge_points-1))*numeric::sq(n); // NOTE: quadratic t-dependece. More points at smaller t.
            else if (t_src>t_trg) { // decreasing t-value
                int m = _edge_points-1-n; // m goes from (N-1)...0   as n goes from 0...(N-1)
                t = t_min + ((t_max-t_min)/numeric::sq(_edge_points-1))*numeric::sq(m);
            }
            Point p = g[edge].point(t);
            MedialPoint pt( p, t );
            point_list.push_back(pt);
        }
    }
    chain.push_back( point_list );
}

/// we are at target(e). find the next suitable edge.
/// \return true if a next-edge was found, false otherwise.
bool MedialAxisWalk::next_edge(HEEdge e, HEEdge& next) {
    HEVertex trg = g.target(e);
    EdgeVector out_edges = g.out_edges(trg);
    
    std::vector<HEEdge> valid_edges; // find all valid edges
    BOOST_FOREACH( HEEdge oe, out_edges) {
        if ( valid_next_edge(oe) ) {
            valid_edges.push_back(oe);
        }
    }
    if (!valid_edges.empty() ) {
        next = valid_edges[0]; // return the first valid one
        return true;
    }
    return false; 
}

    
/// set edge and its twin invalid
void MedialAxisWalk::set_invalid(HEEdge e) {
    g[e].valid = false;
    if (g[e].twin != HEEdge()) {
        g[ g[e].twin ].valid = false;
    }
}

/// loop through all edges and find an edge where we can start
/// valid edges have a source-vertex with exactly one valid out-edge.
bool MedialAxisWalk::find_start_edge(HEEdge& start) {
    start = HEEdge();
    BOOST_FOREACH(HEEdge e, g.edges() ) { 
        if ( valid_next_edge(e) ) {
            if (degree_one_source(e)) { // for an o-shaped medial axis, this doesn't find a start-edge!
                start = e;
                return true;
            }
        }
    }
    
    // if we get here, there are no "dangling" edges where we can start.
    // but there might be an o-shaped feature which is un-machined.
    EdgeVector valid_edges;
    BOOST_FOREACH(HEEdge e, g.edges() ) { 
        if ( valid_next_edge(e) ) {
            valid_edges.push_back(e);
        }
    }
    // just return the first valid edge. (FIXME: better to start with a shallow cut, i.e. return the edge with the smallest clearance-disk?
    if (!valid_edges.empty()) { // for an o-shaped medial axis, this doesn't find a start-edge!
        start = valid_edges[0];
        return true;
    }
    
    
    return false;
}

/// check if the source of the edge is a valid starting-point for a path
bool MedialAxisWalk::degree_one_source(HEEdge e) {
    HEVertex src = g.source(e);
    EdgeVector out_edges = g.out_edges(src);
    int count(0);
    BOOST_FOREACH( HEEdge oe, out_edges) {
        if ( valid_next_edge(oe) ) {
            count++;
        }
    }
    return (count==1);
}

/// find start-edgem then walk
void MedialAxisWalk::do_walk() {
    out = MedialChainList();
    HEEdge start = HEEdge();
    while( find_start_edge(start) ) { // find a suitable start-edge
        medial_axis_walk(start); // from the start-edge, walk as far as possible
    }
}

/// we can follow an edge if it is valid, and not a ::LINESITE or ::NULLEDGE
bool MedialAxisWalk::valid_next_edge(HEEdge e) {
    return ( (g[e].type != LINESITE) && (g[e].type !=NULLEDGE) && (g[e].valid) );
}
    
/// start at source of Edge start, and walk as far as possible
void MedialAxisWalk::medial_axis_walk(HEEdge start) {
    // begin chain with start.
    HEEdge next = start; // why does = HEEdge() cause Wuninitialized ?
    MedialChain chain;
    append_edge(chain, start);
    set_invalid(start);
    while (next_edge(start, next)  ) {
        assert( g.target( start ) == g.source( next ) );             
        append_edge(chain, next);
        start=next; 
        set_invalid(start);
    }
    // end chain
    out.push_back( chain );
}
    
}
