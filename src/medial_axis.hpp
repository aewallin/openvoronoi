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

namespace ovd
{

// (approximate) medial-axis filter
// marks the valid-property true for edges belonging to the medial axis
// and false for other edges.
struct medial_filter {
    medial_filter(HEGraph& gi, double thr=0.8) : g(gi) , _dot_product_threshold(thr) { }
    bool operator()(const HEEdge& e) const {
        if (g[e].type == LINESITE || g[e].type == NULLEDGE) 
            return true; // we keep linesites and nulledges
        if (g[e].type == SEPARATOR)
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
    // return true if this is an internal edge, i.e. both endpoints have a nonzero clearance-disk radius 
    bool both_endpoints_positive(HEEdge e) const {
        HEVertex src = g.source(e);
        HEVertex trg = g.target(e);
        return (g[src].dist()>0) && (g[trg].dist()>0);
    }
    // return true if the segments that connect to the given Edge are nearly parallel
    bool segments_parallel( HEEdge e ) const {
        HEVertex endp1 = find_endpoint(e);
        HEVertex endp2 = find_endpoint( g[e].twin );
        // find the segments
        HEEdge e1 = find_segment(endp1);
        HEEdge e2 = find_segment(endp2);
        e2 = g[e2].twin; // this makes the edges oriented in the same direction 
        double dotprod = edge_dotprod(e1,e2);
        return fabs(dotprod)>_dot_product_threshold;
    }
    
    // calculate the dot-product between unit vectors aligned along edges e1->e2
    // since e1 and e2 are both line-sites the direction is easy to find
    // FIXME: more code needed here for tangent calculation if we have arc-sites
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
    
    // find the linesite edge that connects to v.
    HEEdge find_segment(HEVertex v) const {
        BOOST_FOREACH(HEEdge e, g.out_edges(v)) {
            if ( g[e].type == LINESITE )
                return e;
        }
        assert(0);
        exit(-1);
        return HEEdge();
    }
    
    // find an ENDPOINT vertex that connects to Edge e through a NULLEDGE at either the source or target of e.
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
    double _dot_product_threshold;
};

/// \brief From a voronoi-diagram, generate offset curve(s).
/// argument dp_thr: dot-product threshold used to decide whether the segments
/// that connect to a given Edge are nearly parallel
class MedialAxis {
public:
    MedialAxis(HEGraph& gi, double dp_thr=0.8): g(gi) {
        medial_filter flt(g, dp_thr);
        g.filter_graph(flt);
    }
private:
    MedialAxis(); // don't use.
    HEGraph& g; // original graph
};


/// \brief Medial-axis point and associated clearance-disc radius.
struct MedialPoint {
    Point p;
    double clearance_radius;
    MedialPoint(Point pi, double r): p(pi), clearance_radius(r) {}
};
typedef std::list<MedialPoint> MedialPointList;
typedef std::list<MedialPointList> MedialChain;
typedef std::list<MedialChain> MedialChainList;

// FIXME: MedialAxisWalk could probably be optimized to minimize rapid-traverses

/// \brief Walk along the medial-axis edges of a voronoi-diagram.
/// When we want a toolpath along the medial axis we use this class
/// to walk along the "valid" edges which are left in the diagram.
///
/// argument edge_pts: number of points to subdivide parabolas.
///
/// Algorithm:
/// first find one valid edge that has a degree-1 vertex (i.e. a suitable start point for the path)
/// - if there's only one choice for the next edge, go there
/// - if there are two choices, take one of the choices
/// when done, find another valid start-edge.
class MedialAxisWalk {
public:
    MedialAxisWalk(HEGraph& gi, int edge_pts = 20): g(gi), _edge_points(edge_pts) {}

    void do_walk() {
        out = MedialChainList();
        HEEdge start;
        while( find_start_edge(start) ) { // find a suitable start-edge
            medial_axis_walk(start); // from the start-edge, walk as far as possible
        }
    }

    MedialChainList walk() {
        do_walk();
        return out;
    }
    
    // start at source of Edge start, and walk as far as possible
    void medial_axis_walk(HEEdge start) {
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
    
    // add the given edge to the current list of edges.
    // for line-edges we add only two endpoints
    // for parabolic edges we add many points
    void append_edge(MedialChain& chain, HEEdge edge)  {
        MedialPointList point_list; // the endpoints of each edge
        HEVertex v1 = g.source( edge );
        HEVertex v2 = g.target( edge );
        // these edge-types are drawn as a single line from source to target.
        if (   (g[edge].type == LINELINE)  || (g[edge].type == PARA_LINELINE)) {
            MedialPoint pt1( g[v1].position, g[v1].dist() );
            MedialPoint pt2( g[v2].position, g[v2].dist() );
            point_list.push_back(pt1);
            point_list.push_back(pt2);
        } else if ( (g[edge].type == PARABOLA) || (g[edge].type == LINE) ) { // these edge-types are drawn as polylines with _edge_points number of points
            double t_src = g[v1].dist();
            double t_trg = g[v2].dist();
            double t_min = std::min(t_src,t_trg);
            double t_max = std::max(t_src,t_trg);
            //_edge_points: number of points to subdivide parabolas.
            
            for (int n=0;n< _edge_points;n++) {
                double t(0);
                if (t_src<=t_trg) // increasing t-value
                    t = t_min + ((t_max-t_min)/numeric::sq(_edge_points-1))*numeric::sq(n); // NOTE: quadratic t-dependece. More points at smaller t.
                else if (t_trg<t_src) { // decreasing t-value
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
    // we are at target(e). find the next suitable edge.
    // return true if a next-edge was found, false otherwise.
    bool next_edge(HEEdge e, HEEdge& next) {
        HEVertex trg = g.target(e);
        EdgeVector out_edges = g.out_edges(trg);
        std::vector<HEEdge> valid_edges;
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
    
    void set_invalid(HEEdge e) {
        g[e].valid = false;
        if (g[e].twin != HEEdge()) {
            g[ g[e].twin ].valid = false;
        }
    }
    // loop through all edges and find an edge where we can start
    // valid edges have a source-vertex with exactly one valid out-edge.
    bool find_start_edge(HEEdge& start) {
        BOOST_FOREACH(HEEdge e, g.edges() ) { 
            if ( valid_next_edge(e) ) {
                if (degree_one_source(e)) {
                    start = e;
                    return true;
                }
            }
        }
        return false;
    }
    // we can follow an edge if it is valid, and not a LINESITE or NULLEDGE
    bool valid_next_edge(HEEdge e) {
        return ( (g[e].type != LINESITE) && (g[e].type !=NULLEDGE) && (g[e].valid) );
    }
    // check if the source of the edge is a valid starting-point for a path
    bool degree_one_source(HEEdge e) {
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
protected:
    MedialChainList out;
private:
    MedialAxisWalk(); // don't use.
    HEGraph& g; // original graph
    int _edge_points; // number of points to subdivide parabolas.

};

} // end namespace

// end file medial_axis.hpp
