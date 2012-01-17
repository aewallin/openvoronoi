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


// walk along the "valid" edges which are left in the diagram 
// first find one valid edge that has a degree-1 vertex (i.e. a suitable start point for the path)
// -- if there's only one choice for the next edge, go there
// -- if there are two choices, take one of the choices
// when done, find another valid start-edge

class MedialAxisWalk {
public:
    MedialAxisWalk(HEGraph& gi): g(gi) {
        //medial_filter f(g);
        //medial_filter flt(g);
        //g.filter_graph(flt);
        
    }
    boost::python::list walk() {
        out = boost::python::list();
        
        HEEdge start;
        while( find_start_edge(start) ) { // find a suitable start-edge
            ma_walk(start); // from the start-edge, walk as far as possible
        }
        
        //std::cout << " start edge is "; g.print_edge(start);
        return out;
    }
    
    void ma_walk(HEEdge start) {
        // start at source of start, and walk as far as possible
        // begin chain with start.
        HEEdge next;
        boost::python::list chain;
        append_edge(chain, start);
        set_invalid(start);
        //bool flag;
        while (next_edge(start, next)  ) {
            
            if ( g.target( start ) == g.source( next ) ) { 
                append_edge(chain, next);
                start=next; // g[next].twin; //next;
            } else if ( g.target( start ) == g.target( next ) ) {
                exit(-1);
            } else {
                exit(-1);
            }
            
            //start = next;
            set_invalid(start);
        }
        // end chain
        out.append( chain );
    }
    
    void append_edge(boost::python::list& list, HEEdge edge)  {
        //boost::python::list edge_data;
        
        boost::python::list point_list; // the endpoints of each edge
        HEVertex v1 = g.source( edge );
        HEVertex v2 = g.target( edge );
        // these edge-types are drawn as a single line from source to target.
        if (  (g[edge].type == LINE) || (g[edge].type == LINELINE)  || (g[edge].type == PARA_LINELINE)) {
            boost::python::list pt1;
            pt1.append( g[v1].position ); pt1.append( g[v1].dist() );
            point_list.append(pt1);
            boost::python::list pt2;
            pt2.append( g[v2].position ); pt2.append( g[v2].dist() );
            point_list.append(pt2);
            //point_list.append( g[v2].position );
        } else if ( g[edge].type == PARABOLA ) { // these edge-types are drawn as polylines with edge_points number of points
            double t_src = g[v1].dist();
            double t_trg = g[v2].dist();
            double t_min = std::min(t_src,t_trg);
            double t_max = std::max(t_src,t_trg);
            int _edge_points= 20; // number of points to subdivide parabolas
            
            for (int n=0;n< _edge_points;n++) {
                //double t = t_min + n*(t_max-t_min)/(_edge_points-1); // linear
                double t;
                if (t_src<=t_trg) // increasing t-value
                    t = t_min + ((t_max-t_min)/sq(_edge_points-1))*sq(n);
                else if (t_trg<t_src) // decreasing t-value
                    t = t_max - ((t_max-t_min)/sq(_edge_points))*sq(n);
                else
                    exit(-1);
                Point p = g[edge].point(t);
                boost::python::list pt;
                pt.append( p); pt.append( t );
                point_list.append(pt);
            }
        }
        //edge_data.append( point_list );
        //edge_data.append( g[edge].type );
        //edge_data.append( g[v1].status ); // source status
        //edge_data.append( g[v2].status ); // target status
        list.append( point_list );
    }
    
    bool next_edge(HEEdge e, HEEdge& next) {
        HEVertex trg = g.target(e);
        EdgeVector out_edges = g.out_edges(trg);
        //int count(0);
        std::vector<HEEdge> valid_edges;
        BOOST_FOREACH( HEEdge oe, out_edges) {
            if ( valid_next_edge(oe) ) {
                valid_edges.push_back(oe);
                //count++;
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
        g[ g[e].twin ].valid = false;
    }
    bool find_start_edge(HEEdge& start) {
        BOOST_FOREACH(HEEdge e, g.edges() ) { // loop through all edges and find an edge where we can start
            if ( valid_next_edge(e) ) {
                if (degree_one_source(e)) {
                    start = e;
                    return true;
                } 
                /*else if ( degree_one_endpoint( g[e].twin ) ) {
                    start = g[e].twin;
                    return true;
                }*/
            }
        }
        return false;
    }
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
private:
    boost::python::list out;
    MedialAxisWalk(); // don't use.
    HEGraph& g; // original graph

};

} // end namespace

// end file polygon_interior.hpp
