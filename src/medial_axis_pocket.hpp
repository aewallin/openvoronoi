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
#include <boost/math/tools/roots.hpp>

#include "graph.hpp"
#include "common/numeric.hpp"
#include "site.hpp"

namespace ovd
{

// keep track of data for each edge here
struct edata {
    edata() {
        done = false;
    }
    bool done;
};

typedef std::pair<HEEdge , edata> Edata;

class CutWidthError  {
public:
    CutWidthError(HEGraph& gi, double wmax, HEEdge search_edge, Point cen1, double rad1) 
    : g(gi), w_max(wmax), e(search_edge), c1(cen1), r1(rad1) {
    }
    double operator()(const double x) {
        // w_max = | c2 - c1 | + r2 - r1
        Point c2 = g[e].point(x);
        double r2 = x;
        double w = (c2-c1).norm() + r2 - r1;
        return w-w_max;
    }
private:
    HEGraph& g;
    double w_max;
    HEEdge e;
    Point c1;
    double r1;
    
    
};

// experimental medial-axis pocketing
class medial_axis_pocket {
public:
    medial_axis_pocket(HEGraph& gi): g(gi) {
        BOOST_FOREACH( HEEdge e, g.edges() ) {
            if ( g[e].valid && g[e].type != LINESITE ) {
                ma_edges.push_back(e);
                edata ed;
                edge_data.insert( Edata(e, ed ) );
            }
        }
        current_edge = HEEdge();
        max_width = 0.05;
    }
    void set_width(double w) {max_width=w;}
    
    // return the maximum mic
    // clear this with a spiral-path before proceeding.
    boost::python::list max_mic() {
        boost::python::list out;
        
        double max_mic_radius(-1);
        Point max_mic_pos(0,0);
        HEVertex max_mic_vertex;
        BOOST_FOREACH( HEEdge e, ma_edges ) {
            HEVertex src = g.source(e);
            if ( g[src].dist() > max_mic_radius ) {
                max_mic_radius = g[src].dist();
                max_mic_pos = g[src].position; 
                max_mic_vertex = src;
            }
        }
        out.append(max_mic_pos);
        out.append(max_mic_radius);
        current_radius = max_mic_radius;
        
        // find the edge on which we start machining.
        double max_adj_radius(-1);
        BOOST_FOREACH( HEEdge e, g.out_edge_itr(max_mic_vertex) ) {
            //std::cout << "potential start edge: "; g.print_edge(e);
            if ( g[ g.target(e) ].dist() > max_adj_radius ) {
                max_adj_radius = g[ g.target(e) ].dist();
                current_edge = e;
            }
        }
        std::cout << " start edge is: "; g.print_edge(current_edge);
        return out;
    }
    
    // get the next mic
    boost::python::list nxt_mic() {
        boost::python::list out;
        // find a point on current-edge so that we get the desired 
        // cut-width
        //  w_max = | c2 - c1 | + r2 - r1
        
        // we allways move from source to target.
        double target_radius = g[ g.target(current_edge) ].dist();
        Point c1 = g[current_edge].point(current_radius); 
        double r1 = current_radius;
        Point c2 = g[current_edge].point(target_radius);
        double r2 = target_radius;
        double w_target = ( c2-c1 ).norm() + r2 - r1;
        std::cout << " target width " << w_target << "\n";
        if ( w_target > max_width ) {
            std::cout << " searching on the current ege "; g.print_edge(current_edge);
            // find a point on the current edge
            double next_radius = find_next_radius();
            std::cout << " next_radius = " << next_radius << "\n";
            
            out.append( g[current_edge].point(next_radius) );
            out.append( next_radius );
            current_radius = next_radius;
            
        } else {
            std::cout << "Finding new edge !\n";// g.print_edge(current_edge);
            // move to the next edge
            exit(-1);
        }
        return out;
    }
    
    // on the current edge, move from current_radius towards target_radius
    // and find a radius-value that satisfies the cut-width constraint.
    double find_next_radius() {
        // HEGraph& gi, double wmax, HEEdge search_edge, Point cen1, double rad1)
        CutWidthError t(g,max_width, current_edge, g[current_edge].point(current_radius), current_radius);
        typedef std::pair<double, double> Result;
        boost::uintmax_t max_iter=500;
        boost::math::tools::eps_tolerance<double> tol(30);
        double target_radius = g[ g.target(current_edge) ].dist();
        
        //std::cout << " error at current = " << t(current_radius) << "\n";
        //std::cout << " error at target = " << t(target_radius) << "\n";
        
        double min_r = std::min(current_radius, target_radius);
        double max_r = std::max(current_radius, target_radius);
        Result r1 = boost::math::tools::toms748_solve(t, min_r, max_r, tol, max_iter);
        //std::cout << "root bracketed: [ " << r1.first << " , " << r1.second << " ]" << std::endl;
        //std::cout << "f("<< r1.first << ")=" << t(r1.first) << std::endl;
        //std::cout << "f("<< r1.second << ")=" << t(r1.second) << std::endl;
        //std::cout << "max_iter=" << max_iter << std::endl;
        return r1.first;
    }
    
private:
    std::vector<HEEdge> ma_edges; // the edges of the medial-axis
    std::map<HEEdge, edata> edge_data;
    HEGraph& g; // VD graph

    HEEdge current_edge;
    double current_radius;
    
    double max_width;
};

} // end namespace

// end file 
