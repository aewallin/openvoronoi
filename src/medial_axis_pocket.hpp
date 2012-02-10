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
#include <stack>

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


// branch-data when we backtract to machine an un-machined branch
struct branch_point {
    branch_point(Point p, double r, HEEdge e) {
        current_center = p;
        current_radius = r;
        next_edge = e;
    }
    Point current_center;
    double current_radius;
    HEEdge next_edge;
};

// error functor for numerically finding the next MIC
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
            if ( g[e].valid && g[e].type != LINESITE && g[e].type != NULLEDGE) {
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
        
        // find the vertex with the maximum radius mic
        double max_mic_radius(-1);
        Point max_mic_pos(0,0);
        HEVertex max_mic_vertex = HEVertex();
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
        current_center = max_mic_pos;
        previous_branch_center = max_mic_pos;
        previous_branch_radius = max_mic_radius;
        
        // find the edge on which we start machining.
        // stash the other out-edges for visiting later
        double max_adj_radius(-1);
        BOOST_FOREACH( HEEdge e, g.out_edge_itr(max_mic_vertex) ) {
            //std::cout << "potential start edge: "; g.print_edge(e);
            if ( g[ g.target(e) ].dist() > max_adj_radius ) {
                max_adj_radius = g[ g.target(e) ].dist();
                current_edge = e;
            }
        }
        BOOST_FOREACH( HEEdge e, g.out_edge_itr(max_mic_vertex) ) {
            //std::cout << "potential start edge: "; g.print_edge(e);
            if ( e != current_edge ) {
                unvisited.push( branch_point(current_center, current_radius, e ) );
            }
        }
        std::cout << " start edge is: "; g.print_edge(current_edge);
        new_branch=false;
        return out;
    }
    
    // get the next mic
    boost::python::list nxt_mic() {
        boost::python::list out;
        if ( current_edge == HEEdge() ) {
            std::cout << "nxt_mic() end of operation. Nothing to do.\n";
            return boost::python::list();
        }
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
        std::cout << "nxt_mic() target width " << w_target << "\n";
        if ( w_target > max_width ) {
            std::cout << " searching on the current edge "; g.print_edge(current_edge);
            // find a point on the current edge
            double next_radius = find_next_radius();
            return output_next_mic(next_radius, new_branch);
        } else {
            // mark edge DONE. this means we have machined all MICs on this edge.
            edge_data[current_edge].done = true;
            edge_data[ g[current_edge].twin ].done = true;
            std::cout << "Finding new edge !\n";
            current_edge = find_next_edge();              // move to the next edge
            if ( current_edge == HEEdge() ) { // invalid edge marks end of operation
                std::cout << "nxt_mic() end of operation.\n";
                return boost::python::list();
            }
            double next_radius = find_next_radius();
            if (new_branch) {
                new_branch=false;
                return output_next_mic(next_radius, true);
            } else {
                return output_next_mic(next_radius, false);
            }
            //return out;
        }
        
    }
    boost::python::list output_next_mic(double next_radius, bool branch) {
        boost::python::list out;
        std::cout << " next_radius = " << next_radius << "\n";
        out.append( g[current_edge].point(next_radius) );
        out.append( next_radius );
        Point c1 = current_center;
        Point c2 = g[current_edge].point(next_radius);
        double r1 = current_radius;
        double r2 = next_radius;
        current_radius = next_radius;
        current_center = g[current_edge].point(next_radius);
        
        // find the bi-tangents and return them too.
        // see voronoi_bisectors.py
        double detM = c1.x*c2.y - c2.x*c1.y;
        double m = ( c1.y-c2.y ) / detM;
        double p = ( c2.x-c1.x ) / detM;
        double n = ( c2.y*r1 - c1.y*r2 ) / detM;
        double q = ( c1.x*r2 - c2.x*r1 ) / detM;
        std::vector<double> roots = numeric::quadratic_roots( m*m+p*p, 2*(m*n+p*q),  n*n+q*q-1);
        //BOOST_FOREACH(double r, roots) {
        //    std::cout << " root " << r << "\n";
        //}
        // bi-tangent lines are now
        // ax +  by + c = 0
        // with
        // C = root
        // A = m*C+n
        // B = p*C+q
        double lc1 = roots[0];
        double a1 = m*lc1+n;
        double b1 = p*lc1+q;
        double lc2 = roots[1];
        double a2 = m*lc2+n;
        double b2 = p*lc2+q;
        // the bi-tangent points are given by
        Point tang1 = c1 - r1*Point( a1, b1 );
        Point tang2 = c1 - r1*Point( a2, b2 );
        Point tang3 = c2 - r2*Point( a1, b1 );
        Point tang4 = c2 - r2*Point( a2, b2 );
        out.append(tang1);
        out.append(tang2);
        out.append(tang3);
        out.append(tang4);
        out.append(c1); // previous MIC center
        out.append(r1); // previous MIC radius
        out.append(branch); // true/false flag for new branch
        out.append(previous_branch_center);
        out.append(previous_branch_radius);
        
        return out;
    }
    
    HEEdge find_next_branch() {
        if (unvisited.empty() ) {
            std::cout << "find_next_branch(): no un-machined branches. end operation.\n";
            return HEEdge();
        } else {
            branch_point out = unvisited.top();
            std::cout << "find_next_branch(): next branch is "; g.print_edge(out.next_edge);
            unvisited.pop();
            previous_branch_center = current_center;
            previous_branch_radius = current_radius;
            
            current_center = out.current_center;
            current_radius = out.current_radius;
            new_branch = true;
            return out.next_edge;
        }
    }
    
    HEEdge find_next_edge() {
        HEVertex trg = g.target(current_edge);
        EdgeVector out_edges;
        BOOST_FOREACH(HEEdge e, g.out_edge_itr(trg) ) {
            if ( e != g[current_edge].twin && g[e].valid && g[e].type != NULLEDGE ) {
                out_edges.push_back(e);
            }
        }
        std::cout << "find_next_edge(): " << out_edges.size() << " potential next-edges\n";
        if (out_edges.empty() ) {
            std::cout << "find_next_edge(): no out_edges. end of branch.\n";
            return find_next_branch();
        } else if ( out_edges.size() == 1 ) {
            std::cout << "find_next_edge(): only one out-edge: "; g.print_edge(out_edges[0]);
            return out_edges[0];
        } else if (out_edges.size() == 2 ) {
            std::cout << "find_next_edge(): two out-edges, returning first: "; g.print_edge(out_edges[0]);
            // FIXME: some smarter way of selecting next-edge
            unvisited.push( branch_point(current_center, current_radius, out_edges[1] ) );
            return out_edges[0];
        } else {
            std::cout << "find_next_edge(): too many out-edges. ERROR.\n";
            exit(-1);
            return HEEdge();
        }
    }
    
    // on the current edge, move from current_radius towards target_radius
    // and find a radius-value that satisfies the cut-width constraint.
    double find_next_radius() {
        // HEGraph& gi, double wmax, HEEdge search_edge, Point cen1, double rad1)
        CutWidthError t(g,max_width, current_edge, current_center, current_radius);
        typedef std::pair<double, double> Result;
        boost::uintmax_t max_iter=500;
        boost::math::tools::eps_tolerance<double> tol(30);
        double target_radius = g[ g.target(current_edge) ].dist();
        
        //std::cout << " error at current = " << t(current_radius) << "\n";
        //std::cout << " error at target = " << t(target_radius) << "\n";
        double trg_err = t(target_radius);
        double cur_err = t(current_radius);
        if ( !(trg_err*cur_err < 0) ) {
            std::cout << " current rad = " << current_radius << "\n";
            std::cout << " target rad = " << target_radius << "\n";
            std::cout << " error at current = " << t(current_radius) << "\n";
            std::cout << " error at target = " << t(target_radius) << "\n";
            
            
        }
        
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
    std::stack<branch_point> unvisited;
    
    HEEdge current_edge;
    double current_radius;
    Point current_center;

    // flag for indicating new branch
    bool new_branch;
    Point previous_branch_center;
    double previous_branch_radius;
    
    // the max cutting-width
    double max_width; 
    
};

} // end namespace

// end file 
