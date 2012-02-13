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
        Point c2 = g[e].point(x); // current MIC center
        double r2 = x; // current MIC radius
        double w = (c2-c1).norm() + r2 - r1; // this is the cut-width
        return w-w_max; // error compared to desired cut-width
    }
private:
    HEGraph& g;
    double w_max; // desired cut-width
    HEEdge e;
    Point c1; // previous MIC center
    double r1; // previous MIC radius
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
        debug = false;
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
        if (debug) { std::cout << " start edge is: "; g.print_edge(current_edge); }
        new_branch=false;
        return out;
    }
    
    // get the next mic
    boost::python::list nxt_mic() {
        
        if ( current_edge == HEEdge() ) {
            if (debug) std::cout << "nxt_mic() end of operation. Nothing to do.\n";
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
        if (debug) std::cout << "nxt_mic() target width " << w_target << "\n";
        
        if ( w_target > max_width ) {
            // since moving to the target vertex would give too large cut-width
            // we search on the current edge for the next MIC
            if (debug) {  std::cout << " searching on the current edge "; g.print_edge(current_edge); }
            // find a point on the current edge
            double next_radius = find_next_radius();
            return output_next_mic(next_radius, new_branch);
        } else {
            // moving to the target edge does not give a cut-width that is large enough
            // we need to find a new edge for the next MIC
            
            // mark edge DONE. this means we have machined all MICs on this edge.
            mark_done(current_edge);
            if (debug) std::cout << "Finding new edge !\n";
            
            // check for end-of-branch and output a final MIC at the end of a branch.
            // output a final MIC at the end of the branch
            // so that no uncut-material will remain
            /*
            EdgeVector out_edges = find_out_edges();
            std::cout << " out_edges: " << out_edges.size() << " \n";
            if ( out_edges.empty() ) { // end of branch
                std::cout << " end-of-branch MIC at " << g[g.target(current_edge)].index << " \n";
                std::cout << " current_radius = " << current_radius << "\n";
                std::cout << " target radius =  " << g[ g.target(current_edge) ].dist() << "\n";

            }*/
            
            bool end_branch_mic=false;
            boost::tie( current_edge, end_branch_mic) = find_next_edge(); // move to the next edge
            if ( current_edge == HEEdge() ) { // invalid edge marks end of operation
                if (debug)  std::cout << "nxt_mic() end of operation.\n";
                return boost::python::list();
            }
            
            if ( end_branch_mic )
                return output_next_mic(current_radius, true);
            
            double next_radius = find_next_radius();
            if (new_branch) {
                new_branch=false;
                return output_next_mic(next_radius, true);
            } else {
                return output_next_mic(next_radius, false);
            }
            
        }
        
    }

    // pop an unvisited edge from the stack
    // or end-of-operation if the stack is empty
    HEEdge find_next_branch() {
        if (unvisited.empty() ) {
            if (debug) std::cout << "find_next_branch(): no un-machined branches. end operation.\n";
            return HEEdge();
        } else {
            branch_point out = unvisited.top();
            //debug = true;
            if (debug) { std::cout << "find_next_branch(): next branch is "; g.print_edge(out.next_edge); }
            //debug = false;
            unvisited.pop();
            previous_branch_center = current_center;
            previous_branch_radius = current_radius;
            
            current_center = out.current_center;
            current_radius = out.current_radius;
            new_branch = true;
            return out.next_edge;
        }
    }
    
    EdgeVector find_out_edges() {
        HEVertex trg = g.target(current_edge);
        EdgeVector out_edges;
        BOOST_FOREACH(HEEdge e, g.out_edge_itr(trg) ) {
            if ( e != g[current_edge].twin && g[e].valid && g[e].type != NULLEDGE ) {
                out_edges.push_back(e);
            }
        }
        return out_edges;
    }
    
    std::pair<HEEdge,bool> find_next_edge() {
        //HEVertex trg = g.target(current_edge);
        EdgeVector out_edges = find_out_edges();
        //BOOST_FOREACH(HEEdge e, g.out_edge_itr(trg) ) {
        //    if ( e != g[current_edge].twin && g[e].valid && g[e].type != NULLEDGE ) {
        //        out_edges.push_back(e);
        //    }
        //}
        //std::cout << "find_next_edge(): " << out_edges.size() << " potential next-edges\n";
        if (out_edges.empty() ) {
            if (debug) std::cout << "find_next_edge(): no out_edges. end of branch.\n";
            
            if ( current_radius > g[ g.target(current_edge) ].dist() ) {
                //std::cout << " end-of-branch MIC ! \n";
                current_radius = g[ g.target(current_edge) ].dist();
                return std::make_pair(current_edge,true);
                // how do we force an output here?
                //return output_next_mic( current_radius , false);
            }
            
            HEEdge e = find_next_branch();
            if (e==HEEdge())
                return std::make_pair(e,false); // this is end-of-operation.
                
            if (has_next_radius(e)) {
                if (debug) { std::cout << "find_next_edge(): next-branch is: "; g.print_edge(e); }
                return std::make_pair(e,false);
            } else {
                if (debug) { std::cout << "find_next_edge(): next-branch, but not valid\n"; }
                mark_done( e );
                current_edge = e; // jump to the next edge
                return find_next_edge(); // and try to see if that edge is valid.
            } 
        } else if ( out_edges.size() == 1 ) {
            // only one choice for the next edge
            if (has_next_radius(out_edges[0])) {
                if (debug) { std::cout << "find_next_edge(): only one out-edge: "; g.print_edge(out_edges[0]); }
                return std::make_pair(out_edges[0],false);
            } else {
                if (debug) std::cout << "find_next_edge(): one out-edge, but not valid\n"; 
                mark_done( out_edges[0] );
                current_edge = out_edges[0]; // jump to the next edge
                return find_next_edge(); // and try to see if that edge is valid.
            } 
        } else if (out_edges.size() == 2 ) {
            // two choices for the next edge
            // FIXME: some smarter way of selecting next-edge
            unvisited.push( branch_point(current_center, current_radius, out_edges[1] ) );
            
            if (has_next_radius(out_edges[0])) {
                if (debug) { std::cout << "find_next_edge(): two out-edges, returning first: "; g.print_edge(out_edges[0]); }
                return std::make_pair(out_edges[0],false);
            } else {
                mark_done( out_edges[0] ); //].done = true;
                current_edge = out_edges[0]; // jump to the next edge
                return find_next_edge(); // and try to see if that edge is valid.
            } 
            
            //return out_edges[0];
        } else {
            if (debug) std::cout << "find_next_edge(): too many out-edges. ERROR.\n";
            exit(-1);
            return std::make_pair(HEEdge(),false);
        }
    }
    
    void mark_done(HEEdge e) {
        edge_data[ e ].done = true;
        edge_data[ g[e].twin ].done = true;
    }
    
    // does HEEdge e have the next MIC we want ?
    bool has_next_radius(HEEdge e) {
        // check if the edge e is one where we continue
        Point c1 = g[current_edge].point(current_radius); 
        double r1 = current_radius;
        double r2 = g[ g.target(e) ].dist();
        Point c2 = g[e].point(r2);
        //double r2 = target_radius;
        double w_target = ( c2-c1 ).norm() + r2 - r1;
        if (debug) std::cout << " has_next_radius() " << w_target << "\n";
        if (w_target<=0)
            return false;
            
        if ( w_target > max_width )
            return true;
        else
            return false;
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

        double trg_err = t(target_radius);
        double cur_err = t(current_radius);
        if ( debug ||  ( !(trg_err*cur_err < 0) ) ) {
            std::cout << " current rad = " << current_radius << "\n";
            std::cout << " target rad = " << target_radius << "\n";
            std::cout << " error at current = " << t(current_radius) << "\n";
            std::cout << " error at target = " << t(target_radius) << "\n";
        }
        double min_r = std::min(current_radius, target_radius);
        double max_r = std::max(current_radius, target_radius);
        Result r1 = boost::math::tools::toms748_solve(t, min_r, max_r, tol, max_iter);
        return r1.first;
    }
    
    // output the next MIC, for processing by a downstream algorithm
    // that lays out the pattern: lead-out, rapid, lead-in, bi-tangent, cut-arc, bi-tangent
    boost::python::list output_next_mic(double next_radius, bool branch) {
        boost::python::list out;
        if (debug) std::cout << " next_radius = " << next_radius << "\n";
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
        out.append(tang1); // bi-tangent points
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
    
private:
    bool debug;
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
