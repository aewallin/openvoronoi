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

#include "medial_axis_pocket.hpp"

namespace ovd
{

/// the input graph should be a voronoi diagram passed through medial_axis_filter
medial_axis_pocket::medial_axis_pocket(HEGraph& gi): g(gi) {
    BOOST_FOREACH( HEEdge e, g.edges() ) {
        if ( g[e].valid && 
             g[e].type != LINESITE && 
             g[e].type != NULLEDGE && 
             g[e].type != OUTEDGE   ) {
            ma_edges.push_back(e);
            edata ed;
            edge_data.insert( Edata(e, ed ) );
        }
    }
    current_edge = HEEdge();
    max_width = 0.05;
    debug = false;
}


/// set the maximum cut width
void medial_axis_pocket::set_width(double w) {max_width=w;}

/// run the algorithm
void medial_axis_pocket::run() {
    find_initial_mic();
    while (find_next_mic()) {}
    if (debug) std::cout << "medial_axis_pocket::run() done. generated " << mic_list.size() << " MICs \n";
}

void medial_axis_pocket::set_debug(bool b) {debug=b;}
MICList medial_axis_pocket::get_mic_list() {return mic_list;}
    
/// find the largest MIC and add it to the output
void medial_axis_pocket::find_initial_mic() {
    MIC mic;
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
    mic.c2 = max_mic_pos;
    mic.r2 = max_mic_radius;
    mic.c1 = mic.c2;
    mic.r1 = mic.r1;
    current_radius = max_mic_radius;
    current_center = max_mic_pos;
    previous_branch_center = max_mic_pos;
    previous_branch_radius = max_mic_radius;
    mic.c_prev = max_mic_pos;
    mic.r_prev = max_mic_radius;
    // find the edge on which we start machining.
    double max_adj_radius(-1);
    BOOST_FOREACH( HEEdge e, g.out_edge_itr(max_mic_vertex) ) {
        if ( (g[ g.target(e) ].dist() > max_adj_radius) && g[e].valid && g[e].type != OUTEDGE ) {
            max_adj_radius = g[ g.target(e) ].dist();
            current_edge = e;
        }
    }
    // stash the other out-edges for visiting later
    BOOST_FOREACH( HEEdge e, g.out_edge_itr(max_mic_vertex) ) {
        if ( e != current_edge ) {
            unvisited.push( branch_point(current_center, current_radius, e ) );
        }
    }
    if (debug) { std::cout << " start edge is: "; g.print_edge(current_edge); }
    new_branch=false;
    mic.new_branch = new_branch;
    mic_list.push_back(mic);
}


// return true if next mic was found and added to list
// false means end-of-operation
bool medial_axis_pocket::find_next_mic() {
    if ( current_edge == HEEdge() ) {
        if (debug) std::cout << "find_next_mic() end of operation. Nothing to do.\n";
        return false;
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
    if (debug) std::cout << "find_next_mic() target width " << w_target << "\n";

    if ( w_target > max_width ) {
        // since moving to the target vertex would give too large cut-width
        // we search on the current edge for the next MIC
        if (debug) {  std::cout << " searching on the current edge "; g.print_edge(current_edge); }
        // find a point on the current edge
        double next_radius = find_next_radius();
        if (debug) {  std::cout << " next_radius = " << next_radius << "\n"; }
        output_next_mic(next_radius, new_branch);
        return true;
    } else {
        // moving to the target edge does not give a cut-width that is large enough
        // we need to find a new edge for the next MIC
        
        // mark edge DONE. this means we have machined all MICs on this edge.
        mark_done(current_edge);
        if (debug) std::cout << "Finding new edge !\n";
        
        bool end_branch_mic;
        boost::tie( current_edge, end_branch_mic) = find_next_edge(); // move to the next edge
        if ( current_edge == HEEdge() ) { // invalid edge marks end of operation
            if (debug)  std::cout << "find_next_mic() end of operation.\n";
            return false;
        }
        
        if ( end_branch_mic ) {
            if (debug)  std::cout << "find_next_mic() end-of-branch MIC.\n";
            // this is unreliable, so comment out for now
            //output_next_mic(current_radius, false);
            return true;
        }
        
        double next_radius = find_next_radius();
        if (new_branch) {
            new_branch=false;
            output_next_mic(next_radius, true);
            return true;
        } else {
            output_next_mic(next_radius, false);
            return true;
        }
    }
}

// pop an unvisited edge from the stack
// or end-of-operation if the stack is empty
HEEdge medial_axis_pocket::find_next_branch() {
    if (unvisited.empty() ) {
        if (debug) std::cout << "find_next_branch(): no un-machined branches. end operation.\n";
        return HEEdge();
    } else {
        branch_point out = unvisited.top();
        if (debug) { std::cout << "find_next_branch(): next branch is "; g.print_edge(out.next_edge); }
        unvisited.pop();
        previous_branch_center = current_center;
        previous_branch_radius = current_radius;            
        current_center = out.current_center;
        current_radius = out.current_radius;
        new_branch = true;
        return out.next_edge;
    }
}

EdgeVector medial_axis_pocket::find_out_edges() {
    HEVertex trg = g.target(current_edge);
    EdgeVector out_edges;
    BOOST_FOREACH(HEEdge e, g.out_edge_itr(trg) ) {
        if ( e != g[current_edge].twin && g[e].valid && g[e].type != NULLEDGE && g[e].type != OUTEDGE ) {
            out_edges.push_back(e);
        }
    }
    return out_edges;
}

// output the next edge
// return true if we need a final MIC at the end of a branch
std::pair<HEEdge,bool> medial_axis_pocket::find_next_edge() {
    EdgeVector out_edges = find_out_edges();
    if (out_edges.empty() ) {
        if (debug) std::cout << "find_next_edge(): no out_edges. end of branch.\n";
        
        if ( current_radius > g[ g.target(current_edge) ].dist() ) {
            current_radius = g[ g.target(current_edge) ].dist();
            current_center = g[ g.target(current_edge) ].position; //point(current_radius);
            return std::make_pair(current_edge,true); // this outputs a final MIC at the end of a branch
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
            mark_done( out_edges[0] ); 
            current_edge = out_edges[0]; // jump to the next edge
            return find_next_edge(); // and try to see if that edge is valid.
        }
    } else {
        if (debug) std::cout << "find_next_edge(): too many out-edges. ERROR.\n";
        exit(-1);
        return std::make_pair(HEEdge(),false);
    }
}
    
void medial_axis_pocket::mark_done(HEEdge e) {
    edge_data[ e ].done = true;
    edge_data[ g[e].twin ].done = true;
}
    
// does HEEdge e have the next MIC we want ?
bool medial_axis_pocket::has_next_radius(HEEdge e) {
    // check if the edge e is one where we continue
    double r1 = current_radius;
    Point  c1 = g[current_edge].point(r1); 
    double r2 = g[ g.target(e) ].dist();
    Point  c2 = g[e].point(r2);

    double w_target = ( c2-c1 ).norm() + r2 - r1;
    if (debug) {
        std::cout << "has_next_radius() "; g.print_edge(e);
        std::cout << " taget width = " << w_target << "\n";
    }
    if (w_target<=0)
        return false;
        
    if ( w_target > max_width )
        return true;
    else
        return false;
}
    
// on the current edge, move from current_radius towards target_radius
// and find a radius-value that satisfies the cut-width constraint.
double medial_axis_pocket::find_next_radius() {
    CutWidthError t(g, max_width, current_edge, current_center, current_radius);
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=500;
    boost::math::tools::eps_tolerance<double> tol(30);
    double target_radius = g[ g.target(current_edge) ].dist();

    double trg_err = t(target_radius);
    double cur_err = t(current_radius);
    if ( debug ||  ( !(trg_err*cur_err < 0) ) ) {
        std::cout << "find_next_radius():\n";
        std::cout << " current edge: "; g.print_edge(current_edge);// << current_radius << "\n";
        std::cout << "    edge type: " << g[current_edge].type_str() << "\n";
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
void medial_axis_pocket::output_next_mic(double next_radius, bool branch) {
    MIC mic;
    if (debug) std::cout << "output_next_mic() next_radius = " << next_radius << "\n";

    //mic.c2 = g[current_edge].point(next_radius);
    //mic.r2 = next_radius;
    
    //mic.c1 = current_center;
    //mic.r1 = current_radius;
    Point c1 = current_center;
    Point c2 = g[current_edge].point(next_radius);
    double r1 = current_radius;
    double r2 = next_radius;
    
    mic.c1 = c1;
    mic.r1 = r1;
    mic.c2 = c2;
    mic.r2 = r2;
    
    if (c1 != c2) {
        std::vector<Point> tangents = bitangent_points2(c1,r1,c2,r2);
        mic.t1 = tangents[0]; //tang1; //2
        mic.t2 = tangents[1]; //3
        mic.t3 = tangents[2]; //4
        mic.t4 = tangents[3]; //5
    }
    mic.new_branch = branch;
    mic.c_prev = previous_branch_center;
    mic.r_prev = previous_branch_radius;

    mic_list.push_back(mic);
    
    current_radius = next_radius;
    current_center = g[current_edge].point(next_radius);
    
}

std::vector<Point> medial_axis_pocket::bitangent_points(Point c1, double r1, Point c2, double r2) {
        // find the bi-tangents and return them too.
    // see voronoi_bisectors.py
    double detM = c1.x*c2.y - c2.x*c1.y;
    double m = ( c1.y-c2.y ) / detM;
    double p = ( c2.x-c1.x ) / detM;
    double n = ( c2.y*r1 - c1.y*r2 ) / detM;
    double q = ( c1.x*r2 - c2.x*r1 ) / detM;
    // this quadratic will give the C-values for the lines
    std::vector<double> roots = numeric::quadratic_roots( m*m+p*p, 2*(m*n+p*q),  n*n+q*q-1);
    if (debug) {
        std::cout << " c1 = " << c1 << "\n";
        std::cout << " r1 = " << r1 << "\n";
        std::cout << " c2 = " << c2 << "\n";
        std::cout << " r2 = " << r2 << "\n";
        std::cout << " cutwidth = " << (c1-c2).norm()+r2-r1 << "\n";
        std::cout << " " << roots.size() << " bi-tangent roots\n"; 
    }
    if ( roots.empty() ) {
        std::cout << "bitangent_points() ERROR: " << roots.size() << " bi-tangent roots\n"; 
        assert(0);
        exit(-1);
    }
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
    std::vector<Point> out;
    out.push_back(tang1);
    out.push_back(tang2);
    out.push_back(tang3);
    out.push_back(tang4);
    return out;
}

// try a more robust solution approach

std::vector<Point> medial_axis_pocket::bitangent_points2(Point c1, double r1, Point c2, double r2) {
    std::vector<Point> out;
    Point bd1,bd2;
    if ( r1 == r2 ) {
        Point c1c2 = c2-c1;
        c1c2.normalize();
        bd1 = c1c2.xy_perp();
        bd2 = -1*c1c2.xy_perp();
    } else {
        double c1c2_dist = (c1-c2).norm();
        double dr = fabs(r1-r2);
        double bitang_length = sqrt( c1c2_dist*c1c2_dist + dr*dr );
        double area = dr*bitang_length;
        double height = area / c1c2_dist;
        double bitang_c1c2 = sqrt( bitang_length*bitang_length - height*height );
        Point cdir =  c2-c1;
        if (r1>r2)
            cdir = c1-c2;
        cdir.normalize(); 
        // the bitangents
        Point bit1 = bitang_c1c2*cdir + height*cdir.xy_perp();
        Point bit2 = bitang_c1c2*cdir - height*cdir.xy_perp();
        
        if (r1>r2) {
            bd1 = bit1.xy_perp(); 
            bd2 = -1*bit2.xy_perp(); 
        } else {
            bd1 = -1*bit2.xy_perp(); 
            bd2 = bit1.xy_perp(); 
        }
        bd1.normalize();
        bd2.normalize();
    }
    out.push_back( c1 + r1*bd1 );
    out.push_back( c1 + r1*bd2 );
    out.push_back( c2 + r2*bd1 );
    out.push_back( c2 + r2*bd2 );
    return out;
}


} // end namespace

// end file 
