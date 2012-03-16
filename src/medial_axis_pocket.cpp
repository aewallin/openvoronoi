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
            edge_data.insert( std::make_pair(e, ed ) );
        }
    }
    current_edge = HEEdge();
    max_width = 0.05;
    debug = false;
    //max_mic_count=300; // limit output size in debug mode
}

/// set the maximum cut width
void medial_axis_pocket::set_width(double w) {max_width=w;}
/// set debug mode
void medial_axis_pocket::set_debug(bool b) {debug=b;}
/// return output
//medial_axis_pocket::MICList medial_axis_pocket::get_mic_list() {return mic_list;}
/// return the algorithm output
std::vector<medial_axis_pocket::MICList> medial_axis_pocket::get_mic_components() {return ma_components;}

/*
/// run the algorithm (single connected component)
void medial_axis_pocket::run() {
    find_initial_mic();
    while (find_next_mic()) {}
    if (debug) std::cout << "medial_axis_pocket::run() done. generated " << mic_list.size() << " MICs \n";
}*/

/// many component run
void medial_axis_pocket::run() {
    mic_list.clear();
    while ( find_initial_mic() ) {
        while (find_next_mic()) {}
        if (debug) std::cout << "medial_axis_pocket::run() component done. generated " << mic_list.size() << " MICs \n";
        ma_components.push_back(mic_list);
        mic_list.clear();
    }
    //if (debug) std::cout << "medial_axis_pocket::run() component done. generated " << mic_list.size() << " MICs \n";
    
    if (debug) std::cout << "medial_axis_pocket::run() all done. generated " << ma_components.size() << " components \n";

}

    
/// find the largest MIC and add it to the output
bool medial_axis_pocket::find_initial_mic() {
    MIC mic;
    // find the vertex with the maximum radius mic
    double max_mic_radius(-1);
    Point max_mic_pos(0,0);
    HEVertex max_mic_vertex = HEVertex();
    bool found(false);
    BOOST_FOREACH( HEEdge e, ma_edges ) {
        HEVertex src = g.source(e);
        if ( !edge_data[e].done && g[src].dist() > max_mic_radius ) {
            max_mic_radius = g[src].dist();
            max_mic_pos = g[src].position; 
            max_mic_vertex = src;
            found = true;
        }
    }
    if (!found)
        return false;
        
    if (debug) { std::cout << "find_initial_mic() max mic is c="<< max_mic_pos << " r=" << max_mic_radius << "\n"; }

    mic.c2 = max_mic_pos;
    mic.r2 = max_mic_radius;
    mic.c1 = mic.c2;
    mic.r1 = mic.r1;
    current_u = 0;
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
    if (debug) { std::cout << "find_initial_mic() start edge is: "; g.print_edge(current_edge); }
    new_branch=false;
    mic.new_branch = new_branch;
    mic_list.push_back(mic);
    return true;
}

/// \param p position of degree-3 branch
/// \param r clearance-disk radius
/// \param e edge on which to start machining
medial_axis_pocket::branch_point::branch_point(Point p, double r, HEEdge e) {
    current_center = p;
    current_radius = r;
    next_edge = e;
}
        
/// return true if next mic was found and added to list.
/// false means end-of-operation
bool medial_axis_pocket::find_next_mic() {
    if ( current_edge == HEEdge() ) {
        if (debug) std::cout << "find_next_mic() end of operation. Nothing to do.\n";
        return false;
    }
    //if ( debug && mic_list.size() > max_mic_count ) {
    //    std::cout << " max_mic_count reached. stopping.\n";
    //    return false;
    //}
    // find a point on current-edge so that we get the desired 
    // cut-width
    //  w_max = | c2 - c1 | + r2 - r1
    
    Point c2;
    double r2;
    boost::tie(c2,r2) = edge_point(current_edge, 1.0 );
    double w_target = cut_width(current_center, current_radius, c2, r2); //( c2-c1 ).norm() + r2 - r1;

    if ( w_target > max_width ) {
        // since moving to the target vertex would give too large cut-width
        // we search on the current edge for the next MIC
        if (debug) {  std::cout << " searching on the current edge "; g.print_edge(current_edge); }
        // find a point on the current edge
        double next_u;
        double next_radius; // = find_next_radius();
        boost::tie(next_u,next_radius) = find_next_u();
        if (debug) {  std::cout << " next_radius = " << next_radius << "\n"; }
        output_next_mic(next_u, next_radius, new_branch);
        return true;
    } else {
        // moving to the target edge does not give a cut-width that is large enough
        // we need to find a new edge for the next MIC
        
        // mark edge DONE. this means we have machined all MICs on this edge.
        mark_done(current_edge);
        if (debug) std::cout << "find_next_mic() Finding new edge !\n";
        
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
        
        double next_u;
        double next_radius;
        boost::tie(next_u,next_radius) = find_next_u();
        if (new_branch) {
            new_branch=false;
            output_next_mic(next_u, next_radius, true);
            return true;
        } else {
            output_next_mic(next_u, next_radius, false);
            return true;
        }
    }
}

/// pop an unvisited edge from the stack
/// or end-of-operation if the stack is empty
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
        if ( !edge_data[ out.next_edge ].done )
            return out.next_edge;
        else
            return find_next_branch();
    }
}
/// find out-edges of current_edge
EdgeVector medial_axis_pocket::find_out_edges() {
    HEVertex trg = g.target(current_edge);
    EdgeVector out_edges;
    BOOST_FOREACH(HEEdge e, g.out_edge_itr(trg) ) {
        if ( e != g[current_edge].twin && 
             g[e].valid && !edge_data[e].done && 
             g[e].type != NULLEDGE && 
             g[e].type != OUTEDGE ) {
            out_edges.push_back(e);
        }
    }
    return out_edges;
}

/// find the next edge
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
            current_u = 0;
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
            current_u = 0;
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
            current_u=0;
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

/// mark edge and its twin done
void medial_axis_pocket::mark_done(HEEdge e) {
    edge_data[ e ].done = true;
    edge_data[ g[e].twin ].done = true;
}

/// does HEEdge e have the next MIC we want ?
bool medial_axis_pocket::has_next_radius(HEEdge e) {
    // check if the edge e is one where we continue
    double r2;
    Point  c2;
    boost::tie(c2,r2) = edge_point(e,1.0);

    double w_target = cut_width(current_center, current_radius, c2, r2);
    if (debug) {
        CutWidthError t(this, e,max_width, current_center, current_radius);
        std::cout << "has_next_radius() ?"<< ( w_target > max_width ) <<" "; g.print_edge(e);
        std::cout << "has_next_radius() err src "<< t(0) <<"\n";
        std::cout << "has_next_radius() err trg "<< t(1) <<"\n";
    }
    if (w_target<=0)
        return false;
        
    if ( w_target > max_width )
        return true;
    else
        return false;
}

/// find the next u-value that produces the desired cut-width
std::pair<double,double> medial_axis_pocket::find_next_u() {
    CutWidthError t(this, current_edge, max_width, current_center, current_radius);
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=500;
    boost::math::tools::eps_tolerance<double> tol(30);
    double trg_err = t(1.0);
    double cur_err = t(current_u);
    if ( debug ||  ( !(trg_err*cur_err < 0) ) ) {
        std::cout << "find_next_u():\n";
        std::cout << " current edge: "; g.print_edge(current_edge);// << current_radius << "\n";
        std::cout << "    edge type: " << g[current_edge].type_str() << "\n";
        std::cout << " source c= "<< g[ g.source(current_edge) ].position << " r= " << g[ g.source(current_edge) ].dist() << " err= " <<  t(0.0) <<"\n";
        std::cout << " target c= "<< g[ g.target(current_edge) ].position << " r= " << g[ g.target(current_edge) ].dist() << " err= " <<  t(1.0) <<"\n";
        std::cout << " current c1=" << current_center << " r1=" << current_radius << "\n";
        std::cout << " current u = " << current_u << "\n";
        std::cout << " error at current = " << t(current_u) << "\n";
        std::cout << " error at target = " << t(1.0) << "\n";
    }
    //double min_r = std::min(current_radius, target_radius);
    //double max_r = std::max(current_radius, target_radius);
    Result r1 = boost::math::tools::toms748_solve(t, current_u, 1.0, tol, max_iter);
    double rnext;
    Point pnext;
    boost::tie( pnext, rnext ) = edge_point( current_edge, r1.first );
    return std::make_pair( r1.first, rnext );
}
    
/// \brief output the next MIC
///
/// based on the output here a downstream algorithm
/// will lay out the pattern: lead-out, rapid, lead-in, bi-tangent, cut-arc, bi-tangent
void medial_axis_pocket::output_next_mic(double next_u, double next_radius, bool branch) {
    MIC mic;
    Point c1 = current_center;
    double r1 = current_radius;
    Point c2;
    double r2;
    boost::tie(c2,r2) = edge_point(current_edge, next_u); //g[current_edge].point(next_radius);
    //double r2 = next_radius;
    if (debug) {
        std::cout << "output_next_mic(): \n";
        std::cout << " next_u = " << next_u << "\n";
        std::cout << " next_radius = " << next_radius << "\n";
        std::cout << " c= " << c2 << " r= " << r2 << "\n";
    }
    
    mic.c1 = c1;
    mic.r1 = r1;
    mic.c2 = c2;
    mic.r2 = r2;

    if (c1 != c2) {
        std::vector<Point> tangents = bitangent_points(c1,r1,c2,r2);
        mic.t1 = tangents[0]; //tang1; //2
        mic.t2 = tangents[1]; //3
        mic.t3 = tangents[2]; //4
        mic.t4 = tangents[3]; //5
    }
    mic.new_branch = branch;
    mic.c_prev = previous_branch_center;
    mic.r_prev = previous_branch_radius;
    mic_list.push_back(mic);
    current_radius = r2; //next_radius;
    current_center = c2; //g[current_edge].point(next_radius);
    current_u = next_u;
}

/// bi-tangent points to c1-c2
std::vector<Point> medial_axis_pocket::bitangent_points(Point c1, double r1, Point c2, double r2) {
    std::vector<Point> out;
    Point bd1,bd2;
    if ( r1 == r2 ) {
        Point c1c2 = c2-c1;
        c1c2.normalize();
        bd1 = -1*c1c2.xy_perp();
        bd2 = c1c2.xy_perp();
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


/// \return cut-width from c1/r1 to c2/r2
///
/// - (c1,r1) is the previously machined MIC
/// - (c2,r2) is the new MIC
/// the maximum cut-width when cutting c2 is
///  w_max = | c2 - c1 | + r2 - r1
double medial_axis_pocket::cut_width(Point c1, double r1, Point c2, double r2) {
    return ( c2-c1 ).norm() + r2 - r1;
}

/// find a point on current_edge at u [0,1]
std::pair<Point,double> medial_axis_pocket::edge_point(HEEdge e, double u) {
    Point p;
    double r;
    // on line-type edges return
    // p = p1 + u ( p2 - p1 )
    /*
    enum VoronoiEdgeType {
    LINE,          // PontSite PointSite
    LINELINE,      // LineSite LineSite
    PARA_LINELINE, // LineSite LineSite
    OUTEDGE, 
    PARABOLA,      // LineSite PointSite
    ELLIPSE, 
    HYPERBOLA, 
    SEPARATOR,     // LineSite PointSite(endpoint) 
    NULLEDGE, 
    LINESITE
    };
    */
    if ( g[e].type == LINE ||
        g[e].type == LINELINE ||
        g[e].type == PARA_LINELINE ) // note SEPARATOR shouldn't occur in the medial-axis
        {
        Point src = g[ g.source(e) ].position;
        Point trg = g[ g.target(e) ].position;
        p = src + u * (trg-src);
        HEFace f = g[e].face;
        Site* s = g[f].site;
        Point pa = s->apex_point(p);
        r = (p-pa).norm();
    } else if ( g[e].type == PARABOLA ) {
        // use the existing t-parametrisation (?)
        HEVertex src_v = g.source(e);
        HEVertex trg_v = g.target(e);
        double src_t = g[src_v].dist();
        double trg_t = g[trg_v].dist();
        double u_t = src_t + u*(trg_t-src_t);
        p = g[e].point(u_t);
        r=u_t;
    } else {
        std::cout << "ERROR medial_axis_pocket::edge_point() unsuppoerted edge type.\n";
        std::cout << " type= " << g[e].type_str() << "\n";
        exit(-1);
    }
    
    return std::make_pair(p,r);
}

/// \param ma calling medial_axis_pocket object
/// \param ed edge on which to position MIC
/// \param wmax maximum cut width
/// \param cen1 previous MIC center
/// \param rad1 previous MIC radius
medial_axis_pocket::CutWidthError::CutWidthError(medial_axis_pocket* ma, 
                            HEEdge ed, double wmax, Point cen1, double rad1) 
: m(ma), e(ed), w_max(wmax),  c1(cen1), r1(rad1) {}

/// cut-width if next MIC positioned at \a x
double medial_axis_pocket::CutWidthError::operator()(const double x) {
    // w_max = | c2 - c1 | + r2 - r1
    Point c2; // = m->edge_point(x); //g[e].point(x); // current MIC center
    double r2; // = x; // current MIC radius
    boost::tie(c2,r2) = m->edge_point(e,x);
    double w = (c2-c1).norm() + r2 - r1; // this is the cut-width
    return w-w_max; // error compared to desired cut-width
}

medial_axis_pocket::edata::edata() { done = false; }

} // end namespace

// end file 
