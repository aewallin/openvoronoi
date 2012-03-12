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

#include <boost/math/tools/roots.hpp>

#include "graph.hpp"
#include "common/numeric.hpp"
#include "site.hpp"

namespace ovd
{

/// keep track of bool done true/false flag for each edge
struct edata {
    edata() {
        done = false;
    }
    bool done;
};
typedef std::pair<HEEdge , edata> Edata;


/// branch-data when we backtrack to machine an un-machined branch
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

/// \brief Maximal Inscribed Circle. A combination of a Point and a clearance-disk radius.
///
/// it is the responsibility of a downstream algorithm to lay down a toolpath
/// that machines the area between c2 and c1
/// c1 is the circle that is assumed already cut
/// c2 is the new circle
struct MIC {
    Point c1,c2;  // center
    double r1,r2; // radius
    Point t1,t2,t3,t4; // bi-tangent points
    bool new_branch;
    Point c_prev;
    double r_prev;
};

// the list of MICs from one connected component of the MA
typedef std::vector<MIC> MICList;

/// experimental medial-axis pocketing
class medial_axis_pocket {
public:
    medial_axis_pocket(HEGraph& gi);
    void set_width(double w);
    void run();
    void run2();
    void set_debug(bool b);
    MICList get_mic_list();
    std::vector<MICList> get_mic_components(); // {return ma_components;}
    std::pair<Point,double> edge_point(HEEdge e, double u); // used by the error-functor also. move somewhere else?
protected:
    bool find_initial_mic();
    bool find_next_mic();
    HEEdge find_next_branch();
    EdgeVector find_out_edges();
    std::pair<HEEdge,bool> find_next_edge();
    void mark_done(HEEdge e);
    bool has_next_radius(HEEdge e); 
    std::pair<double,double> find_next_u();
    void output_next_mic(double next_u, double next_radius, bool branch);
    std::vector<Point> bitangent_points(Point c1, double r1, Point c2, double r2);
    double cut_width(Point c1, double r1, Point c2, double r2);
//DATA
    bool debug;
    std::vector<HEEdge> ma_edges; // the edges of the medial-axis
    std::map<HEEdge, edata> edge_data;
    HEGraph& g; // VD graph
    std::stack<branch_point> unvisited;

    HEEdge current_edge;
    double current_radius;
    double current_u;
    Point current_center;
    
    // flag for indicating new branch
    bool new_branch;
    Point previous_branch_center;
    double previous_branch_radius;
    
    // the max cutting-width
    double max_width;
    // the result of the operation is a list of MICs 
    MICList mic_list;
    std::vector<MICList> ma_components;
    //int max_mic_count;
};

/// \brief error-functor for medial_axis_pocket::find_next_u()
class CutWidthError  {
public:
    CutWidthError(medial_axis_pocket* ma, HEEdge ed, double wmax, Point cen1, double rad1) 
    : m(ma), e(ed), w_max(wmax),  c1(cen1), r1(rad1) {}
    double operator()(const double x) {
        // w_max = | c2 - c1 | + r2 - r1
        Point c2; // = m->edge_point(x); //g[e].point(x); // current MIC center
        double r2; // = x; // current MIC radius
        boost::tie(c2,r2) = m->edge_point(e,x);
        double w = (c2-c1).norm() + r2 - r1; // this is the cut-width
        return w-w_max; // error compared to desired cut-width
    }
private:
    medial_axis_pocket* m;
    HEEdge e;
    double w_max; // desired cut-width
    Point c1; // previous MIC center
    double r1; // previous MIC radius
};

} // end namespace

// end file 
