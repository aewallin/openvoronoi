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

namespace ovd {






/// experimental medial-axis pocketing
class medial_axis_pocket {
public:
    medial_axis_pocket(HEGraph& gi);
    void set_width(double w);
    void run();
    //void run2();
    void set_debug(bool b);
    /// \brief Maximal Inscribed Circle. A combination of a Point and a clearance-disk radius.
    ///
    /// it is the responsibility of a downstream algorithm to lay down a toolpath
    /// that machines the area between c2 and c1
    /// c1 is the circle that is assumed already cut
    /// c2 is the new circle
    struct MIC {
        Point c1;       ///< center
        Point c2;       ///< center
        double r1;      ///< radius
        double r2;      ///< radius
        Point t1;       ///< bi-tangent point
        Point t2;       ///< bi-tangent point
        Point t3;       ///< bi-tangent point
        Point t4;       ///< bi-tangent point
        bool new_branch;   ///< is this a new branch?
        Point c_prev;      ///< for a new branch, the previous center
        double r_prev;     ///< for a new branch, the previous radius
    };
    typedef std::vector<MIC> MICList; ///< the list of MIC.s from one connected component of the medial-axis
    //MICList get_mic_list();
    std::vector<MICList> get_mic_components(); // {return ma_components;}
    std::pair<Point,double> edge_point(HEEdge e, double u); // used by the error-functor also. move somewhere else?
    
    
protected:
    /// \brief error-functor for find_next_u()
    /// \sa medial_axis_pocket
    class CutWidthError  {
    public:
        CutWidthError(medial_axis_pocket* ma, HEEdge ed, double wmax, Point cen1, double rad1);
        double operator()(const double x);
    private:
        medial_axis_pocket* m; ///< calling class
        HEEdge e;     ///< current edge
        double w_max; ///< desired cut-width
        Point c1;     ///< previous MIC center
        double r1;    ///< previous MIC radius
    };

    /// keep track of bool done true/false flag for each edge
    struct edata {
        edata();// { done = false; }
        bool done; ///< is edge done?
    };

    /// branch-data when we backtrack to machine an un-machined branch
    struct branch_point {
        branch_point(Point p, double r, HEEdge e);
        Point current_center;  ///< current center
        double current_radius; ///< current radius
        HEEdge next_edge;      ///< edge on which to start machining when we switch to the new branch
    };

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
    bool debug; ///< debug output flag
    std::vector<HEEdge> ma_edges; ///< the edges of the medial-axis
    std::map<HEEdge, edata> edge_data; ///< map from edge-descriptor to bool-flag
    HEGraph& g; ///< VD graph
    std::stack<branch_point> unvisited; ///< stack of unvisited branch_point:s
    HEEdge current_edge; ///< the current edge
    double current_radius; ///< current clearance-disk radius
    double current_u; ///< current position along edge. u is in [0,1]
    Point current_center; ///< current position
    bool new_branch;     ///< flag for indicating new branch
    Point previous_branch_center; ///< prev branch position
    double previous_branch_radius; ///< prev branch radius
    double max_width; ///< the max cutting-width
    MICList mic_list; ///< the result of the operation is a list of MICs 
    std::vector<MICList> ma_components; ///< algorithm output
};

} // end ovd namespace
// end file 
