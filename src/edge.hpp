/* 
 *  Copyright 2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include <cassert>
#include <cmath>
#include <boost/array.hpp>

#include <boost/graph/adjacency_list.hpp>

#include "common/point.hpp"
#include "site.hpp"
#include "solvers/solution.hpp"

namespace ovd {

#define OUT_EDGE_CONTAINER boost::listS 
//#define OUT_EDGE_CONTAINER boost::vecS 


// note: cannot use vecS since remove_vertex invalidates iterators/edge_descriptors (?)
#define VERTEX_CONTAINER boost::listS
#define EDGE_LIST_CONTAINER boost::listS

// type of edge-descriptors in the graph
typedef boost::adjacency_list_traits<OUT_EDGE_CONTAINER, 
                                     VERTEX_CONTAINER, 
                                     boost::bidirectionalS, 
                                     EDGE_LIST_CONTAINER >::edge_descriptor HEEdge;
                                     
typedef unsigned int HEFace;

enum VoronoiEdgeType {LINE, LINELINE, PARA_LINELINE, OUTEDGE, PARABOLA, ELLIPSE, HYPERBOLA, SEPARATOR, NULLEDGE, LINESITE};


/// properties of an edge in the voronoi diagram
/// each edge stores a pointer to the next HEEdge 
/// and the HEFace to which this HEEdge belongs
class EdgeProps {
public:
    EdgeProps();
    EdgeProps(HEEdge n, HEFace f): next(n), face(f), has_null_face(false), valid(true) {}
    /// create edge with given next, twin, and face
    EdgeProps(HEEdge n, HEEdge t, HEFace f): next(n), twin(t), face(f), has_null_face(false), valid(true) {}
    /// the next edge, counterclockwise, from this edge
    HEEdge next; 
    /// the twin edge
    HEEdge twin;
    /// the face to which this edge belongs
    HEFace face; // each face corresponds to an input Site/generator
    HEFace null_face;
    bool has_null_face;
    
    double k; // offset-direction from the adjacent site, either +1 or -1
    VoronoiEdgeType type;
    
    boost::array<double,8> x;
    boost::array<double,8> y;
    bool sign; // choose either +/- in front of sqrt()
    
    Point point(double t) const; 
    double error(Solution& sl) const;
    double minimum_t( Site* s1, Site* s2);
    void copy_parameters(EdgeProps& other);   
    void set_parameters(Site* s1, Site* s2, bool sig);
    void set_sep_parameters(Point& endp, Point& p);
    EdgeProps &operator=(const EdgeProps &p);
    bool valid; // for filtering graph
    bool inserted_direction; // true if linesite-edge inserted in this direction
private:
    Point projection_point(Solution& sl) const;
    double minimum_pp_t(Site* s1, Site* s2);
    double minimum_pl_t(Site* s1, Site* s2);

    void set_pp_parameters(Site* s1, Site* s2);
    void set_pl_parameters(Site* s1, Site* s2);
    void set_ll_parameters(Site* s1, Site* s2);
    void set_ll_para_parameters(Site* s1, Site* s2);
    void print_params() const;
    

};

} // end namespace
