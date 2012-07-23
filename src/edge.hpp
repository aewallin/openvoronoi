/* 
 *  Copyright 2011-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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
#define VERTEX_CONTAINER boost::listS
#define EDGE_LIST_CONTAINER boost::listS

/// edge-descriptors
typedef boost::adjacency_list_traits<OUT_EDGE_CONTAINER, 
                                     VERTEX_CONTAINER, 
                                     boost::bidirectionalS, 
                                     EDGE_LIST_CONTAINER >::edge_descriptor HEEdge;
/// face descriptor                                     
typedef unsigned int HEFace;

/// edge type 
enum EdgeType {
    LINE,          /*!< Line edge between PointSite and PointSite */ 
    LINELINE,      /*!< Line edge between LineSite and LineSite */ 
    PARA_LINELINE, /*!< Line edge between LineSite and LineSite (parallel case) */ 
    OUTEDGE,       /*!< special outer edge set by initialize() */ 
    PARABOLA,      /*!< Parabolic edge between PointSite and LineSite */ 
    ELLIPSE, 
    HYPERBOLA, 
    SEPARATOR,     /*!< Separator edge between PointSite (endpoint) and LineSite or ArcSite */ 
    NULLEDGE,      /*!< zero-length null-edge around a PointSite which is and endpoint */ 
    LINESITE,       /*!< pseudo-edge corresponding to a LineSite */ 
    ARCSITE       /*!< pseudo-edge corresponding to a LineSite */ 
    };
/*
* bisector formulas
* x = x1 - x2 - x3*t +/- x4 * sqrt( square(x5+x6*t) - square(x7+x8*t) )
* (same formula for y-coordinate)
* line (line/line)
* parabola (circle/line)
* hyperbola (circle/circle)
* ellipse (circle/circle)
*/


/// \brief properties of an edge in the VoronoiDiagram
///
/// each edge stores a pointer to the next HEEdge 
/// and the HEFace to which this HEEdge belongs
class EdgeProps {
public:
    EdgeProps();
    /// create edge with given next and face
    EdgeProps(HEEdge n, HEFace f): next(n), face(f), has_null_face(false), valid(true) {}
    /// create edge with given next, twin, and face
    EdgeProps(HEEdge n, HEEdge t, HEFace f): next(n), twin(t), face(f), has_null_face(false), valid(true) {}
    
    HEEdge next;     ///< the next edge, counterclockwise on the face, from this edge
    HEEdge twin; ///< the twin edge
    HEFace face; ///< the face to which this edge belongs
    HEFace null_face; ///< face descriptor of null-face (parallel SEPARATOR case ?)
    bool has_null_face; ///< flag for indicating null-face (this means null_face is set?)
    
    double k; ///< offset-direction from the adjacent site, either +1 or -1
    EdgeType type; ///< the type of this edge
    
    boost::array<double,8> x; ///< 8-parameter parametrization
    boost::array<double,8> y; ///< 8-parameter parametrization
    bool sign; ///< flag to choose either +/- in front of sqrt()

    Point point(double t) const; 
    double minimum_t( Site* s1, Site* s2);
       
    void set_parameters(Site* s1, Site* s2, bool sig);
    void set_sep_parameters(Point& endp, Point& p);
    EdgeProps &operator=(const EdgeProps &p);
    bool valid; ///< flag set by Filter, for use by downstream algorithms
    bool inserted_direction; ///< true if ::LINESITE-edge inserted in this direction
    std::string type_str() const;
private:
    double minimum_pp_t(Site* s1, Site* s2);
    double minimum_pl_t(Site* s1, Site* s2);
    double minimum_pa_t(Site* s1, Site* s2);

    void set_pp_parameters(Site* s1, Site* s2);
    void set_pl_parameters(Site* s1, Site* s2);
    void set_ll_parameters(Site* s1, Site* s2);
    void set_ll_para_parameters(Site* s1, Site* s2);
    void set_pa_parameters(Site* s1, Site* s2);
    void set_la_parameters(Site* s1, Site* s2);
    void print_params() const;
};

} // end ovd namespace
