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

#include "graph.hpp"
#include "site.hpp"

namespace ovd
{

/// \brief Line- or arc-vertex of an offset curve.
///
/// \todo this duplicates the idea of the Ofs class. Remove this or Ofs!
struct OffsetVertex {
    Point p;  ///< position (start)
    double r; ///< arc radius (line-vertex is indicated by radius of -1)
    Point c;  ///< arc center
    bool cw;  ///< clockwise (or not)
    HEFace f; ///< corresponding face in the vd-graph
    /// ctor
    OffsetVertex(Point pi, double ri, Point ci, bool cwi, HEFace fi): p(pi), r(ri), c(ci), cw(cwi), f(fi) {}
    /// ctor
    OffsetVertex(Point pi): p(pi), r(-1.), cw(false), f(0) {}
};

/// a single offset loop
struct OffsetLoop {
    std::list<OffsetVertex> vertices;   ///< list of offsetvertices in this loop
    double offset_distance;             ///< offset distance for this loop
    /// add an offsetvertex to this loop
    void push_back(OffsetVertex v) { vertices.push_back(v); }
};

/// multiple loops. the output of the algorithm
typedef std::vector<OffsetLoop> OffsetLoops;

/// \brief From a voronoi-diagram, generate offsets.
///
/// an offset is allways a closed loop.
/// the loop consists of offset-elements from each face that the loop visits.
/// each face is associated with a Site, and the offset element from
/// - a point-site is a circular arc
/// - a line-site is a line
/// - an arc is a circular arc
///
/// This class produces offsets at the given offset-distance on the entire
/// voronoi-diagram. To produce offsets only inside or outside a given geometry,
/// use a filter first. The filter sets the valid-property of edges, so that offsets
/// are not produced on faces with one or more invalid edge.
class Offset {
public:
    /// \param gi vd-graph
    Offset(HEGraph& gi): g(gi) {
        face_done.clear();
        face_done.assign( g.num_faces(), 1 );
    }
    /// print stats
    void print();
    /// create offsets at offset distance \a t
    OffsetLoops offset(double t);
protected:
    bool find_start_face(HEFace& start);
    void offset_loop_walk(HEFace start, double t);
    OffsetVertex offset_element_from_face(HEFace current_face, HEEdge current_edge, HEEdge next_edge, double t);
    bool edge_mode(HEEdge e, double t);
    bool find_cw(Point start, Point center, Point end);
    HEEdge find_next_offset_edge(HEEdge e, double t, bool mode);
    void set_flags(double t);
    bool t_bracket(double a, double b, double t);
    void print_status();
    
    OffsetLoops offset_list; ///< list of output offsets
private:
    Offset(); // don't use.
    HEGraph& g; ///< vd-graph
    /// hold a 0/1 flag for each face, indicating if an offset for this face has been produced or not.
    std::vector<unsigned char> face_done;
};


} // end ovd namespace
// end file offset.hpp
