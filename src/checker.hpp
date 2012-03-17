/*  
 *  Copyright 2010-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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
#include "graph.hpp"

namespace ovd {

class VoronoiDiagram;

/// this class provides sanity-checks for the VoronoiDiagram class
class VoronoiDiagramChecker {
public:
    /// \param gi input graph
    VoronoiDiagramChecker(HEGraph& gi);// : g(gi) {}
    ~VoronoiDiagramChecker();// {}

    bool is_valid();
    bool face_count_equals_generator_count();
    bool vertex_degree_ok();
    bool all_in( const VertexVector& q);
    bool noUndecidedInFace( HEFace f );
    bool faceVerticesConnected( HEFace f, VertexStatus Vtype );
    bool current_face_equals_next_face( HEEdge e); 
    bool face_ok(HEFace f, bool debug=false);
    bool all_faces_ok();
    bool check_edge(HEEdge e) const;
private:
    HEGraph& g; ///< vd-graph
};

} // end ovd namespace

// end voronoidiagram_checker.hpp
