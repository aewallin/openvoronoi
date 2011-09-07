/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  
 *  This file is part of OpenVoronoi.
 *
 *  OpenCAMlib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  OpenCAMlib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with OpenCAMlib.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VODI_VERTEX_HPP
#define VODI_VERTEX_HPP

#include <map>

#include "point.hpp"


namespace ovd {

/// voronoi-vertices can have one of these four different states
/// as we incrementally construct the diagram the type is updated as follows:
/// OUT-vertices will not be deleted
/// IN-vertices will be deleted
/// UNDECIDED-vertices have not been examied yet
/// NEW-vertices are constructed on OUT-IN edges
enum VoronoiVertexStatus {OUT, IN, UNDECIDED, NEW };
// OUTER vertices are special vertices added in init(), should have degree==4
// VERTEXGEN are vertex generators, should have degree==0
// NORMAL are normal voronoi-vertices, should have degree==6  (degree 3 graph with double-edges)
enum VoronoiVertexType {OUTER, NORMAL, VERTEXGEN};

typedef std::map<VoronoiVertexType, unsigned int> VertexDegreeMap;



/// properties of a vertex in the voronoi diagram
class VoronoiVertex {
public:
    VoronoiVertex();
    /// construct vertex at position p with type t
    VoronoiVertex( Point p, VoronoiVertexStatus st);
    VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t);
    void init();
    void reset();
    void set_generators(const Point& pi, const Point& pj, const Point& pk);
    /// based on precalculated J2, J3, J4, calculate the H determinant (in-circle predicate) for input Point pl
    /// Eq.(20) from Sugihara&Iri 1994
    /// H<0 means point p is inside the clearance-circle
    /// H>0 means point is outside clearance circle
    double detH(const Point& p) const;
    /// index of vertex
    int index;
    /// vertex status. when the incremental algorithm runs
    /// vertices are marked: undecided, in, out, or new
    VoronoiVertexStatus status;
    VoronoiVertexType type;
    bool in_queue;
    /// the position of the vertex
    Point position;
    
    typedef unsigned int HEFace; 
    HEFace face; // the face associated with this vertex, if type==VERTEXGEN
    friend class VoronoiDiagramChecker;
protected:
    /// based on previously calculated J2, J3, and J4, set the position of the vertex
    /// Eq.(24) from Sugihara&Iri 1994
    void set_position();
    /// set the J values
    void set_J(const Point& p1, const Point& p2, const Point& p3);
    /// calculate J2
    /// Eq(21) from Sugihara&Iri 1994
    /// see also Eq(4.6.4), page 256, of Okabe et al book
    double detH_J2(const Point& pi, const Point& pj);
    /// calculate J3
    /// Eq(22) from Sugihara&Iri 1994
    /// see also Eq(4.6.5), page 257, of Okabe et al book
    double detH_J3(const Point& pi, const Point& pj);
    /// calculate J4
    /// Eq(23) from Sugihara&Iri 1994
    /// see also Eq(4.6.6), page 257, of Okabe et al book
    double detH_J4(const Point& pi, const Point& pj);
// DATA
    /// the reference point for J-calculations and detH
    Point pk;
    /// J2 determinant
    double J2;
    /// J3 determinant
    double J3;
    /// J4 determinant
    double J4;
    /// global vertex count
    static int count;
    static VertexDegreeMap expected_degree;
};


} // end ocl namespace
#endif
