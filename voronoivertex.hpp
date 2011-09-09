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



class Site {
public:
    Site() {}
    virtual ~Site() {}
    virtual Point apex_point(const Point& p) = 0;
};

class PointSite : public Site {
public:
    PointSite( const Point& p): _p(p) {}
    ~PointSite() {}
    Point apex_point(const Point& p) {
        return _p;
    }
    void position( const Point& p ) {
        _p = p;
    }
    Point position() const {
        return _p;
    }
private:
    Point _p;
};


/// properties of a vertex in the voronoi diagram
class VoronoiVertex {
public:
    VoronoiVertex();
    /// construct vertex at position p with type t
    VoronoiVertex( Point p, VoronoiVertexStatus st);
    VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t);
    void init();
    void reset();
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
    
    void init_dist(const Point& p) {
        r = dist(p);
    }
    void update_dist(const Point& p) {
        double d = dist(p);
        if (d<r)
            r=d;
    }
    
    double dist(const Point& p) const {
        return (position-p).norm_sq(); 
    }
    double dist() const {
        return r;
    }
    double in_circle(const Point& p) {
        return dist(p) - r;
    }
    
    Site* site;
protected:
    /// global vertex count
    static int count;
    static VertexDegreeMap expected_degree;
    double r;
};


} // end ocl namespace
#endif
