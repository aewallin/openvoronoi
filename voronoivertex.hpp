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


/// As we incrementally construct the diagram voronoi-vertices can have one of these four different states. 
/// The status is updated as follows:
/// OUT-vertices will not be deleted
/// IN-vertices will be deleted
/// UNDECIDED-vertices have not been examied yet
/// NEW-vertices are constructed on OUT-IN edges
enum VoronoiVertexStatus {OUT, IN, UNDECIDED, NEW };

/// This is the permanent type of a vertex in the diagram. 
/// OUTER vertices are special vertices added in init(), should have degree==4
/// VERTEXGEN are vertex generators, should have degree==0
/// NORMAL are normal voronoi-vertices, should have degree==6  (degree 3 graph with double-edges)
enum VoronoiVertexType {OUTER, NORMAL, VERTEXGEN};

typedef std::map<VoronoiVertexType, unsigned int> VertexDegreeMap;

/// Base-class for a voronoi-diagram site, or generator.
class Site {
public:
    Site() {}
    virtual ~Site() {}
    /// return closest point on site to given point p
    virtual Point apex_point(const Point& p) = 0;
    virtual const Point position() const = 0;
};

/// point, or vertex site.
class PointSite : public Site {
public:
    PointSite() {}
    PointSite( const Point& p): _p(p) {}
    ~PointSite() {}
    Point apex_point(const Point& p) {
        return _p;
    }
    void position( const Point& p ) {
        _p = p;
    }
    const Point position() const {
        return _p;
    }
private:
    Point _p;
};


/// properties of a vertex in the voronoi diagram
/// an object of this type is held in the BGL-graph for each vertex.
class VoronoiVertex {
public:
    VoronoiVertex();
    /// construct vertex at position p with type t
    VoronoiVertex( Point p, VoronoiVertexStatus st);
    /// vertex with given position, status, and type
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
    
    /// initialize clerance-disk
    void init_dist(const Point& p) {
        r = dist(p);
    }
    /// update clearance-disk
    void update_dist(const Point& p) {
        double d = dist(p);
        if (d<r)
            r=d;
    }
    /// return distance to a point from this vertex
    double dist(const Point& p) const { return (position-p).norm_sq(); }
    /// return clearance-disk radius
    double dist() const { return r; }
    /// in-circle predicate 
    double in_circle(const Point& p) const { return dist(p) - r; }
    /// if this vertex is a Site (PointSite, or end-point of a line-segment or arc)
    /// then we store a pointer to the site here.
    Site* site;
protected:
    /// global vertex count
    static int count;
    static VertexDegreeMap expected_degree;
    double r;
};


} // end ocl namespace
#endif
