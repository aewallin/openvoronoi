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
#include <cmath>

#include <boost/graph/adjacency_list.hpp>

#include "common/point.hpp"
#include "site.hpp"

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
/// POINTSITE are point sites, should have degree==0
/// NORMAL are normal voronoi-vertices, should have degree==6  (degree 3 graph with double-edges)
/// ENDPOINT vertices are end-points of line-segments or arc-segments
/// APEX vertices split quadratic edges at their apex(closest point to site)
enum VoronoiVertexType {OUTER, NORMAL, POINTSITE, ENDPOINT, SEPPOINT, APEX, SPLIT};

/// a map of this type is used by topology-checker to check that all vertices
/// have the expected (correct) degree (i.e. number of edges)
typedef std::map<VoronoiVertexType, unsigned int> VertexDegreeMap;

typedef unsigned int HEFace;
                                     
/// A vertex in the voronoi diagram
/// an object of this type is held in the BGL-graph for each vertex.
class VoronoiVertex {
public:
    VoronoiVertex();
    /// construct vertex at position p with type t
    VoronoiVertex( Point p, VoronoiVertexStatus st);
    /// vertex with given position, status, and type
    VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t);
    VoronoiVertex( Point pos, VoronoiVertexStatus st, VoronoiVertexType t, Point initDist);

    virtual ~VoronoiVertex() {}
    void init();
    /// reset status
    void reset();
    friend class VoronoiDiagramChecker;
    /// initialize clerance-disk
    void init_dist(const Point& p) { r = dist(p); }
    /// return distance to a point from this vertex
    double dist(const Point& p) const { return (position-p).norm(); }
    /// return clearance-disk radius
    void zero_dist() {r=0;}
    double dist() const { return r; }
    /// in-circle predicate 
    double in_circle(const Point& p) const { return dist(p) - r; }
    static void reset_count() { count = 0; }
// DATA
    int index;
    /// vertex status. when the incremental algorithm runs
    /// vertices are marked: undecided, in, out, or new
    VoronoiVertexStatus status;
    VoronoiVertexType type;
    bool in_queue;
    /// the position of the vertex
    Point position;
    void set_alfa(const Point& dir);
    //Site* site; // pointSite, if this is a point-site (required??)
    double k3; // the offset-direction {-1,+1} to the newly inserted site.
    double alfa; // angle for a null-vertex
    HEFace null_face;
    HEFace face;
protected:
    /// global vertex count
    static int count;
    /// map for checking topology correctness
    static VertexDegreeMap expected_degree;
    /// clearance-disk radius, i.e. the closest site is at this distance
    double r;
    
};

} // end ocl namespace
#endif
