/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#ifndef VODI_VERTEX_HPP
#define VODI_VERTEX_HPP

#include <map>
#include <cmath>

#include "common/point.hpp"
#include "site.hpp"

namespace ovd {

/// As we incrementally construct the diagram voronoi-vertices can have one of these four different states. 
enum VoronoiVertexStatus {
    OUT,          /*!< OUT-vertices will not be deleted */
    IN,           /*!< IN-vertices will be deleted */
    UNDECIDED,    /*!< UNDECIDED-vertices have not been examied yet */
    NEW           /*!< NEW-vertices are constructed on OUT-IN edges */
};

/// This is the permanent type of a vertex in the diagram. 
enum VoronoiVertexType {
    OUTER,      /*!< OUTER vertices are special vertices added in init(), should have degree==4 */
    NORMAL,     /*!< NORMAL are normal voronoi-vertices, should have degree==6  (degree 3 graph with double-edges) */
    POINTSITE,  /*!< POINTSITE are point sites, should have degree==0 */
    ENDPOINT,   /*!< ENDPOINT vertices are end-points of line-segments or arc-segments */
    SEPPOINT,   /*!< separator start-vertices on a null-face */
    APEX,       /*!< APEX vertices split quadratic edges at their apex(closest point to site) */
    SPLIT       /*!< split-vertices of degree==2 to avoid loops in the delete-tree */
};

/// A map of this type is used by VoronoiDiagramChecker to check that all vertices
/// have the expected (correct) degree (i.e. number of edges)
typedef std::map<VoronoiVertexType, unsigned int> VertexDegreeMap;
/// face-descriptor
typedef unsigned int HEFace;
                                     
/// \brief A vertex in the voronoi diagram
///
/// an object of this type is held in the BGL-graph for each vertex.
class VoronoiVertex {
public:
    VoronoiVertex();
    VoronoiVertex( Point p, VoronoiVertexStatus st);
    VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t);
    VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t, double init_radius);
    VoronoiVertex( Point pos, VoronoiVertexStatus st, VoronoiVertexType t, Point initDist);
    VoronoiVertex( Point pos, VoronoiVertexStatus st, VoronoiVertexType t, Point initDist, double k3);
    
    virtual ~VoronoiVertex() {}
    /// reset vertex status
    void reset();
    friend class VoronoiDiagramChecker;
    /// initialize clerance-disk
    void init_dist(const Point& p) { r = dist(p); }
    /// return distance to a point from this vertex
    double dist(const Point& p) const { return (position-p).norm(); }
    /// set clearance-disk to zero
    void zero_dist() {r=0;}
    /// return clearance disk-radius
    double dist() const { return r; }
    /// in-circle predicate 
    double in_circle(const Point& p) const {
        //if ( r==0 && dist(p) == 0 ) 
        //    return -1;
        //else
            return dist(p) - r; 
    }
    /// reset the index count
    static void reset_count() { count = 0; }
// DATA
    /// index of vertex
    int index;
    /// vertex status. when the incremental algorithm runs
    /// vertices are marked: undecided, in, out, or new
    VoronoiVertexStatus status;
    /// The type of the vertex
    VoronoiVertexType type;
    
    /// \todo what is this? remove?
    double max_error;
    
    /// flag for indicating wether vertex is in the vertexQueue
    bool in_queue;
    /// the position of the vertex
    Point position;
    /// set alfa
    void set_alfa(const Point& dir);
    /// the offset-direction {-1,+1} to the newly inserted site.
    double k3; 
    /// angle for a null-vertex
    double alfa;
    /// if this is a null-face, a handle to the null-face 
    HEFace null_face;
    /// the face of this vertex, if the vertex is a point-site
    HEFace face; 
protected:
    void init();
    void init(Point p, VoronoiVertexStatus st);
    void init(Point p, VoronoiVertexStatus st, VoronoiVertexType t);
    void init(Point p, VoronoiVertexStatus st, VoronoiVertexType t, Point initDist);
    void init(Point p, VoronoiVertexStatus st, VoronoiVertexType t, Point initDist, double k3);
    /// global vertex count
    static int count; // hold this in hedigraph instead?
    /// map for checking topology correctness
    static VertexDegreeMap expected_degree;
    /// clearance-disk radius, i.e. the closest site is at this distance
    double r;
};

} // end namespace
#endif
