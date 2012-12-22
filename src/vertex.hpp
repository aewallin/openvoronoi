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

/// As we incrementally construct the diagram voronoi-vertices 
/// can have one of these four different states. 
enum VertexStatus {
    OUT,          /*!< OUT-vertices will not be deleted */
    IN,           /*!< IN-vertices will be deleted */
    UNDECIDED,    /*!< UNDECIDED-vertices have not been examied yet */
    NEW           /*!< NEW-vertices are constructed on OUT-IN edges */
};

/// This is the permanent type of a vertex in the diagram. 
enum VertexType {
    OUTER,      /*!< OUTER vertices are special vertices added in init(), should have degree==4 */
    NORMAL,     /*!< NORMAL are normal voronoi-vertices, should have degree==6  (degree 3 graph with double-edges) */
    POINTSITE,  /*!< POINTSITE are point sites, should have degree==0 */
    ENDPOINT,   /*!< ENDPOINT vertices are end-points of line-segments or arc-segments */
    SEPPOINT,   /*!< separator start-vertices on a null-face */
    APEX,       /*!< APEX vertices split quadratic edges at their apex(closest point to site) */
    SPLIT       /*!< split-vertices of degree==2 to avoid loops in the delete-tree */
};




                                     
/// \brief A vertex in the voronoi diagram
///
/// an object of this type is held in the BGL-graph for each vertex.
class VoronoiVertex {
public:
    VoronoiVertex( Point p, VertexStatus st, VertexType t);
    VoronoiVertex( Point p, VertexStatus st, VertexType t, double init_radius);
    VoronoiVertex( Point pos, VertexStatus st, VertexType t, Point initDist);
    VoronoiVertex( Point pos, VertexStatus st, VertexType t, Point initDist, double k3);
    virtual ~VoronoiVertex();

    typedef unsigned int HEFace; ///< face-descriptor

    void reset_status();
    friend class VoronoiDiagramChecker;
    void init_dist(const Point& p);
    double dist(const Point& p) const;
    void zero_dist();
    double dist() const; 
    double in_circle(const Point& p) const; 
    static void reset_count(); 
    void set_alfa(const Point& dir); ///< set alfa. This is only for debug-drawing of null-face vertices.
// DATA
    
    int index; ///< unique integer index of vertex
    VertexStatus status; ///< vertex status. updated/changed during an incremental graph update
    VertexType type; ///< The type of the vertex. Never(?) changes
    double max_error; ///< \todo what is this? remove?
    bool in_queue; ///< flag for indicating wether vertex is in the vertexQueue
    Point position; ///< the position of the vertex.
    double k3;  ///< the offset-direction {-1,+1} of this vertex to the newly inserted site.
    double alfa; ///< diangle for a null-vertex. only for debug-drawing
    HEFace null_face; ///< if this is a null-face, a handle to the null-face 
    HEFace face; ///< the face of this vertex, if the vertex is a point-site
protected:
    void init();
    void init(Point p, VertexStatus st);
    void init(Point p, VertexStatus st, VertexType t);
    void init(Point p, VertexStatus st, VertexType t, Point initDist);
    void init(Point p, VertexStatus st, VertexType t, Point initDist, double k3);
    static int count; ///< global vertex count \todo hold this in hedigraph instead?
    /// A map of this type is used by VoronoiDiagramChecker to check that all vertices
    /// have the expected (correct) degree (i.e. number of edges)
    typedef std::map<VertexType, unsigned int> VertexDegreeMap;
    static VertexDegreeMap expected_degree; ///< map for checking topology correctness
    double r; ///< clearance-disk radius, i.e. the closest Site is at this distance
private:
    VoronoiVertex();
};

} // end namespace
#endif
