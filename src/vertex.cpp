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

#include <cassert>
#include <limits>

#include <boost/assign.hpp>

#include "vertex.hpp"
#include "common/numeric.hpp"

namespace ovd {

int VoronoiVertex::count = 0;

// the expected degree of a vertex. checked by topology-checker
VoronoiVertex::VertexDegreeMap VoronoiVertex::expected_degree = boost::assign::map_list_of 
    (OUTER,4)     // special outer vertices
    (NORMAL,6)    // normal vertex in the graph
    (POINTSITE,0) // point site
    (ENDPOINT,6)  // end-point of line or arc
    (SEPPOINT,6)  // end-point of separator
    (SPLIT,4)     // split point, to avoid loops in delete-tree
    (APEX,4) ;    // apex point on quadratic bisector
    
/// ctor with given status and type
VoronoiVertex::VoronoiVertex( Point p, VertexStatus st, VertexType t) {
    init(p,st,t);
}
/// ctor with initial apex Point
VoronoiVertex::VoronoiVertex( Point p, VertexStatus st, VertexType t, Point initDist) {
    init(p,st,t,initDist);
}
/// ctor with initial k3-value
VoronoiVertex::VoronoiVertex( Point p, VertexStatus st, VertexType t, Point initDist, double lk3) {   
    init(p,st,t,initDist,lk3);
}
/// ctor with initial clearance-disk radius
VoronoiVertex::VoronoiVertex( Point p, VertexStatus st, VertexType t, double init_radius) {
    init(p,st,t);
    r = init_radius;
}
VoronoiVertex::~VoronoiVertex() {}

/// set index, increase count, initialize in_queue to false.
void VoronoiVertex::init() {
    index = count;
    count++;
    in_queue = false;
    alfa=-1; // invalid/non-initialized alfa value
    null_face = std::numeric_limits<HEFace>::quiet_NaN();    
    type = NORMAL;
    face = std::numeric_limits<HEFace>::quiet_NaN();  
    max_error = 0;
}

/// set position and status
void VoronoiVertex::init(Point p, VertexStatus st) {
    init();
    position=p;
    status=st;
}
/// set position, status and type
void VoronoiVertex::init(Point p, VertexStatus st, VertexType t) {
    init(p,st);
    type = t;
}
/// set position, status, type, and clearance-disk through givem apex-point
void VoronoiVertex::init(Point p, VertexStatus st, VertexType t, Point initDist) {
    init(p,st,t);
    init_dist(initDist);
}
/// set position, status, type, clerance-disk radius, and k3-side
void VoronoiVertex::init(Point p, VertexStatus st, VertexType t, Point initDist, double lk3) {
    init(p,st,t,initDist);
    k3 = lk3;
}

/// set in_queue false, and status to ::UNDECIDED
void VoronoiVertex::reset_status() {
    in_queue = false;
    status = UNDECIDED;
}
void VoronoiVertex::set_alfa(const Point& dir) {
    alfa = numeric::diangle(dir.x,dir.y);
}
/// initialize clerance-disk
void VoronoiVertex::init_dist(const Point& p) { r = dist(p); }
/// return distance to a point from this vertex
double VoronoiVertex::dist(const Point& p) const { return (position-p).norm(); }
/// set clearance-disk to zero
void VoronoiVertex::zero_dist() {r=0;}
/// return clearance disk-radius
double VoronoiVertex::dist() const { return r; }
/// in-circle predicate 
double VoronoiVertex::in_circle(const Point& p) const {
    //if ( r==0 && dist(p) == 0 ) 
    //    return -1;
    //else
        return dist(p) - r; 
}
/// reset the index count
void VoronoiVertex::reset_count() { count = 0; }


} // end ovd namespace
