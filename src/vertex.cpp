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

inline double sq(double x) {return x*x;}

int VoronoiVertex::count = 0;

// the expected degree of a vertex. checked by topology-checker
VertexDegreeMap VoronoiVertex::expected_degree = boost::assign::map_list_of 
    (OUTER,4)     // special outer vertices
    (NORMAL,6)    // normal vertex in the graph
    (POINTSITE,0) // point site
    (ENDPOINT,6)  // end-point of line or arc
    (SEPPOINT,6)  // end-point of separator
    (SPLIT,4)     // split point, to avoid loops in delete-tree
    (APEX,4) ;    // apex point on quadratic bisector

VoronoiVertex::VoronoiVertex() {
    init();
    status = UNDECIDED;
    type = NORMAL;
}

/// construct vertex at position p with type t
VoronoiVertex::VoronoiVertex( Point p, VoronoiVertexStatus st) {
    init();
    position=p;
    status=st;
    type = NORMAL;
}
    
VoronoiVertex::VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t) {
    init();
    position=p;
    status=st;
    type=t;
}

/// set index, increase count, initialize in_queue to false.
void VoronoiVertex::init() {
    index = count;
    count++;
    in_queue = false;
    alfa=-1; // invalid/non-initialized alfa value
    null_face = std::numeric_limits<HEFace>::quiet_NaN();
    
}
/// set in_queue false, and status to UNDECIDED
void VoronoiVertex::reset() {
    in_queue = false;
    status = UNDECIDED;
}
void VoronoiVertex::set_alfa(const Point& dir) {
    alfa = numeric::diangle(dir.x,dir.y);
}

} // end ocl namespace
