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

#ifndef VERTEX_POSITIONER_HPP
#define VERTEX_POSITIONER_HPP

#include "voronoidiagram_graph.hpp"


namespace ovd {

class VoronoiDiagram;

class VertexPositioner {
public:
    inline double sq(double x) {return x*x;}
    Point position(const Point& p1, const Point& p2, const Point& p3);

    // signature: edge, new_site 
    Point position(VoronoiDiagram* vd, HEEdge e, HEVertex v);

};

}
#endif
