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
    VertexPositioner(VoronoiDiagram* vodi): vd(vodi) {}
    
    // signature: edge, new_site 
    Point position( HEEdge e, HEVertex v);
private:
    Point position(const Point& p1, const Point& p2, const Point& p3);


    bool check_on_edge(HEEdge e, const Point& p);
    bool check_in_edge(HEEdge e, const Point& p, HEVertex v);
    bool check_far_circle(const Point& p);
    bool check_dist(HEEdge e, const Point& p, HEVertex v);
    bool equal(double d1, double d2);
    double error(HEEdge e, const Point& p, HEVertex v);
    
    VoronoiDiagram* vd;
    inline double sq(double x) {return x*x;}
};

}
#endif
