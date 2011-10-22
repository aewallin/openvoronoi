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

#include "voronoivertex.hpp"


namespace ovd {

class VoronoiDiagram;
//class Site;

/// Calculates the (x,y) position of vertices in a voronoi diagram
class VertexPositioner {
public:
    VertexPositioner(VoronoiDiagram* vodi): vd(vodi) {}
    /// calculate the position of a new voronoi-vertex lying on the given edge
    /// with the given new site 
    Point position( HEEdge e, Site* s);
private:
    inline double sq(double x) {return x*x;}
    /// point-point-point positioner
    Point ppp_solver(const Point& p1, const Point& p2, const Point& p3);
    int lll_solver(Site* s1, Site* s2, Site* s3); // linear 3x3 system
    
    int solver(Site* s1, Site* s2, Site* s3,  double solns[][3] ); 
    int qqq_solver( double l0[], double l1[], int xi, int yi, int ti, double xk, double yk, int kk, double rk, double solns[][3]);
    
    // int qll_solver( double l0[], double l1[], int xi, int yi, int ti, double xk, double yk, int kk, double rk);
    int qll_solve( double a0, double b0, double c0, double d0, 
                      double e0, double f0, double g0, 
                      double a1, double b1, 
                      double a2, double b2, 
                      double soln[][3]);
      
    int quadratic_roots(double a, double b, double c, double roots[]);


    Point position(Site* p1, Site* p2, Site* p3);
    //virtual Point position(PointSite* s1, PointSite* s2, PointSite* s3);
    
// geometry-checks
    bool check_on_edge(HEEdge e, const Point& p);
    bool check_in_edge(HEEdge e, const Point& p, HEVertex v);
    bool check_far_circle(const Point& p);
    bool check_dist(HEEdge e, const Point& p, HEVertex v);
    bool equal(double d1, double d2);
    double error(HEEdge e, const Point& p, HEVertex v);
// DATA
    VoronoiDiagram* vd;
    
};

}
#endif
