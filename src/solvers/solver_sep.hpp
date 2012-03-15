/*  
 *  Copyright 2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#pragma once

#include "common/point.hpp"
#include "common/numeric.hpp"

using namespace ovd::numeric; // sq() chop() determinant()

namespace ovd {
namespace solvers {
    
// this solver is called when we want to position a vertex on a SEPARATOR edge
// a SEPARATOR edge exists between a LineSite and one of its PointSite end-points
// the input sites are thus s1=LineSite and s2=PointSite  (if arcs are supported in the future then s1=ArcSite is possible)
// s3 can be either a LineSite or a PointSite (of arcs are supported, s3=ArcSite is possible)
//
//  s1 (LineSite) offset eq. is     a1 x + b1 y + c1 + k1 t = 0   (1)
//  s2 (PointSite) offset eq. is    (x-x2)^2 + (y-y2)^2 = t^2     (2)
// 
// Two possibilities for s3:
//   s3 (LineSite)   a3 x + b3 y + c3 + k3 t = 0
//   s3 (PointSite)  (x-x3)^2 + (y-y3)^2 = t^2 
// 
// This configuration constrains the solution to lie on the separator edge.
// The separator is given by
// SEP = p2 + t* sv
// where p2 is the location of s2, and the separator direction sv is 
// sv = (-a1,-b1)   if k1=-1
// sv = (a1,b1)   if k1=-1 
// thus points on the separator are located at:
//
//  x_sep = x2 + t*sv.x
//  y_sep = y2 + t*sv.y
//
//  This can be inserted into (1) or (2) above, which leads to a linear equation in t.
//
//  Insert into (1):
//    a3 (x2 + t*sv.x) + b3 (y2 + t*sv.y) + c3 + k3 t = 0
//      ==>
//          t = -( a3*x2 + b3*y2 + c3 ) / (sv.x*a3 + sv.y*b3 + k3)
//
//  Insert into (2):
//    (x2 + t*sv.x-x3)^2 + (y2 + t*sv.y-y3)^2 = t^2
//    ==> (using dx= x2-x3 and dy = x2-x3)
//   t^2 (sv.x^2 + sv.y^1 - 1)  + t (2*dx*sv.x + 2*dy*sv.y) + dx^2 + dy^2 = 0
//    ==>  (since sv is a unit-vector sv.x^2 + sv.y^1 - 1 = 0)
//         t = - (dx^2+dy^2) / (2*(dx*sv.x + dy*sv.y))
//
//  FIXME: what happens if we get a divide by zero situation ??
//

/// \brief ::SEPARATOR Solver
class SEPSolver : public Solver {
public:

int solve( Site* s1, double k1, 
           Site* s2, double k2, 
           Site* s3, double k3, std::vector<Solution>& slns ) {
    assert( s1->isLine() && s2->isPoint() );
    assert(s3->isLine());
    if (debug) 
        std::cout << "SEPSolver.\n";
    
    // separator direction
    Point sv(-s1->a(),-s1->b());
    if (debug) std::cout << " SEPSolver sv= "<< sv << "\n";
    
    double tsln = -(s3->a()*s2->x()+s3->b()*s2->y()+s3->c()) / ( sv.x*s3->a() + sv.y*s3->b() + k3  );

    Point psln = Point(s2->x(), s2->y() ) + tsln * sv;
    slns.push_back( Solution( psln, tsln, k3 ) );
    return 1;
}

};

} // solvers
} // ovd
