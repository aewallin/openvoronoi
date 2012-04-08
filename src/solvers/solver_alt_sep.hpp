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
//  s1 (LineSite) offset eq. is     a1 x + b1 y + c1 + k1 t = 0   
//  s2 (PointSite) offset eq. is    (x-x2)^2 + (y-y2)^2 = t^2     
// 
// Two possibilities for s3:
//   s3 (LineSite)   a3 x + b3 y + c3 + k3 t = 0           (1)
//   s3 (PointSite)  (x-x3)^2 + (y-y3)^2 = t^2             (2)
// 
// This configuration constrains the solution to lie on the separator edge.
// The separator is given by
// SEP = p2 + t* sv
// where p2 is the location of s2, and the separator direction sv is 
// sv = (-a1,-b1)   if k1=-1
// sv = (a1,b1)   if k1=+1 
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
/// \brief alternative ::SEPARATOR Solver
class ALTSEPSolver : public Solver {
public:
//virtual void set_type(int t) {type=t;}
int solve( Site* s1, double k1, 
           Site* s2, double k2, 
           Site* s3, double k3, std::vector<Solution>& slns ) {
    if (debug && !silent) 
        std::cout << "ALTSEPSolver.\n";
    Site* lsite;
    Site* psite;
    Site* third_site;
    double lsite_k,  third_site_k;
    
    if ( type == 0 ) {
        lsite = s3; lsite_k = k3;
        psite = s1; // psite_k = k1;    l3 / p1 form a separator
        third_site = s2;      third_site_k = 1;
    } else if ( type == 1 ) {
        lsite = s3; lsite_k = k3;
        psite = s2; // psite_k = k2;    l3 / p2 form a separator
        third_site = s1; third_site_k = 1; 
    } else {
        std::cout << "ALTSEPSolver FATAL ERROR! type not known.\n";
        exit(-1);
        return 0;
    }
    // separator direction
    Point sv = (k3 == - 1) ? Point(lsite->a(),lsite->b()) : Point(-lsite->a(),-lsite->b());
    
    if (debug && !silent) {
        std::cout << "ALTSEPSolver type="<< type <<"\n";
        std::cout << " s1= " << s1->str2() << "(k=" << k1<< ")\n";
        std::cout << " s2= " << s2->str2() << "(k=" << k2<< ")\n";
        std::cout << " s3= " << s3->str2() << "(k=" << k3<< ")\n";
        std::cout << " lsite_k=" << lsite_k << "\n";
        std::cout << " sv= " << sv << "\n";
    }

    
    // now we should have this:
    assert( lsite->isLine() && psite->isPoint() );

    double tsln(0);

    if ( third_site->isPoint() ) {
        double dx = psite->x() - third_site->x();
        double dy = psite->y() - third_site->y();
        if ( fabs(2*( dx*sv.x+dy*sv.y  )) > 0 ) {
            tsln = -(dx*dx+dy*dy) / (2*( dx*sv.x+dy*sv.y  )); // check for divide-by-zero?
        } else {
            //std::cout << " no solutions. (isPoint)\n";
            return 0;
        }
    } else if (third_site->isLine()) {
        if ( fabs(( sv.x*third_site->a() + sv.y*third_site->b() + third_site_k )) > 0 ) {
            tsln = -(third_site->a()*psite->x()+third_site->b()*psite->y()+third_site->c()) / 
                ( sv.x*third_site->a() + sv.y*third_site->b() + third_site_k );
        } else {
            //std::cout << " no solutions. (isLine)\n";
            return 0;
        }
    } else {
        assert(0);
        exit(-1);
    }
    Point psln = Point(psite->x(), psite->y() ) + tsln * sv;
    slns.push_back( Solution( psln, tsln, k3 ) );
    return 1;
}

};

} // solvers
} // ovd
