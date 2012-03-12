/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 * 
 *  Idea and code for point/line/arc voronoi-vertex positioning code by
 *  Andy Payne, andy "at" payne "dot" org, November, 2010
 *  see: http://www.payne.org/index.php/Calculating_Voronoi_Nodes
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
    
/// \brief line-line-line solver (parallel line-segment case)
///
/// solves 3x3 system.
class LLLPARASolver : public Solver {
public:

// parallel linesegment edge case.
//  a1 x + b1 y + c1 + k1 t = 0
//  a2 x + b2 y + c2 + k2 t = 0
//  a3 x + b3 y + c3 + k3 t = 0
//
// s1 and s2 are parallel, so they have a PARA_LINELINE edge between them
//
// this constrains the solution to lie on a line parallel to s1/s2
// passing through a point equidistant from s1/s2
// 
// equation of bisector is:
// ab x + bb y + cb = 0
// ab = a1
// bb = b1
// cb = (c1+c2)2
// all points on the bisector have a t value
// tb = fabs(c1-c2)/2
//
// find intersection of bisector and offset of third site
//  ab x + bb y + cb = 0
//  a3 x + b3 y + c3 + k3 tb = 0
//  or
//  ( ab  bb ) ( x ) = ( -cb )
//  ( a3  b3 ) ( y ) = ( -c3-k3*tb )
//


int solve( Site* s1, double k1, 
           Site* s2, double k2, 
           Site* s3, double k3, std::vector<Solution>& slns ) {
    if (debug)
        std::cout << "LLLPARASolver.\n";    
    assert( s1->isLine() && s2->isLine() && s3->isLine() );
    Eq<double> bis;
    bis.a = s1->a(); //eq[0].a;
    bis.b = s1->b(); //eq[0].b;
    double s2c = s2->c();
    // if s1 and s2 have opposite (a,b) normals, flip the sign of s2c
    if ( Point(s1->a(),s1->b()) == -1*Point(s2->a(),s2->b()) ) {
        s2c = -s2c;
    }
    
    bis.c = (s1->c() + s2c)/2;
    double tb = fabs( (s1->c() - s2c )/2 );
    
    if (debug) {
        std::cout << " s1 : " << s1->a() << " " << s1->b() << " " << s1->c() << " " << s1->k() << "\n";
        std::cout << " s2 : " << s2->a() << " " << s2->b() << " " << s2->c() << " " << s2->k() << "\n";
        std::cout << " s3 : " << s3->a() << " " << s3->b() << " " << s3->c() << " " << s3->k() << "\n";
        std::cout << " bis: " << bis.a << " " << bis.b << " " << bis.c << " \n";
    }
    double x,y;
    if ( two_by_two_solver(bis.a, bis.b, s3->a(), s3->b(), -bis.c, -s3->c()-k3*tb, x,y) ) {
        if (debug) std::cout << " Solution: t=" << tb << " " << Point( x, y ) << " k3=" << k3 << " \n";
        slns.push_back( Solution( Point( x, y ) , tb, k3 ) ); 
        return 1;
    } else {
        if (debug) std::cout << "LLLPARASolver. NO Solution!\n";    
        return 0;
    }
}

private:

/// solve 2z2 system Ax = y by inverting A
/// x = Ainv * y
/// returns false if det(A)==0, i.e. no solution found
bool two_by_two_solver( double a, 
                        double b, 
                        double c,
                        double d,
                        double e,
                        double f,
                        double& u,
                        double& v) {
    //  [ a  b ] [u] = [ e ]
    //  [ c  d ] [v] = [ f ]
    // matrix inverse is
    //          [ d  -b ]
    //  1/det * [ -c  a ]
    //  so
    //  [u]              [ d  -b ] [ e ]
    //  [v]  =  1/det *  [ -c  a ] [ f ]
    double det = a*d-c*b;
    if ( fabs(det) < 1e-15 ) // TOLERANCE!!
        return false;
    u = (1.0/det) * (d*e - b*f);
    v = (1.0/det) * (-c*e + a*f);
    return true;
}

};

} // solvers
} // ovd
