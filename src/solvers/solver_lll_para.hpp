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
    
/// \brief line-line-line Solver (parallel line-segment case)
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
    
    Eq<double> bisector;
    bisector.a = s1->a();
    bisector.b = s1->b();
    double s2c = s2->c();

    // if s1 and s2 have opposite (a,b) normals, flip the sign of s2c
    Point n0(s1->a(), s1->b());
    Point n1(s2->a(), s2->b());
    if (n0.dot(n1) < 0.f) 
	{
        s2c = -s2c;
    }
    
    bisector.c = (s1->c() + s2c)*0.5;
    double tb = 0.5*fabs(s1->c() - s2c); // bisector offset distance
    
    if (debug) {
        std::cout << " s1 : " << s1->a() << " " << s1->b() << " " << s1->c() << " " << s1->k() << "\n";
        std::cout << " s2 : " << s2->a() << " " << s2->b() << " " << s2->c() << " " << s2->k() << "\n";
        if ( s3->isLine() )
            std::cout << " s3 : " << s3->a() << " " << s3->b() << " " << s3->c() << " " << s3->k() << "\n";
        std::cout << " bisector: " << bisector.a << " " << bisector.b << " " << bisector.c << " \n";
    }
    //if ( s3->isLine() ) {
        double x,y;
        if ( two_by_two_solver(bisector.a, bisector.b, s3->a(), s3->b(), -bisector.c, -s3->c()-k3*tb, x,y) ) {
            Point psln(x, y);
            if (debug) std::cout << " Solution: t=" << tb << " " << psln << " k3=" << k3 << " \n";
            if ((s1->end() - s1->start()).cross(psln - s1->start()) * k1 < 0 ||
                (s2->end() - s2->start()).cross(psln - s2->start()) * k2 < 0 ||
                (s3->end() - s3->start()).cross(psln - s3->start()) * k3 < 0) {
                // solution lies on the wrong side from one of the lines
                return 0;
            } else {
                slns.push_back( Solution( Point( x, y ) , tb, k3 ) ); 
                return 1;
            }
        } else {
            if (debug) std::cout << "LLLPARASolver. NO Solution!\n";    
            return 0;
        }
    //}
    /*
    if ( s3->isArc() ) {
        // bisector ax + by + c = 0
        // all points are at offset tb from s1 and s2
        // find a point on the bisector which is also a distance tb from the circle
        // this point is a distance r+tb from the circle center (?)
        circle_line_intersection(bis.a, bis.b, bis.c, s3->x(), s3->y(), s3->r(), tb,slns);
        return 1;
    }*/
    return 0;
}

private:
/*
void circle_line_intersection(double a, double b, double c, 
         double cx, double cy, double r, double tb, std::vector<Solution>& slns) {
    // line ax+by+c = 0
    // circle (cx,cy) radius r
}*/
    
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
