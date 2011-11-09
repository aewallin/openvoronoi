
#pragma once

#include <cmath>
#include "point.hpp"

namespace ovd {
/// solve system Ax = y by inverting A
/// x = Ainv * y
/// returns false if det(A)==0, i.e. no solution found
bool two_by_two_solver( const double& a, 
                        const double& b, 
                        const double& c,
                        const double& d,
                        const double& e,
                        const double& f,
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
    if ( fabs(det) < 1e-10 )
        return false;
    u = (1.0/det) * (d*e - b*f);
    v = (1.0/det) * (-c*e + a*f);
    return true;
}


/// find an intersection point in the XY-plane between two lines
/// first line:   p1 + v*(p2-p1)
/// second line:  p3 + t*(p4-p3)
/// sets (v,t) to the intersection point and returns true if an intersection was found 
bool xy_line_line_intersection( const Point& p1, const Point& p2, double& v,
                                const Point& p3, const Point& p4, double& t) {
    // p1 + v*(p2-p1) = p3 + t*(p4-p3)
    // =>
    // [ (p2-p1).x  -(p4-p3).x ] [ v ]  = [ (p3-p1).x ]
    // [ (p2-p1).y  -(p4-p3).y ] [ t ]  = [ (p3-p1).y ]
    return two_by_two_solver( (p2-p1).x , -(p4-p3).x , (p2-p1).y , -(p4-p3).y,  (p3-p1).x, (p3-p1).y, v, t);
}

}
