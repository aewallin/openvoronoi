/*  
 *  Copyright 2010 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include <cmath> // sqrt, sin, cos, fabs
#include <cassert>
#include <sstream>

#include "point.hpp"

namespace ovd {
    
/// create a point at (0,0)
Point::Point(): x(0.0), y(0.0) {}
/// create a point at (x,y)
Point::Point(double xi, double yi) : x(xi), y(yi) {}
/// create a point at p
Point::Point(const Point &p): x(p.x), y(p.y) {}
/// dtor
Point::~Point() {}

/// norm of vector, or distance from (0,0,0) to *this
double Point::norm() const {
    return sqrt( x*x + y*y  );
}
/// squared norm (avoiding sqrt() might be faster in some cases)
double Point::norm_sq() const {
    return x*x + y*y;
}

/// dot product
double Point::dot(const Point &p) const {
    return x * p.x + y * p.y ; 
}
/// 2D cross-product
double Point::cross(const Point &p) const {
    return x * p.y - y * p.x ; 
}
/// scales vector so that norm()==1.0
void Point::normalize() {
    if (this->norm() != 0.0)
        *this *=(1/this->norm());
}
/// return perpendicular in the xy plane, rotated 90 degree to the left
Point Point::xy_perp() const {
    return Point(-y, x);
    // 2D rotation matrix:
    //   cos   -sin
    //   sin   cos
    // for theta = 90
    //   0   -1   ( x )
    //   1    0   ( y )  = ( -y  x )
    
}

/// is this Point right of line through points \a p1 and \a p2 ?
bool Point::is_right(const Point &p1, const Point &p2) const {
    // this is an ugly way of doing a determinant
    // should be prettyfied sometime...
    /// \todo FIXME: what if p1==p2 ? (in the XY plane)
    double a1 = p2.x - p1.x;
    double a2 = p2.y - p1.y;
    double t1 = a2;
    double t2 = -a1;
    double b1 = x - p1.x;
    double b2 = y - p1.y;

    double t = t1 * b1 + t2 * b2;
    if (t > 0.0) 
        return true;
    else
        return false;    
}

/* **************** Operators ***************  
 *  see
 *  http://www.cs.caltech.edu/courses/cs11/material/cpp/donnie/cpp-ops.html
*/

Point& Point::operator=(const Point &p) {
    if (this == &p)
        return *this;
    x=p.x;
    y=p.y;
    return *this;
}

// Point*scalar multiplication
Point& Point::operator*=(const double &a) {
    x*=a;
    y*=a;
    return *this;
}

Point& Point::operator+=(const Point &p) {
    x+=p.x;
    y+=p.y;
    return *this;
}

Point& Point::operator-=(const Point &p) {
    x-=p.x;
    y-=p.y;
    return *this;
}

const Point Point::operator+(const Point &p) const {
    return Point(*this) += p;
}

const Point Point::operator-(const Point &p) const {
    return Point(*this) -= p;
}

const Point Point::operator*(const double &a) const {
    return Point(*this) *= a;
}

/// scalar*Point
const Point operator*(const double &a, const Point &p) {
    return Point(p) *= a;
}

bool Point::operator==(const Point &p) const {
    return (this == &p) || (x==p.x && y==p.y ); 
}

bool Point::operator!=(const Point &p) const {
    return !(*this == p);
}

std::string Point::str() const {
    std::ostringstream o;
    o << *this;
    return o.str();
}

/// string output for Point
std::ostream& operator<<(std::ostream &stream, const Point& p) {
  stream << "(" << p.x << ", " << p.y << ")"; 
  return stream;
}

} // end namespace
// end file point.cpp
