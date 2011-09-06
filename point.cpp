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


namespace ovd
{

/*
Point::Point(const Point &p) {
    x=p.x;
    y=p.y;
}*/

double Point::norm() const {
    return sqrt( x*x + y*y  );
}

double Point::norm_sq() const {
    return x*x + y*y;
}


double Point::dot(const Point &p) const {
    return x * p.x + y * p.y ; 
}

void Point::normalize() {
    if (this->norm() != 0.0)
        *this *=(1/this->norm());
}

Point Point::xyPerp() const {
        return Point(-y, x);
}

double Point::xyDistanceToLine(const Point &p1, const Point &p2) const {
    // see for example
    // http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    if ((p1.x == p2.x) && (p1.y == p2.y)) {// no line in xy plane
        std::cout << "point.cpp: xyDistanceToLine ERROR!: can't calculate distance from \n";
        std::cout << "point.cpp: xyDistanceToLine ERROR!: *this ="<<*this <<" to line through\n";
        std::cout << "point.cpp: xyDistanceToLine ERROR!: p1="<<p1<<" and \n";
        std::cout << "point.cpp: xyDistanceToLine ERROR!: p2="<<p2<< "\n";
        std::cout << "point.cpp: xyDistanceToLine ERROR!: in the xy-plane\n";
        return -1;
    }
    else {
        Point v = Point(p2.y-p1.y, -(p2.x-p1.x) );
        v.normalize();
        Point r = Point(p1.x - x, p1.y - y );
        return fabs( v.dot(r));
    }
}


bool Point::isRight(const Point &p1, const Point &p2) const {
    // is Point right of line through points p1 and p2 ?, in the XY plane.
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

// scalar*Point
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

std::ostream& operator<<(std::ostream &stream, const Point& p) {
  stream << "(" << p.x << ", " << p.y << ")"; 
  return stream;
}

} // end namespace
// end file point.cpp
