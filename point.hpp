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
#ifndef POINT_HPP
#define POINT_HPP

#include <string>
#include <iostream>

namespace ovd
{

/// \brief a point or vector in 2D with coordinates (x, y)
class Point {
    public:
        /// create a point at (0,0)
        Point(): x(0.0), y(0.0) {}
        /// create a point at (x,y)
        Point(double xi, double yi): x(xi), y(yi) {}
        /// create a point at p
        Point(const Point &p): x(p.x), y(p.y) {}
        virtual ~Point() {}        
        /// dot product
        double dot(const Point &p) const;
        /// 2D cross-product
        double cross(const Point& p) const;

        /// norm of vector, or distance from (0,0,0) to *this
        double norm() const; 
        /// squared norm (avoiding sqrt() might be faster in some cases)
        double norm_sq() const; 
        /// scales vector so that norm()==1.0
        void normalize();
        /// return perpendicular in the xy plane, rotated 90 degree to the left
        Point xyPerp() const;
        /// distance from this to p1-p2 line
        double distance_to_line(const Point &p1, const Point &p2) const;
        bool is_right(const Point &p1, const Point &p2) const;
        /// retrun true if Point within line segment p1-p2
        bool isInside(const Point& p1, const Point& p2) const;
        /// assignment
        Point &operator=(const Point &p);
        /// addition
        Point &operator+=(const Point &p);
        /// subtraction
        Point &operator-=(const Point &p);
        /// addition
        const Point operator+(const Point &p)const;
        /// subtraction
        const Point operator-(const Point &p) const;
        /// scalar multiplication
        Point &operator*=(const double &a);  // scalar multiplication with Point *= scalar
        /// Point * scalar
        const Point operator*(const double &a)const;     // Point*scalar 
        /// equality
        bool operator==(const Point &p) const;
        /// inequality
        bool operator!=(const Point &p) const;
        /// string repr
        friend std::ostream& operator<<(std::ostream &stream, const Point &p);
        /// string repr
        std::string str() const;
        /// X coordinate
        double x;
        /// Y coordinate
        double y;
};

/// scalar multiplication   scalar*Point
const Point operator*(const double &a, const Point &p);

} // end namespace
#endif
// end file point.hpp
