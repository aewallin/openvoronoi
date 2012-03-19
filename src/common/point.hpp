/*  
 *  Copyright 2010-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include <string>
#include <iostream>

namespace ovd
{

/// \brief a point or vector in 2D with coordinates (x, y)
class Point {
    public:
        Point();
        Point(double xi, double yi);
        Point(const Point &p);
        virtual ~Point();
        
        double dot(const Point &p) const;
        double cross(const Point& p) const;
        double norm() const; 
        double norm_sq() const; 
        void normalize();
        Point xy_perp() const;
        bool is_right(const Point &p1, const Point &p2) const;
        
        Point &operator=(const Point &p);           ///< assignment
        Point &operator+=(const Point &p);          ///< addition
        Point &operator-=(const Point &p);          ///< subtraction
        const Point operator+(const Point &p)const; ///< addition
        const Point operator-(const Point &p) const;///< subtraction
        Point &operator*=(const double &a);         ///< scalar multiplication with Point *= scalar
        const Point operator*(const double &a)const;///< Point*scalar 
        bool operator==(const Point &p) const;      ///< equality
        bool operator!=(const Point &p) const;      ///< inequality
        
        friend std::ostream& operator<<(std::ostream &stream, const Point &p); ///< string repr
        std::string str() const; ///< string repr
        
        double x;   ///< X coordinate
        double y;   ///< Y coordinate
};

/// scalar multiplication   scalar*Point
const Point operator*(const double &a, const Point &p);

} // end ovd namespace

// end file point.hpp
