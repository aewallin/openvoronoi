/* 
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include "numeric.hpp"

#include <qd/qd_real.h> 
#include <cmath>

namespace ovd {

namespace numeric {
    
    double diangle(double x, double y) {
        if (y >= 0)
            return (x >= 0 ? y/(x+y) : 1-x/(-x+y));
        else
            return (x < 0 ? 2-y/(-x-y) : 3+x/(x-y));
    }
    double diangle_x(double a) {
        return (a < 2 ? 1-a : a-3);
    }
    double diangle_y(double a) {
        return (a < 3 ? ((a > 1) ? 2-a : a) : a-4);
    }
    std::pair<double,double> diangle_xy(double a) {
        double x = diangle_x(a);
        double y = diangle_y(a);
        double norm = sqrt(x*x+y*y);
        return std::make_pair(x/norm,y/norm);
    }
    
} // numeric
} // ovd

