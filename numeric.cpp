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
 
    double chop8(double a) {
        double eps = 1e-8;
        if (fabs(a) < eps)
            return 0.0;
        else
            return a;
    }
    
    double chop(double val) {
        double _epsilon = 1e-10;
        if (fabs(val) < _epsilon) 
            return 0;
        else
            return val;
    }
    qd_real chop(qd_real val) {
        qd_real _epsilon = 1e-10;
        if (fabs(val) < _epsilon) 
            return qd_real(0);
        else
            return val;
    }

} // numeric
} // ovd

