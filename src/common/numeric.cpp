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
#include <cassert>

namespace ovd {

namespace numeric {
    
    /*
    double chop8(double a) {
        double eps = 1e-8;
        if (fabs(a) < eps)
            return 0.0;
        else
            return a;
    }*/
    double chop(double val, double tol) {
        double _epsilon = tol;
        if (fabs(val) < _epsilon) 
            return 0;
        else
            return val;
    }
    
    double chop(double val) {
        double _epsilon = 1e-10;
        if (fabs(val) < _epsilon) 
            return 0;
        else
            return val;
    }
    qd_real chop(qd_real val) {
        qd_real _epsilon = 1e-20; // should leave 47bits of precision
        if (fabs(val) < _epsilon) 
            return qd_real(0);
        else
            return val;
    }


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
    // return true if a lies in [less,more]
    bool diangle_bracket(double less, double a, double more) {
        if (less==more) {
            return false;
        }else if (less<=more) { // normal case..
            if ( (less<=a) && (a<more) )
                return true;
        } else if (less>more) { // we cross zero
            if ( ((less<=a) && (a<=4)) || ((0<=a) && (a<more)) )
                return true;
        } else {
            assert(0);
            return false;
        }
        
        return false;
    }
    // return average of input angles
    double diangle_mid(double alfa1, double alfa2) {
        if (alfa1<=alfa2)
            return (alfa1+alfa2)/2;
        else {
            double opposite_mid = alfa2 + (alfa1-alfa2)/2;
            double mid = opposite_mid + 2;
            if (mid>4)
                mid=mid-4;
            assert( (0<=mid) && (mid<=4) );
            return mid;
        }
    }
} // numeric
} // ovd

