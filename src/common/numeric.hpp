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

#ifndef NUMERIC_HPP
#define NUMERIC_HPP

#include <vector>
#include <qd/qd_real.h> 

namespace ovd {

// this namespace holds general numerical functions that are not specific
// to voronoi-diagrams and may be useful elsewhere too
namespace numeric {
    //double chop8(double a);
    double chop(double val);
    double chop(double val, double tolerance);
    qd_real chop(qd_real val);
    
    template<class Scalar>
    Scalar sq( Scalar x) {return x*x;}
    
    /// solve quadratic eqn: a*x*x + b*x + c = 0
    /// returns real roots (0, 1, or 2) as vector
    template<class Scalar>
    std::vector<Scalar>  quadratic_roots(Scalar a, Scalar b, Scalar c) {
        std::vector<Scalar> roots;
        if ((a == 0) and (b == 0)) {
            //std::cout << " quadratic_roots() a == b == 0. no roots.\n";
            return roots;
        }
        if (a == 0) {
            roots.push_back( -c / b );
            return roots;
        }
        if (b == 0) {
            Scalar sqr = -c / a;
            if (sqr > 0) {
                roots.push_back( sqrt(sqr) );
                roots.push_back( -roots[0] );
                return roots;
            } else if (sqr == 0) {
                roots.push_back( Scalar(0) );
                return roots;
            } else {
                //std::cout << " quadratic_roots() b == 0. no roots.\n";
                return roots;
            }
        }
        Scalar disc = chop(b*b - 4*a*c); // discriminant, chop!
        if (disc > 0) {
            Scalar q;
            if (b > 0)
                q = (b + sqrt(disc)) / -2;
            else
                q = (b - sqrt(disc)) / -2;
            roots.push_back( q / a );
            roots.push_back( c / q ); 
            return roots;
        } else if (disc == 0) {
            roots.push_back( -b / (2*a) );
            return roots;
        }
        //std::cout << " quadratic_roots() disc < 0. no roots. disc= " << disc << "\n";
        return roots;
    }
    
    template <class Scalar>
    inline Scalar determinant( Scalar a, Scalar b, Scalar c,
                        Scalar d, Scalar e, Scalar f,
                        Scalar g, Scalar h, Scalar i ) {
        return a*(e*i-h*f)-b*(d*i-g*f)+c*(d*h-g*e);
    }
    
    double diangle(double x, double y);
    double diangle_x(double a);
    double diangle_y(double a);
    std::pair<double,double> diangle_xy(double a);
    bool diangle_bracket(double less, double a, double more);
    double diangle_mid(double alfa1, double alfa2);
    
} // numeric
} // ovd

#endif
