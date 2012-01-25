/*  
 *  Copyright 2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

// Separator solver
class SEPSolver : public Solver {
public:

int solve( Site* s1, double k1, 
           Site* s2, double k2, 
           Site* s3, double k3, std::vector<Solution>& slns ) {

    assert( s1->isLine() && s2->isPoint() );
    // separator direction
    Point sv(0,0);
    if (k2 == -1) {
        sv.x = s1->a(); //l1.a
        sv.y = s1->b(); //l1.b
    } else {
        sv.x = -s1->a();
        sv.y = -s1->b();
    }

    double tsln(0);
    if ( s3->isPoint() ) {
        double dx = s2->x() - s3->x();
        double dy = s2->y() - s3->y();
        tsln = -(dx*dx+dy*dy) / (2*( dx*sv.x+dy*sv.y  )); // check for divide-by-zero?

    } else if (s3->isLine()) {
        tsln = -(s3->a()*s2->x()+s3->b()*s2->y()+s3->c()) / ( sv.x*s3->a() + sv.y*s3->b() + k3  );
    } else {
        assert(0);
        exit(-1);
    }
    Point psln = Point(s2->x(), s2->y() ) + tsln * sv;
    slns.push_back( Solution( psln, tsln, k3 ) );
    return 1;
}

};


} // ovd
