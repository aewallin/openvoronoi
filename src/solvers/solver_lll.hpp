/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 * 
 *  Idea and code for point/line/arc voronoi-vertex positioning code by
 *  Andy Payne, andy "at" payne "dot" org, November, 2010
 *  see: http://www.payne.org/index.php/Calculating_Voronoi_Nodes
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

#ifndef LLL_SOLVER_HPP
#define LLL_SOLVER_HPP

#include "common/point.hpp"
#include "common/numeric.hpp"

using namespace ovd::numeric; // sq() chop()

namespace ovd {

/// line-line-line vertex positioner.
/// solves 3x3 system.
class LLLSolver : public Solver {
public:

int solve(Site* s1, double k1, 
                           Site* s2, double k2, 
                           Site* s3, double k3, std::vector<Solution>& slns ) {
    assert( s1->isLine() && s2->isLine() && s3->isLine() );
    
    std::vector< Eq<qd_real> > eq; // equation-parameters, in quad-precision
    boost::array<Site*,3> sites = {{s1,s2,s3}};
    boost::array<double,3> kvals = {{k1,k2,k3}};
    for (unsigned int i=0;i<3;i++)
        eq.push_back( sites[i]->eqp_qd( kvals[i] ) );
    
    //std::cout << " lll_solver() k3= " << k3 << "\n";
    unsigned int i = 0, j=1, k=2;
    qd_real d = chop( ( eq[i].a*eq[j].b - eq[j].a*eq[i].b)*eq[k].k + 
                      (-eq[i].a*eq[k].b + eq[k].a*eq[i].b)*eq[j].k +  
                      ( eq[j].a*eq[k].b - eq[k].a*eq[j].b)*eq[i].k   ); // determinant
    if (d != 0) {
        qd_real t = ( (-eq[i].a*eq[j].b + eq[j].a*eq[i].b)*eq[k].c + 
                      ( eq[i].a*eq[k].b - eq[k].a*eq[i].b)*eq[j].c + 
                      (-eq[j].a*eq[k].b + eq[k].a*eq[j].b)*eq[i].c   )/d;
        if (t >= 0) {
            qd_real sol_x = ( ( eq[i].b*eq[j].c - eq[j].b*eq[i].c)*eq[k].k + 
                              (-eq[i].b*eq[k].c + eq[k].b*eq[i].c)*eq[j].k + 
                              ( eq[j].b*eq[k].c - eq[k].b*eq[j].c)*eq[i].k   )/d;
            qd_real sol_y = ( (-eq[i].a*eq[j].c + eq[j].a*eq[i].c)*eq[k].k + 
                              ( eq[i].a*eq[k].c - eq[k].a*eq[i].c)*eq[j].k + 
                              (-eq[j].a*eq[k].c + eq[k].a*eq[j].c)*eq[i].k   )/d;
            slns.push_back( Solution( Point( to_double(sol_x), to_double(sol_y) ), to_double(t), k3 ) ); // kk3 just passes through without any effect!?
            return 1;
        }
    }
    return 0; // no solution if determinant zero, or t-value negative
}

};


} // ovd
#endif
