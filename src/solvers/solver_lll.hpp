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

using namespace ovd::numeric; // sq() chop() determinant()

namespace ovd {

/// line-line-line vertex positioner.
/// solves 3x3 system.
class LLLSolver : public Solver {
public:

//  a1 x + b1 y + c1 + k1 t = 0
//  a2 x + b2 y + c2 + k2 t = 0
//  a3 x + b3 y + c3 + k3 t = 0
//
// or in matrix form
//
//  ( a1 b1 k1 ) ( x )    ( c1 )
//  ( a2 b2 k2 ) ( y ) = -( c2 )          Ax = b
//  ( a3 b3 k3 ) ( t )    ( c3 )
//
//  Cramers rule x_i = det(A_i)/det(A)
//  where A_i is A with column i replaced by b
            
int solve( Site* s1, double k1, 
           Site* s2, double k2, 
           Site* s3, double k3, std::vector<Solution>& slns ) {
    assert( s1->isLine() && s2->isLine() && s3->isLine() );
    
    std::vector< Eq<qd_real> > eq; // equation-parameters, in quad-precision
    boost::array<Site*,3> sites = {{s1,s2,s3}};    
    boost::array<double,3> kvals = {{k1,k2,k3}};
    for (unsigned int i=0;i<3;i++)
        eq.push_back( sites[i]->eqp_qd( kvals[i] ) );
    
    unsigned int i = 0, j=1, k=2;
    qd_real d = chop( determinant( eq[i].a, eq[i].b, eq[i].k, 
                                            eq[j].a, eq[j].b, eq[j].k, 
                                            eq[k].a, eq[k].b, eq[k].k ) ); 
    if (d != 0) {
        qd_real t = determinant(  eq[i].a, eq[i].b, -eq[i].c, 
                                  eq[j].a, eq[j].b, -eq[j].c, 
                                  eq[k].a, eq[k].b, -eq[k].c ) / d ; 
        if (t >= 0) {
            qd_real sol_x = determinant(  -eq[i].c, eq[i].b, eq[i].k, 
                                          -eq[j].c, eq[j].b, eq[j].k, 
                                          -eq[k].c, eq[k].b, eq[k].k ) / d ; 
            qd_real sol_y = determinant(  eq[i].a, -eq[i].c, eq[i].k, 
                                          eq[j].a, -eq[j].c, eq[j].k, 
                                          eq[k].a, -eq[k].c, eq[k].k ) / d ; 
            slns.push_back( Solution( Point( to_double(sol_x), to_double(sol_y) ), to_double(t), k3 ) ); // kk3 just passes through without any effect!?
            return 1;
        }
    } else {
        std::cout << "WARNING: LLLSolver determinant==0 ! no solutions. d= " << d <<"\n";
        std::cout << " 0 : " << eq[0].a << " " << eq[0].b << " " << eq[0].c << " " << eq[0].k << "\n";
        std::cout << " 1 : " << eq[1].a << " " << eq[1].b << " " << eq[1].c << " " << eq[1].k << "\n";
        std::cout << " 2 : " << eq[2].a << " " << eq[2].b << " " << eq[2].c << " " << eq[2].k << "\n";
        std::cout << " 0==1? " << (eq[0]==eq[1]) << "\n";
        std::cout << "da = " << (eq[0].a-eq[1].a) << "\n";
        std::cout << "db = " << (eq[0].b-eq[1].b) << "\n";
        std::cout << "dc = " << (eq[0].c-eq[1].c) << "\n";
        
    }
    return 0; // no solution if determinant zero, or t-value negative
}

};


} // ovd
#endif
