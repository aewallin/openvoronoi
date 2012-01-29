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

#pragma once



namespace ovd {

class Site;


// abstract base-class for voronoi vertex-solvers
// the input to the solver is three Sites (s1,s2,s3) and three offset-directions (k1,k2,k3) .
class Solver {
public:
    virtual int solve(Site* s1, double k1, 
                           Site* s2, double k2, 
                           Site* s3, double k3, std::vector<Solution>& slns ) {return 0;}
    
    // this is used by PPPSolver, when no k-values are needed
    virtual int solve(Site* s1, 
                           Site* s2,  
                           Site* s3,  std::vector<Solution>& slns ) {return 0;}
    
    // used by alt_sep_solver
    virtual void set_type(int) {}
};


} //end namespace
