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
    
namespace solvers {
/*! 
 * \namespace ovd::solvers
 * \brief Voronoi vertex position solvers
 */

/// \brief abstract base-class for voronoi vertex position solvers
///
/// The input to the solver is three Sites (s1,s2,s3) and three offset-directions (k1,k2,k3).
/// The optput is a vector with one or more Solution.
class Solver {
public:
    Solver() { debug = false; type = 0; silent = true; }

    /// virtual dtor required for correct destruction of derived classes
    virtual ~Solver() { }
    
    /// \brief solve for position of VoronoiVertex with given adjacent sites and directions
    ///
    /// \param s1 first adjacent Site
    /// \param k1 direction from \a s1 to new VoronoiVertex
    /// \param s2 second adjacent Site
    /// \param k2 direction from \a s2 to new VoronoiVertex
    /// \param s3 third adjacent Site
    /// \param k3 direction from \a s3 to new VoronoiVertex
    /// \param slns Solution vector, will be updated by Solver
    virtual int solve(Site* s1, double k1, 
                           Site* s2, double k2, 
                           Site* s3, double k3, std::vector<Solution>& slns ) =0;

    /// used by alt_sep_solver
    virtual void set_type(int t) {type=t;}
    /// set the debug mode to \a b
    void set_debug(bool b) {debug=b;}
    /// no warnings/messages to stdout will be written, if silent is set true.
    void set_silent(bool b) {silent=b;}
protected:
    /// flag for debug output
    bool debug;
    /// separator case type.
    /// - type = 0 means l3 / p1 form a separator
    /// - type = 1 means l3 / p2 form a separator
    int  type;
    bool silent; ///< suppress all warnings or other stdout output
};

} // solvers
} //end ovd namespace
