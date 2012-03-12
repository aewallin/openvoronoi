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

#include "common/point.hpp"

namespace ovd {
namespace solvers {
    
/// \brief a new vertex position solution (position, offset-distance, side)
///
/// includes the offset-distamce t, and the offset direction k3
struct Solution {
    /// \param pt vertex position
    /// \param tv clearance-disk radius
    /// \param k offset direction
    Solution(Point pt, double tv, double k) : p(pt), t(tv), k3(k) {}
    /// position
    Point p;
    /// clearance-disk radius
    double t;
    /// offfset direction to third adjacent Site
    double k3;
};

} // solvers
} //end ovd namespace
