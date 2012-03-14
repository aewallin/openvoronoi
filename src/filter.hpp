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

#include <string>
#include <iostream>

#include "graph.hpp"
#include "common/numeric.hpp"
#include "site.hpp"

namespace ovd
{

/// \brief base-class for voronoi-diagram filters
///
/// concrete sub-classes of Filter provide a predicate
/// for determining if the edge belongs to the filtered graph. 
class Filter {
public:
    Filter() {}
    /// set graph
    void set_graph(HEGraph* gi) {g=gi;}
    /// does this edge belong to the filtered graph?
    virtual bool operator()(const HEEdge& e) const =0; //{exit(-1); return true;}
protected:
    HEGraph* g; ///< vd-graph
};

} // end ovd namespace
// end file 
