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

#ifdef BUILD_PYTHON_LIB
#include <boost/python.hpp>
#endif

#include "graph.hpp"
#include "site.hpp"
#include "filter.hpp"

namespace ovd
{
/// \brief Filter that retains only edges between distinct polygons
///
/// \todo FIXME: not done yet!
struct island_filter : public Filter {
    island_filter() { }
    virtual bool operator()(const HEEdge& e) const {

        if (both_endpoints_positive(e)) // these are interior edges which we keep.
            return true;

        return false; // otherwise we keep the edge
    }
private:
    /// return true if this is an internal edge, i.e. both endpoints have a nonzero clearance-disk radius 
    bool both_endpoints_positive(HEEdge e) const {
        HEVertex src = g->source(e);
        HEVertex trg = g->target(e);
        return ( (*g)[src].dist()>0) && ( (*g)[trg].dist()>0);
    }
};



} // end namespace

// end file island_filter.hpp
