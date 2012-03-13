/*  
 *  Copyright 2010-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include <boost/python.hpp>

#include "medial_axis_walk.hpp"
#include "voronoidiagram.hpp"

namespace ovd {
namespace pyovd {

/// \brief python wrapper for MedialAxisWalk
class MedialAxisWalk_py : public MedialAxisWalk {
public:
    /// create walk
    MedialAxisWalk_py(HEGraph& gi, int edge_pts = 20): MedialAxisWalk(gi, edge_pts) { }

    /// return list of medial-axis edges
    boost::python::list walk_py() {
        do_walk();
        boost::python::list py_walk;
        BOOST_FOREACH( MedialChain chain, out ) { // loop through each chain
            boost::python::list py_chain;
            BOOST_FOREACH( MedialPointList pt_list, chain ) { //loop through each medial-point list
                boost::python::list py_pt_list;
                BOOST_FOREACH( MedialPoint pt, pt_list ) { //loop through each medial-point
                    boost::python::list py_pt;
                    py_pt.append( pt.p );
                    py_pt.append( pt.clearance_radius );
                    py_pt_list.append( py_pt );
                }
                py_chain.append( py_pt_list );
            }
            py_walk.append( py_chain );
        }
        return py_walk;
    }
private:
    MedialAxisWalk_py(); // don't use.
};


} // pyovd
} // end ovd namespace

// end medial_axis_walk_py.hpp
