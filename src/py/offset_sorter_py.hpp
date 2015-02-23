/*  
 *  Copyright 2015 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include "offset_sorter.hpp"
#include <boost/python.hpp>

namespace ovd
{


/// Python wrapper for OffsetSorter
class OffsetSorter_py : public OffsetSorter {
public:
    /// create walk
    OffsetSorter_py(HEGraph& gi): OffsetSorter(gi) { }
    
    /// return list of offsets as python list
    boost::python::list offset_list_py() {
        boost::python::list py_offsets;
        BOOST_FOREACH( MGVertex v, boost::vertices(g) ) { // loop through each loop
            boost::python::list py_loop;
            bool first = true;
            int vdeg =boost::out_degree( v, g );
            /*
            BOOST_FOREACH( Edge e, boost::out_edges( v, g ) ) {
                vdeg++;
            }*/
            
            
            BOOST_FOREACH( OffsetVertex lpt, g[v].vertices ) { //loop through each line/arc
                boost::python::list py_lpt;
                double offset_distance = g[v].offset_distance;
                if (first) {
                    first = false;
                    py_lpt.append( lpt.p ); // 0
                    py_lpt.append( -1 ); // 1
                    py_lpt.append( offset_distance ); // 2
                    py_lpt.append( vdeg ); // 3
                } else {
                    py_lpt.append( lpt.p ); // 0, position
                    py_lpt.append( lpt.r ); // 1, radius
                    py_lpt.append( lpt.c ); // 2, center
                    py_lpt.append( lpt.cw ); // 3, cw or ccw
                    py_lpt.append( lpt.f ); // 4, face
                    py_lpt.append( offset_distance ); // 5
                    py_lpt.append( vdeg ); // 6
                }
                py_loop.append( py_lpt );
            }
            py_offsets.append( py_loop );
        }
        return py_offsets;
    }
private:
    OffsetSorter_py(); // don't use.
};


} // end ovd namespace
// end file offset_sorter.hpp
