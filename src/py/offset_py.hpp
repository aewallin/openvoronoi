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

//#ifndef PYTHON_DWA2002810_HPP
#include <boost/python.hpp>
//#endif

//#ifndef OFFSET_HPP
#include "offset.hpp"
//#endif
//#ifndef VORONOI_DIAGRAM_HPP
#include "voronoidiagram.hpp"
//#endif

namespace ovd
{

/// \brief python wrapper for Offset
class Offset_py : public Offset {
public:
    /// offset of graph \a gi
    Offset_py(HEGraph& gi): Offset(gi) { }
    /// return list of offsets at given offset distance \a t
    boost::python::list offset_py(double t) {
        offset(t);
        boost::python::list py_offsets;
        BOOST_FOREACH( OffsetLoop loop, offset_list ) { // loop through each loop
            boost::python::list py_loop;
            bool first = true;
            BOOST_FOREACH( OffsetVertex lpt, loop ) { //loop through each line/arc
                boost::python::list py_lpt;
                if (first) {
                    first = false;
                    py_lpt.append( lpt.p );
                    py_lpt.append( -1 );
                } else {
                    py_lpt.append( lpt.p );
                    py_lpt.append( lpt.r );
                    py_lpt.append( lpt.c );
                    py_lpt.append( lpt.cw );
                }
                py_loop.append( py_lpt );
            }
            py_offsets.append( py_loop );
        }
        return py_offsets;
    }
private:
    Offset_py(); // don't use.
};

} // end namespace
// end offset_py.hpp
