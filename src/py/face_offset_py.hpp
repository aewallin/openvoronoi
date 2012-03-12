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

#include "offset2.hpp"
#include "voronoidiagram.hpp"

namespace ovd
{

/// \brief python wrapper for FaceOffset
class FaceOffset_py : public FaceOffset {
public:
    /// ctor
    FaceOffset_py(HEGraph& gi): FaceOffset(gi) { }
    /// return list of offset elements
    boost::python::list offset_py(double t) {
        offset(t);
        boost::python::list py_offsets;
        BOOST_FOREACH( FaceOffsetLoop loop, offset_list ) { // loop through each loop
            HEEdge current_edge = loop.start_edge;
            boost::python::list py_loop;
            // add the first point to the loop.
            boost::python::list pt;
            pt.append( g[current_edge].point(loop.t) );
            pt.append( -1 ); // radius, center, cw
            py_loop.append(pt);
                
            BOOST_FOREACH( FaceEdge face_edge, loop.face_edges ) { //loop through each face
                Site* s = g[face_edge.f].site;
                Ofs* o = s->offset( g[current_edge].point(loop.t), g[face_edge.next_edge].point(loop.t) ); // ask the Site for offset-geometry here.
                bool cw(true);
                if (!s->isLine() ) // point and arc-sites produce arc-offsets, for which cw must be set.
                    cw = find_cw( o->start(), o->center(), o->end() ); // figure out cw or ccw arcs?
                boost::python::list lpt;
                    lpt.append( g[face_edge.next_edge].point(loop.t) );
                    lpt.append( o->radius() );
                    lpt.append( o->center() );
                    lpt.append( cw );
                py_loop.append(lpt);
                current_edge = face_edge.next_edge;
            }
            py_offsets.append( py_loop );
        }
        return py_offsets;
    }
private:
    FaceOffset_py(); // don't use.
};


} // end namespace

// end face_offset_py.hpp

