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

#include <boost/python.hpp>

#include "voronoidiagram_py.hpp"  

#include "common/point.hpp"
#include "face_offset_py.hpp"
#include "island_filter.hpp"
#include "medial_axis_walk_py.hpp"
#include "offset_py.hpp"
#include "polygon_interior.hpp"
#include "utility/vd2svg.hpp"
#include "version.hpp"

/*
 *  Boost::Python wrapping of voronoi diagram and related classes.
 */

using namespace ovd;

namespace bp = boost::python;

BOOST_PYTHON_MODULE(openvoronoi) {
    bp::def("version", version);
    bp::def("build_type", build_type);
    bp::def("vd2svg", vd2svg);
    
    bp::class_<VoronoiDiagram >("VoronoiDiagram_base", bp::no_init)
    ;
    bp::class_<HEGraph>("Graph")
    ;
    bp::class_< VoronoiDiagram_py, bp::bases<VoronoiDiagram> >("VoronoiDiagram", bp::no_init)
        .def(bp::init<double, unsigned int>())
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site1 ) // (point)
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site2 ) // (point, step)
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site2 ) // takes two arguments
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site3 ) // takes three arguments (idx1, idx2, step)
        .def("getGenerators",  &VoronoiDiagram_py::getGenerators)
        .def("getEdgesGenerators",  &VoronoiDiagram_py::getEdgesGenerators)
        .def("getVoronoiVertices",  &VoronoiDiagram_py::getVoronoiVertices)
        .def("getFaceVertices",  &VoronoiDiagram_py::get_face_vertices) 
        .def("getFarVoronoiVertices",  &VoronoiDiagram_py::getFarVoronoiVertices)
        .def("getFarRadius",  &VoronoiDiagram_py::get_far_radius)
        .def("getVoronoiEdges",  &VoronoiDiagram_py::getVoronoiEdges)
        .def("getVoronoiEdgesOffset",  &VoronoiDiagram_py::getVoronoiEdgesOffset)
        .def("numPointSites", &VoronoiDiagram_py::num_point_sites)
        .def("numLineSites", &VoronoiDiagram_py::num_line_sites)
        .def("numVertices", &VoronoiDiagram_py::num_vertices)
        .def("numFaces", &VoronoiDiagram_py::num_faces)
        .def("numSplitVertices", &VoronoiDiagram_py::num_split_vertices)
        .def("__str__", &VoronoiDiagram_py::print)
        .def("reset_vertex_count", &VoronoiDiagram_py::reset_vertex_count)
        .def("setEdgePoints", &VoronoiDiagram_py::set_edge_points)
        .def("setEdgeOffset", &VoronoiDiagram_py::set_null_edge_offset)
        .def("debug_on", &VoronoiDiagram_py::debug_on)
        .def("check", &VoronoiDiagram_py::check)
        .staticmethod("reset_vertex_count")
        .def("getStat", &VoronoiDiagram_py::getStat)
        .def("filterReset", &VoronoiDiagram_py::filter_reset)
        .def("getFaceStats", &VoronoiDiagram_py::getFaceStats)
        .def("getGraph", &VoronoiDiagram_py::get_graph_reference, bp::return_value_policy<bp::reference_existing_object>())
    ;
    
    bp::enum_<VoronoiVertexStatus>("VoronoiVertexStatus")
        .value("OUT", OUT)   
        .value("IN", IN)
        .value("UNDECIDED", UNDECIDED)
        .value("NEW", NEW)
    ;
    bp::enum_<VoronoiVertexType>("VoronoiVertexType") // OUTER, NORMAL, POINTSITE, ENDPOINT, SEPPOINT, APEX, SPLIT};
        .value("OUTER", OUTER)   
        .value("NORMAL", NORMAL)
        .value("POINTSITE", POINTSITE)
        .value("ENDPOINT", ENDPOINT)
        .value("SEPPOINT", SEPPOINT)
        .value("APEX", APEX)
        .value("SPLIT", SPLIT)
    ;
    bp::enum_<VoronoiFaceStatus>("VoronoiFaceStatus")
        .value("INCIDENT", INCIDENT)
        .value("NONINCIDENT", NONINCIDENT)
    ;
    bp::enum_<VoronoiEdgeType>("VoronoiEdgeType")
        .value("LINE", LINE)
        .value("LINELINE", LINELINE)
        .value("PARA_LINELINE", PARA_LINELINE)
        .value("OUTEDGE", OUTEDGE)
        .value("PARABOLA", PARABOLA)
        .value("ELLIPSE", ELLIPSE)
        .value("HYPERBOLA", HYPERBOLA)
        .value("SEPARATOR", SEPARATOR)
        .value("LINESITE", LINESITE)
        .value("NULLEDGE", NULLEDGE)
    ;
    bp::class_<Point>("Point") 
        .def(bp::init<double, double>())
        .def(bp::init<Point>())
        .def(bp::other<double>() * bp::self)
        .def(bp::self *  bp::other<double>())
        .def(bp::self -= bp::other<Point>())
        .def(bp::self -  bp::other<Point>())
        .def(bp::self += bp::other<Point>())
        .def(bp::self +  bp::other<Point>())
        .def("norm", &Point::norm)
        .def("normalize", &Point::normalize)
        .def("dot", &Point::dot)
        .def("cross", &Point::cross)
        .def("is_ight", &Point::is_right)
        .def("__str__", &Point::str)
        .def_readwrite("x", &Point::x)
        .def_readwrite("y", &Point::y)
        .def_pickle(point_pickle_suite())
    ;
    bp::class_<Offset_py, boost::noncopyable >("Offset", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def("str", &Offset_py::print )
        .def("offset", &Offset_py::offset_py )
    ; 
    bp::class_<FaceOffset_py, boost::noncopyable >("FaceOffset", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def("offset", &FaceOffset_py::offset_py )
        .def("str", &FaceOffset_py::print )
    ; 
    bp::class_<PolygonInterior, boost::noncopyable >("PolygonInterior", bp::no_init)
        .def(bp::init<HEGraph&, bool>())
    ;
    bp::class_<MedialAxis, boost::noncopyable >("MedialAxis", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def(bp::init<HEGraph&, double>())
    ; 
    bp::class_<MedialAxisWalk_py, boost::noncopyable >("MedialAxisWalk", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def(bp::init<HEGraph&, int>())
        .def("walk", &MedialAxisWalk_py::walk_py)
    ;
    bp::class_<IslandFilter, boost::noncopyable >("IslandFilter", bp::no_init)
        .def(bp::init<HEGraph&>())
    ; 
}
