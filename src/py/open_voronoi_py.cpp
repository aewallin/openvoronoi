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
#include "medial_axis_walk_py.hpp"
#include "offset_py.hpp"
#include "offset_sorter_py.hpp"

#include "utility/vd2svg.hpp"
#include "version.hpp"

// filters:
#include "polygon_interior_filter.hpp"
#include "island_filter.hpp"
#include "medial_axis_filter.hpp"

#include "medial_axis_pocket_py.hpp"

/*
 *  Boost::Python wrapping of voronoi diagram and related classes.
 */

namespace ovd {
namespace pyovd {
/*! 
 * \namespace ovd::pyovd
 * \brief Python wrappers for OpenVoronoi
 */
 
namespace bp = boost::python;
/// create openvoronoi python module
BOOST_PYTHON_MODULE(openvoronoi) {
    bp::def("version", version);
    bp::def("build_type", build_type);
    bp::def("vd2svg", vd2svg);
    
    bp::class_<VoronoiDiagram , boost::noncopyable >("VoronoiDiagram_base", bp::no_init)
    ;
    bp::class_<HEGraph , boost::noncopyable >("Graph")
    ;
    bp::class_< VoronoiDiagram_py, boost::noncopyable, bp::bases<VoronoiDiagram> >("VoronoiDiagram", bp::no_init)
        .def(bp::init<double, unsigned int>())
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site1 ) // (point)
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site2 ) // (point, step)
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site2 ) // takes two arguments
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site3 ) // takes three arguments (idx1, idx2, step)
        .def("addArcSite",  &VoronoiDiagram_py::insert_arc_site ) // arc-site (idx1,idx2, center, cw?, step) 
        .def("addArcSite",  &VoronoiDiagram_py::insert_arc_site4 ) // arc-site (idx1,idx2, center, cw?, step) 
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
        .def("numArcSites", &VoronoiDiagram_py::num_arc_sites)
        .def("numVertices", &VoronoiDiagram_py::num_vertices)
        .def("numFaces", &VoronoiDiagram_py::num_faces)
        .def("numSplitVertices", &VoronoiDiagram_py::num_split_vertices)
        .def("__str__", &VoronoiDiagram_py::print)
        .def("reset_vertex_count", &VoronoiDiagram_py::reset_vertex_count)
        .def("setEdgePoints", &VoronoiDiagram_py::set_edge_points)
        .def("setEdgeOffset", &VoronoiDiagram_py::set_null_edge_offset)
        .def("debug_on", &VoronoiDiagram_py::debug_on)
        .def("set_silent", &VoronoiDiagram_py::set_silent)
        .def("check", &VoronoiDiagram_py::check)
        .staticmethod("reset_vertex_count")
        .def("getStat", &VoronoiDiagram_py::getStat)
        .def("filterReset", &VoronoiDiagram_py::filter_reset)
        .def("filter_graph", &VoronoiDiagram_py::filter) // "filter" is a built-in function in Python!
        .def("getFaceStats", &VoronoiDiagram_py::getFaceStats)
        .def("getGraph", &VoronoiDiagram_py::get_graph_reference, bp::return_value_policy<bp::reference_existing_object>())
    ;
    
    bp::enum_<VertexStatus>("VertexStatus")
        .value("OUT", OUT)   
        .value("IN", IN)
        .value("UNDECIDED", UNDECIDED)
        .value("NEW", NEW)
    ;
    bp::enum_<VertexType>("VertexType")
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
    bp::enum_<EdgeType>("EdgeType")
        .value("LINE", LINE)
        .value("LINELINE", LINELINE)
        .value("PARA_LINELINE", PARA_LINELINE)
        .value("OUTEDGE", OUTEDGE)
        .value("PARABOLA", PARABOLA)
        .value("ELLIPSE", ELLIPSE)
        .value("HYPERBOLA", HYPERBOLA)
        .value("SEPARATOR", SEPARATOR)
        .value("LINESITE", LINESITE)
        .value("ARCSITE", ARCSITE)
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
        .def("is_right", &Point::is_right)
        .def("xy_perp", &Point::xy_perp)
        .def("__str__", &Point::str)
        .def_readwrite("x", &Point::x)
        .def_readwrite("y", &Point::y)
        .def_pickle(point_pickle_suite())
    ;
// Offsetting
    bp::class_<Offset_py, boost::noncopyable >("Offset", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def("str", &Offset_py::print )
        .def("offset", &Offset_py::offset_py )
        .def("offset_loop_list", &Offset_py::offset_loop_list )
    ; 
    bp::class_< OffsetLoop  >("OffsetLoop")
    ;  
    bp::class_< OffsetSorter_py , boost::noncopyable  >("OffsetSorter", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def("add_loop", &OffsetSorter_py::add_loop )
        .def("sort_loops", &OffsetSorter_py::sort_loops )
        .def("get_loops", &OffsetSorter_py::offset_list_py )
    ;  
  
// Filters
    bp::class_< Filter, boost::noncopyable >(" Filter_base", bp::no_init) // pure virtual base class!
    ;
    bp::class_<polygon_interior_filter, bp::bases<Filter> >("PolygonInterior")
        .def(bp::init<bool>())
    ;
    bp::class_<island_filter,  bp::bases<Filter> >("IslandFilter")
    ; 
    
    bp::class_< medial_axis_filter, bp::bases<Filter> >("MedialAxis")
        .def(bp::init<double>())
    ; 
    
    
    bp::class_<MedialAxisWalk_py, boost::noncopyable >("MedialAxisWalk", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def(bp::init<HEGraph&, int>())
        .def("walk", &MedialAxisWalk_py::walk_py)
    ;
    
    bp::class_<medial_axis_pocket_py, boost::noncopyable >("MedialAxisPocket", bp::no_init)
        .def(bp::init<HEGraph&>())
        .def("run", &medial_axis_pocket_py::run)
        //.def("run2", &medial_axis_pocket_py::run2) 
        //.def("get_mic_list", &medial_axis_pocket_py::py_get_mic_list)
        .def("get_mic_components", &medial_axis_pocket_py::py_get_mic_components)
        .def("setWidth", &medial_axis_pocket_py::set_width)
        .def("debug", &medial_axis_pocket_py::set_debug)
    ;
}

} // pyovd namespace
} // ovd namespace
