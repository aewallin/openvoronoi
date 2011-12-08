/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  
 *  This file is part of OpenCAMlib.
 *
 *  OpenCAMlib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  OpenCAMlib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with OpenCAMlib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <boost/python.hpp>

#include "voronoidiagram_py.hpp"  
#include "common/point.hpp"

/*
 *  Python wrapping of voronoi diagram
 */

using namespace ovd;

namespace bp = boost::python;

std::string ovd_revision() {
    return VERSION_STRING;
}

BOOST_PYTHON_MODULE(openvoronoi) {
    bp::def("revision", ovd_revision);
    
    bp::class_<VoronoiDiagram >("VoronoiDiagram_base")
    ;
    bp::class_< VoronoiDiagram_py, bp::bases<VoronoiDiagram> >("VoronoiDiagram")
        .def(bp::init<double, unsigned int>())
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site1 )
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site2 )
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site2 ) // takes one argument
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site3 ) // takes two arguments
        .def("getGenerators",  &VoronoiDiagram_py::getGenerators)
        .def("getEdgesGenerators",  &VoronoiDiagram_py::getEdgesGenerators)
        .def("getVoronoiVertices",  &VoronoiDiagram_py::getVoronoiVertices)
        .def("getFarVoronoiVertices",  &VoronoiDiagram_py::getFarVoronoiVertices)
        .def("getFarRadius",  &VoronoiDiagram_py::get_far_radius)
        .def("getVoronoiEdges",  &VoronoiDiagram_py::getVoronoiEdges)
        .def("numPointSites", &VoronoiDiagram_py::num_point_sites)
        .def("numLineSites", &VoronoiDiagram_py::num_line_sites)
        .def("numVertices", &VoronoiDiagram_py::num_vertices)
        .def("numSplitVertices", &VoronoiDiagram_py::num_split_vertices)
        .def("__str__", &VoronoiDiagram_py::print)
        .def("version", &VoronoiDiagram_py::version)
        .def("reset_vertex_count", &VoronoiDiagram_py::reset_vertex_count)
        .def("setEdgePoints", &VoronoiDiagram_py::set_edge_points)
        .def("debug_on", &VoronoiDiagram_py::debug_on)
        .def("check", &VoronoiDiagram_py::check)
        .staticmethod("reset_vertex_count")
        .def("getStat", &VoronoiDiagram_py::getStat)
        .def("getFaceStats", &VoronoiDiagram_py::getFaceStats)
    ;
    bp::enum_<VoronoiVertexStatus>("VoronoiVertexStatus")
        .value("OUT", OUT)   
        .value("IN", IN)
        .value("UNDECIDED", UNDECIDED)
        .value("NEW", NEW)
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
    ;
    bp::class_<Point>("Point") 
        .def(bp::init<double, double>())
        .def(bp::init<Point>())
        .def(bp::other<double>() * bp::self)
        .def(bp::self * bp::other<double>())
        .def(bp::self -= bp::other<Point>())
        .def(bp::self - bp::other<Point>())
        .def(bp::self += bp::other<Point>())
        .def(bp::self + bp::other<Point>())
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
}

