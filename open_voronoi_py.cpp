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
#include "point.hpp"

/*
 *  Python wrapping of voronoi diagram
 */

using namespace ovd;

namespace bp = boost::python;


BOOST_PYTHON_MODULE(openvoronoi) {

    bp::class_<VoronoiDiagram >("VoronoiDiagram_base")
    ;
    bp::class_< VoronoiDiagram_py, bp::bases<VoronoiDiagram> >("VoronoiDiagram")
        .def(bp::init<double, unsigned int>())
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site)
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site)
        .def("addLineSite1",  &VoronoiDiagram_py::insert_line_site1) // only find seed vertex
        .def("addLineSite2",  &VoronoiDiagram_py::insert_line_site2) // find seed, augment tree
        .def("addLineSite3",  &VoronoiDiagram_py::insert_line_site3) // find seed, augment tree
        .def("getGenerators",  &VoronoiDiagram_py::getGenerators)
        .def("getEdgesGenerators",  &VoronoiDiagram_py::getEdgesGenerators)
        .def("getVoronoiVertices",  &VoronoiDiagram_py::getVoronoiVertices)
        .def("getFarVoronoiVertices",  &VoronoiDiagram_py::getFarVoronoiVertices)
        .def("getFarRadius",  &VoronoiDiagram_py::get_far_radius)
        .def("getVoronoiEdges",  &VoronoiDiagram_py::getVoronoiEdges)
        //.def("getClosestFaceGenerator",  &VoronoiDiagram_py::getClosestFaceGenerator)
        //.def("getSeedVertex",  &VoronoiDiagram_py::getSeedVertex) 
        //.def("getSeedVertexLine",  &VoronoiDiagram_py::getSeedVertexLine) 
        //.def("getDeleteSet",  &VoronoiDiagram_py::getDeleteSet) //getDeleteSetLine(
        //.def("getDeleteSetLine",  &VoronoiDiagram_py::getDeleteSetLine)
        //.def("getDeleteEdges",  &VoronoiDiagram_py::getDeleteEdges)
        //.def("getModEdges",  &VoronoiDiagram_py::getModEdges)
        .def("__str__", &VoronoiDiagram_py::print)
        .def("version", &VoronoiDiagram_py::version)
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
        .value("PARABOLA", PARABOLA)
        .value("ELLIPSE", ELLIPSE)
        .value("HYPERBOLA", HYPERBOLA)
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
        .def("isRight", &Point::isRight)
        .def("__str__", &Point::str)
        .def_readwrite("x", &Point::x)
        .def_readwrite("y", &Point::y)
    ;
}

