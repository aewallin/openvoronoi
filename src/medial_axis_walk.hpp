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

/// \brief Medial-axis point and associated clearance-disc radius.
struct MedialPoint {
    Point p; ///< position
    double clearance_radius; ///< clearance-disk radius
    /// \param pi position
    /// \param r radius
    MedialPoint(Point pi, double r): p(pi), clearance_radius(r) {}
};
typedef std::list<MedialPoint> MedialPointList; ///< list of points on the medial-axis
typedef std::list<MedialPointList> MedialChain; ///< a list of several lists
typedef std::list<MedialChain> MedialChainList; ///< a list of lists

// FIXME: MedialAxisWalk could probably be optimized to minimize rapid-traverses

/// \brief Walk along the medial-axis edges of a voronoi-diagram.
///
/// When we want a toolpath along the medial axis we first filter down the voronoi-diagram
/// with to a medial-axis with MedialAxisFilter and then use this class
/// to walk along the "valid" edges which remain.
///
/// Algorithm:
/// - first find one valid edge that has a degree-1 vertex (i.e. a suitable start point for the path)
/// - if there's only one choice for the next edge, go there
/// - if there are two choices, take one of the choices
/// - when done, find another valid start-edge.
class MedialAxisWalk {
public:
    /// \param gi vd-graph
    /// \param edge_pts subdivision for non-line edges
    MedialAxisWalk(HEGraph& gi, int edge_pts = 20): g(gi), _edge_points(edge_pts) {}

    /// run algorithm
    MedialChainList walk() {
        do_walk();
        return out;
    }
protected:
    void do_walk();
    void medial_axis_walk(HEEdge start);    
    bool valid_next_edge(HEEdge e);
    bool degree_one_source(HEEdge e);
    void append_edge(MedialChain& chain, HEEdge edge);
    bool next_edge(HEEdge e, HEEdge& next);
    void set_invalid(HEEdge e);
    bool find_start_edge(HEEdge& start);
    MedialChainList out; ///< output of algorithm
private:
    MedialAxisWalk(); // don't use.
    HEGraph& g; ///< original graph
    int _edge_points; ///< number of points to subdivide parabolas (non-line edges).

};

} // end namespace

// end file medial_axis.hpp
