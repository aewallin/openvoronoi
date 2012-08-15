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
#include "site.hpp"
#include "offset.hpp"

namespace ovd
{

/*struct OffsetVertex {
    Point p;  ///< position (start)
    double r; ///< arc radius (line-vertex is indicated by radius of -1)
    Point c;  ///< arc center
    bool cw;  ///< clockwise (or not)
    HEFace f; ///< corresponding face in the vd-graph
    /// ctor
    OffsetVertex(Point pi, double ri, Point ci, bool cwi, HEFace fi): p(pi), r(ri), c(ci), cw(cwi), f(fi) {}
    /// ctor
    OffsetVertex(Point pi): p(pi), r(-1.), cw(false), f(0) {}
};

/// a single offset loop
struct OffsetLoop {
    std::list<OffsetVertex> vertices;
    double offset_distance;
    void push_back(OffsetVertex v) {
        vertices.push_back(v);
    }
};
*/

//typedef std::list<OffsetLoop> OffsetLoops;

class OffsetLoopCompare {
public:
    bool operator() (OffsetLoop l1, OffsetLoop l2) {
        return (l1.offset_distance > l2.offset_distance);
    }
};

class OffsetSorter {
public:
    OffsetSorter() {}
    void add_loop(OffsetLoop l) {
        all_loops.push_back(l);
    }
    void sort_loops() {
        BOOST_FOREACH( const OffsetLoop l, all_loops ) {
            dist_sorted_loops.insert(l);
        }
        BOOST_FOREACH( OffsetLoop l, dist_sorted_loops ) {
            std::cout << "loop at " << l.offset_distance << "\n";
        }
    }
protected:
    std::set<OffsetLoop, OffsetLoopCompare> dist_sorted_loops;
    OffsetLoops all_loops;
};


} // end ovd namespace
// end file offset_sorter.hpp
