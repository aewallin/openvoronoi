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

#include <boost/python.hpp>

#include "graph.hpp"
#include "site.hpp"

namespace ovd
{

struct interior_filter {
    interior_filter(HEGraph& gi) : g(gi) { }
    bool operator()(const HEEdge& e) const {
        if (g[e].type == LINESITE || g[e].type == NULLEDGE) 
            return true;
        
        // if polygon inserted ccw  as (id1->id2), then the linesite should occur on valid faces as id1->id2
        // for islands and the outside the edge is id2->id1
        HEFace f = g[e].face;
        Site* s = g[f].site;
        if ( s->isLine() && linesite_ccw(f) ) 
            return true;
        else if ( s->isPoint() ) {
            // we need to search for an adjacent linesite. (? can we have a situation where this fails?)
            HEEdge linetwin = find_adjacent_linesite(f);
            HEEdge twin = g[linetwin].twin;
            HEFace twin_face = g[twin].face;
            if (linesite_ccw(twin_face))
                return true;
        } 
        return false;
    }
    HEEdge find_adjacent_linesite( HEFace f ) const {
        HEEdge current = g[f].edge;
        HEEdge start = current;
        
        do {
            HEEdge twin = g[current].twin;
            if (twin != HEEdge() ) {
                HEFace twf = g[twin].face;
                if ( g[twf].site->isLine() )
                    return current;
            }
            current = g[current].next;
        } while(current!=start);
        return HEEdge();
    }
        
    bool linesite_ccw( HEFace f ) const {
        HEEdge current = g[f].edge;
        HEEdge start = current;
        do {
            if ( g[current].type == LINESITE && g[current].inserted_direction )
                return true;
            current = g[current].next;
        } while(current!=start);
        return false;
    }
    
    HEGraph& g;
};

/// \brief From a voronoi-diagram, generate offset curve(s).
class PolygonInterior {
public:
    PolygonInterior(HEGraph& gi): g(gi) {
        interior_filter f(g);
        interior_filter flt(g);
        g.filter_graph(flt);
    }


private:
    PolygonInterior(); // don't use.
    HEGraph& g; // original graph
    

};


} // end namespace

// end file polygon_interior.hpp
