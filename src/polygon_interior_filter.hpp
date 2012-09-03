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

#include "filter.hpp"
#include "graph.hpp"
#include "site.hpp"

namespace ovd
{

/// \brief Filter for retaining voronoi-diagram inside a polygon
///
/// this filter sets the valid-property of edges
/// all interior edges are marked valid=true
/// all exterior edges are marked valid=false
///
/// a polygon/pocket boundary shoud be specified in CW order
/// islands within the polygon should be specified in CCW order
class polygon_interior_filter : public Filter {
public:
    /// \brief create a polygon interior Filter with given \a side
    /// \param side set true (false) for polygons inserted in CW (CCW) order and islands inserted in CCW (CW) order.
    polygon_interior_filter( bool side=true) : _side(side) { }
    
    /// determine if an edge is valid or not
    virtual bool operator()(const HEEdge& e) const {
        
        if ( (*g)[e].type == LINESITE || (*g)[e].type == NULLEDGE) 
            return true;
        
        // if polygon inserted ccw  as (id1->id2), then the linesite should occur on valid faces as id1->id2
        // for islands and the outside the edge is id2->id1
        
        HEFace f = (*g)[e].face;
        Site* s = (*g)[f].site;
        if ( s->isLine() && linesite_ccw(f) ) 
            return true;
        else if ( s->isPoint() ) {
            //HEVertex site_vertex = s->vertex();
            //std::cout << "PointSite type: " << (*g)[site_vertex].type << "\n";
            //if ( (*g)[site_vertex].type == OUTER ) {
            //    std::cout << " OUTER face edge \n";
            //    return false;
            //}
            // we need to search for an adjacent linesite. 
            // (? can we have a situation where this fails?)
            HEEdge linetwin = find_adjacent_linesite(f);
            if (linetwin != HEEdge()) {
                HEEdge twin = (*g)[linetwin].twin;
                HEFace twin_face = (*g)[twin].face;
                if (linesite_ccw(twin_face))
                    return true;
            } else {
                return false;
            }
        } 
        return false;
    }
private:
    /// on the face f, find the adjacent linesite
    HEEdge find_adjacent_linesite(  HEFace f ) const {
        HEEdge current = (*g)[f].edge;
        HEEdge start = current;
        
        do {
            HEEdge twin = (*g)[current].twin;
            if (twin != HEEdge() ) {
                //std::cout << (*g)[ (*g).source(current) ].index << " - " << (*g)[ (*g).target(current) ].index;
                //std::cout << " tw: " << (*g)[ (*g).source(twin) ].index << " - " << (*g)[ (*g).target(twin) ].index;
                //std::cout << " t= " << (*g)[ twin ].type << "\n";
                
                HEFace twf = (*g)[twin].face;
                if ( (*g)[twf].site->isLine() ) {
                    //std::cout << "  returning: " << (*g)[ (*g).source(current) ].index << " - " << (*g)[ (*g).target(current) ].index << "\n";
                    return current;
                }
            } else {
                //std::cout << (*g)[ (*g).source(current) ].index << " - " << (*g)[ (*g).target(current) ].index;
                //std::cout << " t= " << (*g)[ current ].type << " has no twin!\n";
            }
            current = (*g)[current].next;
        } while(current!=start);
        
        return HEEdge();
    }
    /// return true if linesite was inserted in the direction indicated by _side
    bool linesite_ccw(  HEFace f ) const {
        HEEdge current = (*g)[f].edge;
        HEEdge start = current;
        do {
            if ( (_side && (*g)[current].type == LINESITE && (*g)[current].inserted_direction) ||
                  (!_side && (*g)[current].type == LINESITE && !(*g)[current].inserted_direction)  )
                return true;
                
            current = (*g)[current].next;
        } while(current!=start);
        return false;
    }
    /// CW / CCW flag
    bool _side;
};

} // end ovd namespace
// end file polygon_interior.hpp
