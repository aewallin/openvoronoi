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

#include "medial_axis_walk.hpp"

namespace ovd
{


/// \brief add the given edge to the current list of edges.
///
/// for line-edges we add only two endpoints
/// for parabolic edges we add many points
void MedialAxisWalk::append_edge(MedialChain& chain, HEEdge edge)  {
    MedialPointList point_list; // the endpoints of each edge
    HEVertex v1 = g.source( edge );
    HEVertex v2 = g.target( edge );
    // these edge-types are drawn as a single line from source to target.
    if ( (g[edge].type == LINELINE)  || (g[edge].type == PARA_LINELINE) ) {
        MedialPoint pt1( g[v1].position, g[v1].dist() );
        MedialPoint pt2( g[v2].position, g[v2].dist() );
        point_list.push_back(pt1);
        point_list.push_back(pt2);
    } else if ( (g[edge].type == PARABOLA) || (g[edge].type == LINE) ) { // these edge-types are drawn as polylines with _edge_points number of points
        double t_src = g[v1].dist();
        double t_trg = g[v2].dist();
        double t_min = std::min(t_src,t_trg);
        double t_max = std::max(t_src,t_trg);
        
        for (int n=0;n< _edge_points;n++) {
            double t(0);
            if (t_src<=t_trg) // increasing t-value
                t = t_min + ((t_max-t_min)/numeric::sq(_edge_points-1))*numeric::sq(n); // NOTE: quadratic t-dependece. More points at smaller t.
            else if (t_src>t_trg) { // decreasing t-value
                int m = _edge_points-1-n; // m goes from (N-1)...0   as n goes from 0...(N-1)
                t = t_min + ((t_max-t_min)/numeric::sq(_edge_points-1))*numeric::sq(m);
            }
            Point p = g[edge].point(t);
            MedialPoint pt( p, t );
            point_list.push_back(pt);
        }
    }
    chain.push_back( point_list );
}
    

}
