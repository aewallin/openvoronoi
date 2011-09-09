/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  
 *  This file is part of OpenVoronoi.
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


#include "vertex_positioner.hpp"
#include "voronoidiagram.hpp"

namespace ovd {


    Point VertexPositioner::position(const Point& p1, const Point& p2, const Point& p3) {
        Point pi(p1),pj(p2),pk(p3);
        if ( pi.isRight(pj,pk) ) 
            std::swap(pi,pj);
        assert( !pi.isRight(pj,pk) );
        // 2) point pk should have the largest angle. largest angle is opposite longest side.
        double longest_side = (pi - pj).norm();
        while (  ((pj - pk).norm() > longest_side) || (((pi - pk).norm() > longest_side)) ) { 
            std::swap(pi,pj); // cyclic rotation of points until pk is opposite the longest side pi-pj
            std::swap(pi,pk);  
            longest_side = (pi - pj).norm();
        }
        assert( !pi.isRight(pj,pk) );
        assert( (pi - pj).norm() >=  (pj - pk).norm() );
        assert( (pi - pj).norm() >=  (pk - pi).norm() );
        double J2 = (pi.y-pk.y)*( sq(pj.x-pk.x)+sq(pj.y-pk.y) )/2.0 - (pj.y-pk.y)*( sq(pi.x-pk.x)+sq(pi.y-pk.y) )/2.0;
        double J3 = (pi.x-pk.x)*( sq(pj.x-pk.x)+sq(pj.y-pk.y) )/2.0 - (pj.x-pk.x)*( sq(pi.x-pk.x)+sq(pi.y-pk.y) )/2.0;
        double J4 = (pi.x-pk.x)*(pj.y-pk.y) - (pj.x-pk.x)*(pi.y-pk.y);
        assert( J4 != 0.0 );
        return Point( -J2/J4 + pk.x, J3/J4 + pk.y );
    }

    // signature: edge, new_site 
    Point VertexPositioner::position(VoronoiDiagram* vd, HEEdge e, HEVertex v) {
        HEFace face = vd->g[e].face;     assert(  vd->g[face].status == INCIDENT);
        HEEdge twin = vd->g[e].twin;
        HEFace twin_face = vd->g[twin].face;      assert( vd->g[twin_face].status == INCIDENT);
        
        return position( vd->g[face].generator  , vd->g[twin_face].generator  , vd->g[v].position );
    }



}
