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
    Point VertexPositioner::position(HEEdge e, HEVertex v) {
        HEFace face = vd->g[e].face;     assert(  vd->g[face].status == INCIDENT);
        HEEdge twin = vd->g[e].twin;
        HEFace twin_face = vd->g[twin].face;      assert( vd->g[twin_face].status == INCIDENT);
        
        Point p = position( vd->g[face].generator  , vd->g[twin_face].generator  , vd->g[v].position );
        assert( check_far_circle(p) );
        return p;
    }


bool VertexPositioner::check_far_circle(const Point& p) {
    if (!(p.norm() < 18*vd->far_radius)) {
        std::cout << "WARNING check_far_circle() new vertex outside far_radius! \n";
        std::cout << p << " norm=" << p.norm() << " far_radius=" << vd->far_radius << "\n"; 
        return false;
    }
    return true;
}

/*    
void VoronoiDiagram::check_vertex_on_edge(HEVertex q, HEEdge e) {
    // sanity check on new vertex

    assert( g[q].position.norm() < 18*far_radius);
    
    HEVertex trg = g.target(e);
    HEVertex src = g.source(e);
    Point trgP = g[trg].position;
    Point srcP = g[src].position;
    Point newP = g[q].position;
    
    {
        double dtl_orig = g[q].position.xyDistanceToLine(srcP, trgP);
        //double dtl(dtl_orig);
        if (dtl_orig > 1e-3* ( trgP - srcP ).norm() ) {
            double t = ( g[q].position - srcP).dot( trgP - srcP ) / ( trgP - srcP ).dot( trgP - srcP ) ;
            //g[q].position = srcP + t*( trgP-srcP);
            //dtl = g[q].position.xyDistanceToLine(srcP, trgP);
            std::cout << "WARNING!! check_vertex_on_edge()  dtl= " << dtl_orig << " t= " << t  << " edgelength= " << ( trgP - srcP ).norm()   <<"\n";
        }
    }
    
    
    if (( trgP - srcP ).norm() <= 0 ) {
        std::cout << "WARNING check_vertex_on_edge() zero-length edge! \n";
        g[q].position = srcP;
    } else {
        assert( ( trgP - srcP ).norm() > 0.0 ); // edge has finite length
        assert( ( trgP - srcP ).dot( trgP - srcP ) > 0.0 ); // length squared
        double torig = ((newP - srcP).dot( trgP - srcP )) / ( trgP - srcP ).dot( trgP - srcP ) ;
        bool warn = false;
        double t(torig);
        if (torig < 0.0) { // clamp the t-parameter to [0,1]
            warn = true;
            t=0.0;
        } else if (torig> 1.0) {
            warn = true;
            t=1.0;
        }
        if ( warn ) {
            std::cout << "WARNING: check_vertex_on_edge() t_old= " << torig << " CORRECTED t_new= " << t << "\n";
            std::cout << "src= " << srcP << " new= " << newP << " trg= " << trgP << "\n";
            std::cout << " (src-trg).norm()= " << (srcP-trgP).norm() << "\n";
            g[q].position = srcP + t*( trgP-srcP);
            t = ( g[q].position - srcP).dot( trgP - srcP ) / ( trgP - srcP ).dot( trgP - srcP ) ;
            
        }
        // now we are clamped:
        assert( t >= 0.0 );
        assert( t <= 1.0 );        
        double dtl_orig = g[q].position.xyDistanceToLine(srcP, trgP);
        double dtl(dtl_orig);
        if (dtl_orig > 1e-3* ( trgP - srcP ).norm() ) {
            t = ( g[q].position - srcP).dot( trgP - srcP ) / ( trgP - srcP ).dot( trgP - srcP ) ;
            g[q].position = srcP + t*( trgP-srcP);
            dtl = g[q].position.xyDistanceToLine(srcP, trgP);
            std::cout << "WARNING check_vertex_on_edge()  old_dtl= " << dtl_orig << " new_dtl= " << dtl  << " edgelength= " << ( trgP - srcP ).norm()   <<"\n";
        }
        assert( dtl < 1e-3* ( trgP - srcP ).norm() );
    }
}*/

    
}
