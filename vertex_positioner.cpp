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
    check_far_circle(p);
    if ( !check_in_edge(e, p, v) ) {
        std::cout << " gen1= " << vd->g[face].generator << "\n";
        std::cout << " gen2= " << vd->g[twin_face].generator << "\n";
        std::cout << " gen3= " << vd->g[v].position << "\n";
    }
    
    check_on_edge(e, p);
    check_dist(e, p, v);
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

bool VertexPositioner::check_in_edge(HEEdge e, const Point& p, HEVertex v) {
    HEVertex trg = vd->g.target(e);
    HEVertex src = vd->g.source(e);
    Point trgP = vd->g[trg].position;
    Point srcP = vd->g[src].position;
    Point newP = p; //vd->g[q].position;
    if (( trgP - srcP ).norm() <= 0 ) {
        std::cout << "WARNING check_vertex_in_edge() zero-length edge! \n";
    } else {
        assert( ( trgP - srcP ).norm() > 0.0 ); // edge has finite length
        assert( ( trgP - srcP ).dot( trgP - srcP ) > 0.0 ); // length squared
        double t = ((newP - srcP).dot( trgP - srcP )) / ( trgP - srcP ).dot( trgP - srcP ) ;
        if ( t < 0.0 || t > 1.0  ) {
            std::cout << "WARNING: check_vertex_In_edge() t= " << t << "\n";
            std::cout << "    edge= " << vd->g[src].index << " - " << vd->g[trg].index << "\n";
            std::cout << "    src= " << vd->g[src].index << "  " << srcP << " error= " << error(e,srcP,v) << "\n";
            std::cout << "    new= " <<  newP << " error= " << error(e,p,v) << "\n";
            std::cout << "    trg= " << vd->g[trg].index << "  " << trgP << " error= " << error(e,trgP,v) << "\n";
            std::cout << "    (src-trg).norm()= " << (srcP-trgP).norm() << "\n";
            int Nmax = 100;
            double dt =1.0/(double)(Nmax-1);
            for (int n=0;n<Nmax;++n) {
                double mt = n*dt;
                Point mp = srcP + mt*(trgP-srcP);
                std::cout << "    mid= " <<  mp << " error= " << error(e, mp ,v) << "\n";
            }
            return false;
        }
    }
    return true;
}

double VertexPositioner::error(HEEdge e, const Point& p, HEVertex v) {
    //HEVertex trg = vd->g.target(e);
    //HEVertex src = vd->g.source(e);
    HEFace face = vd->g[e].face;     
    HEEdge twin = vd->g[e].twin;
    HEFace twin_face = vd->g[twin].face; 
    // distance from point p to all three generators
    double d1 = (p - vd->g[face].generator).norm_sq();
    double d2 = (p - vd->g[twin_face].generator).norm_sq();  
    double d3 = (p - vd->g[v].position).norm_sq(); 
    return sq(d1-d2)+sq(d1-d3)+sq(d2-d3);
}

bool VertexPositioner::check_on_edge(HEEdge e, const Point& p) {
    HEVertex trg = vd->g.target(e);
    HEVertex src = vd->g.source(e);
    Point trgP = vd->g[trg].position;
    Point srcP = vd->g[src].position;
    Point newP = p;
    double dtl = p.xyDistanceToLine(srcP, trgP);
    if (dtl > 1e-3* ( trgP - srcP ).norm() ) {
        std::cout << "WARNING!! check_vertex_on_edge()  dtl= " << dtl << "\n";
        std::cout << "    edge= " << vd->g[src].index << " - " << vd->g[trg].index << "\n";
        std::cout << "    (src-trg).norm()= " << (srcP-trgP).norm() << "\n";
        return false;
    }
    return true;
}

// distance to adjacent sites should be equal
bool VertexPositioner::check_dist(HEEdge e, const Point& p, HEVertex v) {
    HEVertex trg = vd->g.target(e);
    HEVertex src = vd->g.source(e);
    HEFace face = vd->g[e].face;     
    HEEdge twin = vd->g[e].twin;
    HEFace twin_face = vd->g[twin].face;      
    
    double d1 = (p - vd->g[face].generator).norm_sq();
    double d2 = (p - vd->g[twin_face].generator).norm_sq();  
    double d3 = (p - vd->g[v].position).norm_sq(); 
        
    if ( !equal(d1,d2) || !equal(d1,d3) || !equal(d2,d3) ) {
        std::cout << "WARNING check_dist() ! \n";
        std::cout << "  src.dist= " << vd->g[src].dist() << "\n";
        std::cout << "  trg.dist= " << vd->g[trg].dist() << "\n";
    
        std::cout << "  d1= " << d1 << "\n"; 
        std::cout << "  d2= " << d2 << "\n";
        std::cout << "  d3= " << d3 << "\n";
        return false;
    }
    return true;
}

bool VertexPositioner::equal(double d1, double d2) {
    bool tol = 1e-3;
    if ( fabs(d1-d2) > tol*std::max(d1,d2) )
        return false;
    return true;
}
    
    
} // end namespace
