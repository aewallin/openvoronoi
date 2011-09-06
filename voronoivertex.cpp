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

#define NDEBUG  // this turns off assertions
#include <cassert>

#include "voronoivertex.hpp"

namespace ovd {

inline double sq(double x) {return x*x;}

int VoronoiVertex::count = 0;

    VoronoiVertex::VoronoiVertex() {
        init();
        status = UNDECIDED;
        type = NORMAL;
    }
    /// construct vertex at position p with type t
    VoronoiVertex::VoronoiVertex( Point p, VoronoiVertexStatus st) {
        position=p;
        status=st;
        type = NORMAL;
        init();
    }
    
    VoronoiVertex::VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t) {
        position=p;
        status=st;
        type=t;
        init();
    }
    void VoronoiVertex::init() {
        index = count;
        count++;
        in_queue = false;
    }
    void VoronoiVertex::reset() {
        in_queue = false;
        status = UNDECIDED;
    }
    
    // set the generators and position the vertex
    void VoronoiVertex::set_generators(const Point& p1, const Point& p2, const Point& p3) {
        set_J(p1,p2,p3);
        set_position();
    }
    
    /// based on precalculated J2, J3, J4, _pk, calculate the H determinant for input Point pl
    /// Eq.(20) from Sugihara&Iri 1994
    /// H<0 means pl is inside the circle
    /// H==0 on the edge of the circle
    /// H>9 outside the circle
    double VoronoiVertex::detH(const Point& pl) const {
        return J2*(pl.x- pk.x) - J3*(pl.y-pk.y) + 0.5*J4*( sq(pl.x-pk.x) + sq(pl.y-pk.y) );
    }
    void VoronoiVertex::set_position() {
        position.x =  -J2/J4 + pk.x; 
        position.y =   J3/J4 + pk.y;
    }

    // this allows sorting points
    /*
    bool operator<(const VertexProps& other) const {
        return ( abs(this->H) < abs(other.H) );
    }*/
    /// set the J values
    /// pi, pj, pk define the three PointGenerators that position this vertex
    void VoronoiVertex::set_J(const Point& p1, const Point& p2, const Point& p3) { 
        // 1) i-j-k should come in CCW order
        Point pi(p1),pj(p2);
        pk=p3;
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
        
        J2 = detH_J2( pi, pj);
        J3 = detH_J3( pi, pj);
        J4 = detH_J4( pi, pj);
        assert( J4 != 0.0 ); // we need to divide by J4 later, so it better not be zero...
    }
    double VoronoiVertex::detH_J2(const Point& pi, const Point& pj) { // eq.21 of Sugihara&Iri 1994
        return (pi.y-pk.y)*( sq(pj.x-pk.x)+sq(pj.y-pk.y) )/2.0 - (pj.y-pk.y)*( sq(pi.x-pk.x)+sq(pi.y-pk.y) )/2.0;
    }
    double VoronoiVertex::detH_J3(const Point& pi, const Point& pj) { // eq.22 of Sugihara&Iri 1994
        return (pi.x-pk.x)*( sq(pj.x-pk.x)+sq(pj.y-pk.y) )/2.0 - (pj.x-pk.x)*( sq(pi.x-pk.x)+sq(pi.y-pk.y) )/2.0;
    }
    double VoronoiVertex::detH_J4(const Point& pi, const Point& pj) {
        return (pi.x-pk.x)*(pj.y-pk.y) - (pj.x-pk.x)*(pi.y-pk.y); // eq.23 of Sugihara&Iri 1994
    }

} // end ocl namespace
