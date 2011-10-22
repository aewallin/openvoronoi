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

#include "voronoivertex.hpp"

namespace ovd {


/// point-point-point vertex positioner based on Sugihara & Iri paper
Point VertexPositioner::ppp_solver(const Point& p1, const Point& p2, const Point& p3) {
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


Point VertexPositioner::position(Site* s1, Site* s2, Site* s3) {
    int count;
    double solns[2][3];
    if ( s1->is_linear() && s2->is_linear() && s3->is_linear() )
        count =  lll_solver(s1,s2,s3);
    else
        count = solver(s1,s2,s3, solns);
    
    //std::cout << count << " solutions found by solver \n";
    for (int m=0;m<count;m++)
        std::cout << " new: " << m << " :  ( " << solns[m][0] << " , " << solns[m][1] << " , " << solns[m][2] << " )\n";
    
    Point pt = ppp_solver( s1->position(), s2->position(), s3->position() );
    std::cout << " old: " << pt << "\n";
    return pt;
    //ppp_solver( s1->position(), s2->position(), s3->position() );
}

int VertexPositioner::solver(Site* s1, Site* s2, Site* s3, double solns[][3] ) {
    double vectors[3][4]; // hold eqn data here. three equations with parameters (a,b,k,c) each
    // indexes and count of linear/quadratic eqns
    int linear[3],    linear_count = 0;
    int quadratic[3], quadratic_count = 0;
    
    double xk, yk, rk;
    int kk;
    
    // populate vectors
    vectors[0][0] = s1->eqp().a;
    vectors[0][1] = s1->eqp().b;
    vectors[0][2] = s1->eqp().k;
    vectors[0][3] = s1->eqp().c;
    if ( s1->is_linear() )
        linear[linear_count++] = 0;
    else {
        quadratic[quadratic_count++] = 0;
        xk = s1->x(); // point: (x,y, r=0, kk=1)  circle( x,y,r,k[i] )
        yk = s1->y();
        rk = s1->r();
        kk = s1->k();
    }
    vectors[1][0] = s2->eqp().a;
    vectors[1][1] = s2->eqp().b;
    vectors[1][2] = s2->eqp().k;
    vectors[1][3] = s2->eqp().c;
    if ( s2->is_linear() )
        linear[linear_count++] = 1;
    else {
        quadratic[quadratic_count++] = 1;
        xk = s2->x(); // point: (x,y, r=0, kk=1)  circle( x,y,r,k[i] )
        yk = s2->y();
        rk = s2->r();
        kk = s2->k();
    }
    
    vectors[2][0] = s3->eqp().a;
    vectors[2][1] = s3->eqp().b;
    vectors[2][2] = s3->eqp().k;
    vectors[2][3] = s3->eqp().c;
    if ( s3->is_linear() )
        linear[linear_count++] = 2;
    else {
        quadratic[quadratic_count++] = 2;
        xk = s3->x(); // point: (x,y, r=0, kk=1)  circle( x,y,r,k[i] )
        yk = s3->y();
        rk = s3->r();
        kk = s3->k();
    }
        
    assert( linear_count < 3 ); // ==3 should be caught above by lll_solve()
    assert( quadratic_count > 0); // we should have one or more quadratic
    
    //std::cout << " linear_count = " << linear_count << "\n";
    //std::cout << " quadratic_count = " << quadratic_count << "\n";
    
    // now subtract one quadratic from another to obtain a system
    // of one quadratic and two linear eqns
    int v0 = quadratic[0]; // the one we subtract from the others
    int j;
    // Subtract one quadratic equation from all the others,
    // making new linear equations.
    for(int i=1; i<quadratic_count; i++) {
        int vi = quadratic[i]; // the one to do subtraction on
        for(j=0; j<4; j++) // four parameters: a,b,k,c
            vectors[vi][j] -= vectors[v0][j];
        linear[linear_count++] = vi; // now we have a new linear eqn
    }
    assert(linear_count == 2); // At this point, we should have exactly two linear equations.
    
    /*
    std::cout <<" vectors:\n";
    for (int m=0;m<3;m++) {
        std::cout << vectors[m][0] << "\t" << vectors[m][1] << "\t" << vectors[m][2] << "\t" << vectors[m][3] << "\n";
    }*/
    
    // TODO:  pick the solution appraoch with the best numerical stability.
    
    // index shuffling determines if we solve:
    // x and y in terms of t
    // y and t in terms of x
    // t and x in terms of y
    //double solns[2][3]; 
    
    
    int scount = qqq_solver(vectors[linear[0]], vectors[linear[1]], 0, 1, 2, xk, yk, kk, rk, solns);
    if (scount < 0) { // negative scoung when discriminant is zero, so shuffle around coord-indexes:
        scount = qqq_solver(vectors[linear[0]], vectors[linear[1]], 2, 0, 1, xk, yk, kk, rk, solns);
        if (scount < 0) {
            scount = qqq_solver(vectors[linear[0]], vectors[linear[1]], 1, 2, 0, xk, yk, kk, rk, solns);
        }
    }
    //std::cout << " solver() found " << scount << " roots\n";
    return scount;
}

// l0 first linear eqn
// l1 second linear eqn
// xi,yi,ti  indexes to shuffle around
// solns = output solution triplets (x,y,t) or (u,v,t)
int VertexPositioner::qqq_solver( double l0[], double l1[], int xi, int yi, int ti, 
      double xk, double yk, int kk, double rk , double solns[][3] ) {
    
    double isolns[2][3];
    double aargs[3][2];
    double ai = l0[xi]; // first linear 
    double bi = l0[yi];
    double ki = l0[ti];
    double ci = l0[3];
    double aj = l1[xi]; // second linear
    double bj = l1[yi];
    double kj = l1[ti];
    double cj = l1[3];
    double d = ai*bj - aj*bi; // chop!
    if (d == 0) {// no solution can be found!
        //assert(0);
        
        /*
        std::cout << " qqq_solver() ZERO discriminant -> no solution! \n";
        for (int m=0;m<4;m++)
            std::cout << l0[m] << "   " << l1[m] << "\n";
        */
        return -1;
    }
    double a0 =  (bi*kj - bj*ki) / d;
    double a1 = -(ai*kj - aj*ki) / d;
    double b0 =  (bi*cj - bj*ci) / d;
    double b1 = -(ai*cj - aj*ci) / d;
    
    aargs[0][0] = 1.0;
    aargs[0][1] = -2*xk;
    aargs[1][0] = 1.0;
    aargs[1][1] = -2*yk;
    aargs[2][0] = -1.0;
    aargs[2][1] = -2*rk; // * kk!! (kk == sign of quadratic offset ?)
    
    // this solves for w, and returns either 0, 1, or 2 triplets of (u,v,t) in isolns
    int scount = qll_solve( aargs[xi][0], aargs[xi][1],
                            aargs[yi][0], aargs[yi][1],
                            aargs[ti][0], aargs[ti][1],
                            xk*xk + yk*yk - rk*rk,
                            a0, b0, 
                            a1, b1, isolns);
                            
    for (int i=0; i<scount; i++) {
        solns[i][xi] = isolns[i][0];       // u       x
        solns[i][yi] = isolns[i][1];       // v       y
        solns[i][ti] = isolns[i][2]; // t       t  chop!
    }
    //std::cout << " qqq_solver() found " << scount << " roots\n";
    return scount;
}

/// Solve a system of one quadratic equation, and two linear equations.
/// 
/// a0 u^2 + b0 u + c0 v^2 + d0 v + e0 w^2 + f0 w + g0 = 0
/// u = a1 w + b1
/// v = a2 w + b2
/// solve for w
/// then u and v, t can be calculated from w
int VertexPositioner::qll_solve( double a0, double b0, double c0, double d0, 
                      double e0, double f0, double g0, 
                      double a1, double b1, 
                      double a2, double b2, 
                      double soln[][3])
{
    //std::cout << " qll_solver() \n";
    double a, b, c;
    double roots[2], t;
    // TODO:  optimize using abs(a0) == abs(c0) == abs(d0) == 1
    a = (a0*(a1*a1) + c0*(a2*a2) + e0); // chop
    b = (2*a0*a1*b1 + 2*a2*b2*c0 + a1*b0 + a2*d0 + f0); // chop
    c = a0*(b1*b1) + c0*(b2*b2) + b0*b1 + b2*d0 + g0;
    int rcount = quadratic_roots(a, b, c, roots); // solves a*w^2 + b*w + c = 0
    if (rcount == 0)    // No roots, no solutions
        return 0;

    for (int i=0; i<rcount; i++) {
        t = roots[i];
        soln[i][0] = a1*t + b1; // u
        soln[i][1] = a2*t + b2; // v
        soln[i][2] = t;         // t
    }
    //std::cout << " qll_solver() found " << rcount << " roots\n";
    return rcount;
}

int VertexPositioner::quadratic_roots(double a, double b, double c, double roots[]) {
    if ((a == 0) and (b == 0))
        return 0;
    if (a == 0) {
        roots[0] = -c / b;
        return 1;
    }
    if (b == 0) {
        double sq = -c / a;
        if (sq > 0) {
            roots[0] = sqrt(sq);
            roots[1] = -roots[0];
            return 2;
        } else if (sq == 0) {
            roots[0] = 0;
            return 1;
        } else {
            return 0;
        }
    }
    double disc = (b*b - 4*a*c); // discriminant, chop!
    if (disc > 0) {
        double q;
        if (b > 0)
            q = (b + sqrt(disc)) / -2;
        else
            q = (b - sqrt(disc)) / -2;
        roots[0] = q / a;
        roots[1] = c / q; // (?)
        return 2;
    } else if (disc == 0) {
        roots[0] = (-b / (2*a));
        return 1;
    }
    return 0;
}

int VertexPositioner::lll_solver(Site* s1, Site* s2, Site* s3) {
    std::cout << " lll_solver() \n";
    assert(0); // NOT implemented yet!
    //return Point(0,0);
    return 0;
}


// signature: edge, new_site 
Point VertexPositioner::position(HEEdge e, Site* s) {
    HEFace face = vd->g[e].face;     assert(  vd->g[face].status == INCIDENT);
    HEEdge twin = vd->g[e].twin;
    HEFace twin_face = vd->g[twin].face;      assert( vd->g[twin_face].status == INCIDENT);
    
    // Point p = position( vd->g[face].generator  , vd->g[twin_face].generator  , vd->g[v].position );
    
    Point p = position( vd->g[face].site  , vd->g[twin_face].site  , s );
    
    check_far_circle(p);
    /*
    if ( !check_in_edge(e, p, v) ) {
        std::cout << " gen1= " << vd->g[face].generator << "\n";
        std::cout << " gen2= " << vd->g[twin_face].generator << "\n";
        std::cout << " gen3= " << vd->g[v].position << "\n";
    }*/
    
    check_on_edge(e, p);
    //check_dist(e, p, v);
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
    double d1 = (p - vd->g[face].site->position() ).norm_sq();
    double d2 = (p - vd->g[twin_face].site->position() ).norm_sq();  
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
    
    double d1 = (p - vd->g[face].site->position() ).norm_sq();
    double d2 = (p - vd->g[twin_face].site->position() ).norm_sq();  
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
