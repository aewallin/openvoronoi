/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 * 
 *  Idea and code for point/line/arc voronoi-vertex positioning code by
 *  Andy Payne, andy "at" payne "dot" org, November, 2010
 *  see: http://www.payne.org/index.php/Calculating_Voronoi_Nodes
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

#include <algorithm> // std::erase()

#include "vertex_positioner.hpp"
#include "voronoidiagram.hpp"
#include "numeric.hpp"

using namespace ovd::numeric; // sq() chop()

namespace ovd {

// calculate the position of a new vertex on the given edge e
// the edge e holds information about which face it belongs to.
// each face holds information about which site created it
// so the three sites defining the position of the vertex are:
// - site to the left of HEEdge e
// - site to the right of HEEdge e
// - given new Site s
Point VertexPositioner::position(HEEdge e, Site* s) {
    edge = e;
    HEFace face = vd->g[e].face;     
    HEEdge twin = vd->g[e].twin;
    HEFace twin_face = vd->g[twin].face;      
    assert(  vd->g[face].status == INCIDENT);
    assert( vd->g[twin_face].status == INCIDENT);
    
    HEVertex src = vd->g.source(e);
    HEVertex trg = vd->g.target(e);
    double t_src = vd->g[src].dist();
    double t_trg = vd->g[trg].dist();
    t_min = std::min( t_src, t_trg );
    t_max = std::max( t_src, t_trg );

    std::cout << "  sites: " << vd->g[face].site->str() << "(k="<< vd->g[e].k<< ") ";
    std::cout << vd->g[twin_face].site->str() << "(k="<< vd->g[twin].k;
    std::cout << ") new= " << s->str() << "\n";
    std::cout << " t-vals t_min= " << t_min << " t_max= " << t_max << "\n";
        
    Solution sl = position( vd->g[face].site  , vd->g[e].k, vd->g[twin_face].site  , vd->g[twin].k, s );
    
    std::cout << " new vertex positioned at " << sl.p << " t=" << sl.t << " k3=" << sl.k3;
    std::cout << " err=" << edge_error(edge,sl) << "\n";
    assert( solution_on_edge(sl) );
    check_far_circle(sl.p);
    //check_on_edge(e, p);
    //check_dist(e, p, v);
    k3 = sl.k3;
    return sl.p;
}

// find vertex that is equidistant from s1, s2, s3
// should lie on the k1 side of s1, k2 side of s2
// we try both k3=-1 and k3=+1 for s3
Solution VertexPositioner::position(Site* s1, double k1, Site* s2, double k2, Site* s3) {
    assert( (k1==1) || (k1 == -1) );
    assert( (k2==1) || (k2 == -1) );
    std::vector<Solution> solutions;
    solver(s1,k1,s2,k2,s3,+1, solutions);
    if (!s3->isPoint()) // for points k3=+1 allways
        solver(s1,k1,s2,k2,s3,-1, solutions); // for lineSite or ArcSite we try k3=-1 also    

    // choose only t_min < t < t_max solutions 
    solutions.erase( std::remove_if(solutions.begin(),solutions.end(), t_filter(t_min,t_max) ), solutions.end() );
    // choose only in_region() solutions
    solutions.erase( std::remove_if(solutions.begin(),solutions.end(), in_region_filter(s3) ), solutions.end() );
    
    if ( solutions.size() == 1)
        return solutions[0];
    else if (solutions.size()>1) {
        // two or more points remain so we must further filter here!
        // filter further using edge_error
        double min_error=100;
        Solution min_solution(Point(0,0),0,0);
        std::cout << " edge_error filter: \n";
        BOOST_FOREACH(Solution s, solutions) {
            double err = edge_error(edge,s);
            std::cout << s.p << " k3=" << s.k3 << " t=" <<  s.t << " err=" << err << "\n";
            if ( err < min_error) {
                min_solution = s;
                min_error = err;
            }
        }
        if (min_error >= 1e-6) {
            std::cout << " sln=" << min_solution.p << " err=" << min_error << "\n";
            std::cout << " edge: " << vd->g[ vd->g.source(edge) ].position << " - " << vd->g[ vd->g.target(edge) ].position;
            std::cout << " edge-point(t="<<min_solution.t << ")= " << vd->g[edge].point(min_solution.t) << "\n";
        }
        //assert( min_error < 1e-6 );
        return min_solution;
    } 

    // either 0, or >= 2 solutions found. error.
    std::cout << " None, or too many solutions found! candidates are:\n";
    BOOST_FOREACH(Solution s, solutions ) {
        std::cout << s.p << " t=" << s.t << " k3=" << s.k3 << " tr=" << s3->in_region(s.p) << " e_err=" << edge_error(edge,s) <<"\n";
    }
    std::cout << " edge: " << vd->g[ vd->g.source(edge) ].position << " - " << vd->g[ vd->g.target(edge) ].position << "\n";
    assert(0);
    return Solution( Point(0,0), -1, 1 );
}

int VertexPositioner::solver(Site* s1, double k1, Site* s2, double k2, Site* s3, double kk3, std::vector<Solution>& solns) {
    qd_real vectors[3][4]; // hold eqn data here. three equations with four parameters (a,b,k,c) each
    // indexes and count of linear/quadratic eqns
    int linear[3],    linear_count = 0;
    int quadratic[3], quadratic_count = 0;
    
    std::vector< Eq<qd_real> > equations(3); // equation-parameters, in quad-precision
    equations[0] = s1->eqp(k1);
    equations[1] = s2->eqp(k2);
    equations[2] = s3->eqp(kk3);
    int lin_cnt = 0, quad_cnt=0;
    BOOST_FOREACH( Eq<qd_real>& eqn, equations ) {
        if ( eqn.q )
            quad_cnt++;
        else
            lin_cnt++;
    }
    
    qd_real xk, yk, rk, kk; // holds center, radius, direction of 'last' quadratic, if present
    
    // populate vectors
    vectors[0][0] = s1->eqp().a;
    vectors[0][1] = s1->eqp().b;
    vectors[0][2] = s1->eqp().k * k1;
    vectors[0][3] = s1->eqp().c;
    if ( s1->is_linear() )
        linear[linear_count++] = 0;
    else {
        quadratic[quadratic_count++] = 0;
        xk = s1->x(); // point: (x,y, r=0, kk=1)  circle( x,y,r,k[i] )
        yk = s1->y();
        rk = s1->r();
        kk = s1->k()*k1;
        if (s1->isPoint())
            kk=k1;
    }
    vectors[1][0] = s2->eqp().a;
    vectors[1][1] = s2->eqp().b;
    vectors[1][2] = s2->eqp().k * k2;
    vectors[1][3] = s2->eqp().c;
    if ( s2->is_linear() )
        linear[linear_count++] = 1;
    else {
        quadratic[quadratic_count++] = 1;
        xk = s2->x(); // point: (x,y, r=0, kk=1)  circle( x,y,r,k[i] )
        yk = s2->y();
        rk = s2->r();
        kk = s2->k()*k2;
        if (s2->isPoint())
            kk=k2;
    }
    
    vectors[2][0] = s3->eqp().a;
    vectors[2][1] = s3->eqp().b;
    vectors[2][2] = s3->eqp().k * kk3;
    vectors[2][3] = s3->eqp().c;
    if ( s3->is_linear() )
        linear[linear_count++] = 2;
    else {
        quadratic[quadratic_count++] = 2;
        xk = s3->x(); // point: (x,y, r=0, kk=1)  circle( x,y,r,k[i] )
        yk = s3->y();
        rk = s3->r();
        kk = s3->k()*kk3;
        if (s3->isPoint())
            kk=kk3;
    }
    
    assert( linear_count == lin_cnt );
    assert( quadratic_count == quad_cnt );
    
    if (linear_count == 3)
        return lll_solver(equations, kk3, solns); //kk3 just passes through, has no effect in kk3
    
    assert( linear_count < 3 ); // ==3 should be caught above by lll_solve()
    assert( quadratic_count > 0); // we should have one or more quadratic
    
    if (linear_count != 2 ) { // 3 caught above. 2 can skip this. 1 needs this code
        // now subtract one quadratic from another to obtain a system
        // of one quadratic and two linear eqns
        int v0 = quadratic[0]; // the one we subtract from the others
        // Subtract one quadratic equation from all the others,
        // making new linear equations.
        for(int i=1; i<quadratic_count; i++) {
            int vi = quadratic[i]; // the one to do subtraction on
            for(int j=0; j<4; j++) // four parameters: a,b,k,c
                vectors[vi][j] -= vectors[v0][j];
            linear[linear_count++] = vi; // now we have a new linear eqn
        }
        assert(linear_count == 2); // At this point, we should have exactly two linear equations.
    }
    
    // TODO:  pick the solution appraoch with the best numerical stability.

    // index shuffling determines if we solve:
    // x and y in terms of t
    // y and t in terms of x
    // t and x in terms of y
    int scount = qqq_solver(vectors[linear[0]], vectors[linear[1]], 0, 1, 2, xk, yk, kk, rk, kk3, solns);
    if (scount <= 0) { // negative scount when discriminant is zero, so shuffle around coord-indexes:
        scount = qqq_solver(vectors[linear[0]], vectors[linear[1]], 2, 0, 1, xk, yk, kk, rk, kk3, solns);
        if (scount <= 0) {
            scount = qqq_solver(vectors[linear[0]], vectors[linear[1]], 1, 2, 0, xk, yk, kk, rk, kk3, solns);
        }
    }
    //std::cout << " solver() found " << scount << " roots\n";
    return scount;
}

// l0 first linear eqn
// l1 second linear eqn
// xi,yi,ti  indexes to shuffle around
// xk, yk, kk, rk = params of one ('last') quadratic site (point or arc)
// solns = output solution triplets (x,y,t) or (u,v,t)
// returns number of solutions found
int VertexPositioner::qqq_solver( qd_real l0[], qd_real l1[], int xi, int yi, int ti, 
      qd_real xk, qd_real yk, qd_real kk, qd_real rk , qd_real kk3, std::vector<Solution>& solns) { 

    qd_real ai = l0[xi]; // first linear 
    qd_real bi = l0[yi];
    qd_real ki = l0[ti];
    qd_real ci = l0[3];
    qd_real aj = l1[xi]; // second linear
    qd_real bj = l1[yi];
    qd_real kj = l1[ti];
    qd_real cj = l1[3];

    qd_real d = chop( ai*bj - aj*bi ); // chop! (determinant for 2 linear eqns (?))
    if (d == 0) // no solution can be found!
        return -1;
    // these are the w-equations for qll_solve()
    // (2) u = a1 w + b1
    // (3) v = a2 w + b2
    qd_real a0 =  (bi*kj - bj*ki) / d;
    qd_real a1 = -(ai*kj - aj*ki) / d;
    qd_real b0 =  (bi*cj - bj*ci) / d;
    qd_real b1 = -(ai*cj - aj*ci) / d;
    // based on the 'last' quadratic of (s1,s2,s3)
    qd_real aargs[3][2];
    aargs[0][0] = 1.0;
    aargs[0][1] = -2*xk;   // quad.a
    aargs[1][0] = 1.0;
    aargs[1][1] = -2*yk;   // quad.b
    aargs[2][0] = -1.0;
    aargs[2][1] = -2*rk*kk; // quad.k*kk
    
    qd_real isolns[2][3];
    // this solves for w, and returns either 0, 1, or 2 triplets of (u,v,t) in isolns
    // NOTE: indexes of aargs shuffled depending on (xi,yi,ti) !
    int scount = qll_solve( aargs[xi][0], aargs[xi][1],
                            aargs[yi][0], aargs[yi][1],
                            aargs[ti][0], aargs[ti][1],
                            xk*xk + yk*yk - rk*rk, // quad.c
                            a0, b0, 
                            a1, b1, isolns);
    double tsolns[2][3];
    for (int i=0; i<scount; i++) {
        tsolns[i][xi] = to_double(isolns[i][0]);       // u       x
        tsolns[i][yi] = to_double(isolns[i][1]);       // v       y
        tsolns[i][ti] = to_double(isolns[i][2]);       // t       t  chop!
        solns.push_back( Solution( Point( tsolns[i][0], tsolns[i][1] ), tsolns[i][2], to_double(kk3) ) );
    }
    std::cout << " k3="<<kk3<<" qqq_solve found " << scount << " roots\n";
    return scount;
}

/// Solve a system of one quadratic equation, and two linear equations.
/// 
/// (1) a0 u^2 + b0 u + c0 v^2 + d0 v + e0 w^2 + f0 w + g0 = 0
/// (2) u = a1 w + b1
/// (3) v = a2 w + b2
/// solve (1) for w (can have 0, 1, or 2 roots)
/// then substitute into (2) and (3) to find (u, v, t)
int VertexPositioner::qll_solve( qd_real a0, qd_real b0, qd_real c0, qd_real d0, 
                      qd_real e0, qd_real f0, qd_real g0, 
                      qd_real a1, qd_real b1, 
                      qd_real a2, qd_real b2, 
                      qd_real soln[][3])
{
    // TODO:  optimize using abs(a0) == abs(c0) == abs(d0) == 1
    qd_real a = chop( (a0*(a1*a1) + c0*(a2*a2) + e0) ); 
    qd_real b = chop( (2*a0*a1*b1 + 2*a2*b2*c0 + a1*b0 + a2*d0 + f0) ); 
    qd_real c = a0*(b1*b1) + c0*(b2*b2) + b0*b1 + b2*d0 + g0;
    std::vector<qd_real> roots = quadratic_roots(a, b, c); // solves a*w^2 + b*w + c = 0
    if ( roots.empty() ) { // No roots, no solutions
        std::cout << " qll_solve no w roots. no solutions.\n";
        return 0;
    }
    for (unsigned int i=0; i<roots.size(); i++) {
        qd_real w = roots[i];
        soln[i][0] = a1*w + b1; // u
        soln[i][1] = a2*w + b2; // v
        soln[i][2] = w;         // t
    }
    return roots.size();
}

int VertexPositioner::lll_solver(std::vector< Eq<qd_real> >& eq, double kk3, std::vector<Solution>& slns ) {
    std::cout << " lll_solver() k3= " << k3 << "\n";
    unsigned int i = 0, j=1, k=2;
    qd_real d = chop( ( eq[i].a*eq[j].b - eq[j].a*eq[i].b)*eq[k].k + 
                      (-eq[i].a*eq[k].b + eq[k].a*eq[i].b)*eq[j].k +  
                      ( eq[j].a*eq[k].b - eq[k].a*eq[j].b)*eq[i].k   ); // determinant
    if (d != 0) {
        qd_real t = ( (-eq[i].a*eq[j].b + eq[j].a*eq[i].b)*eq[k].c + 
                      ( eq[i].a*eq[k].b - eq[k].a*eq[i].b)*eq[j].c + 
                      (-eq[j].a*eq[k].b + eq[k].a*eq[j].b)*eq[i].c   )/d;
        if (t >= 0) {
            qd_real sol_x = ( ( eq[i].b*eq[j].c - eq[j].b*eq[i].c)*eq[k].k + 
                              (-eq[i].b*eq[k].c + eq[k].b*eq[i].c)*eq[j].k + 
                              ( eq[j].b*eq[k].c - eq[k].b*eq[j].c)*eq[i].k   )/d;
            qd_real sol_y = ( (-eq[i].a*eq[j].c + eq[j].a*eq[i].c)*eq[k].k + 
                              ( eq[i].a*eq[k].c - eq[k].a*eq[i].c)*eq[j].k + 
                              (-eq[j].a*eq[k].c + eq[k].a*eq[j].c)*eq[i].k   )/d;
            slns.push_back( Solution( Point( to_double(sol_x), to_double(sol_y) ), to_double(t), kk3 ) ); // kk3 just passes through without any effect!?
            return 1;
        }
    }
    return 0; // no solution if determinant zero, or t-value negative
}

/// Old. not used anymore!
/// point-point-point vertex positioner based on Sugihara & Iri paper
Point VertexPositioner::ppp_solver(const Point& p1, const Point& p2, const Point& p3) {
    Point pi(p1),pj(p2),pk(p3);
    if ( pi.is_right(pj,pk) ) 
        std::swap(pi,pj);
    assert( !pi.is_right(pj,pk) );
    // 2) point pk should have the largest angle. largest angle is opposite longest side.
    double longest_side = (pi - pj).norm();
    while (  ((pj - pk).norm() > longest_side) || (((pi - pk).norm() > longest_side)) ) { 
        std::swap(pi,pj); // cyclic rotation of points until pk is opposite the longest side pi-pj
        std::swap(pi,pk);  
        longest_side = (pi - pj).norm();
    }
    assert( !pi.is_right(pj,pk) );
    assert( (pi - pj).norm() >=  (pj - pk).norm() );
    assert( (pi - pj).norm() >=  (pk - pi).norm() );
    double J2 = (pi.y-pk.y)*( sq(pj.x-pk.x)+sq(pj.y-pk.y) )/2.0 - (pj.y-pk.y)*( sq(pi.x-pk.x)+sq(pi.y-pk.y) )/2.0;
    double J3 = (pi.x-pk.x)*( sq(pj.x-pk.x)+sq(pj.y-pk.y) )/2.0 - (pj.x-pk.x)*( sq(pi.x-pk.x)+sq(pi.y-pk.y) )/2.0;
    double J4 = (pi.x-pk.x)*(pj.y-pk.y) - (pj.x-pk.x)*(pi.y-pk.y);
    assert( J4 != 0.0 );
    return Point( -J2/J4 + pk.x, J3/J4 + pk.y );
}


bool VertexPositioner::solution_on_edge(Solution& s) {
    double err = edge_error(edge,s);
    errstat.push_back(err);
    double limit = 5E-5;
    if ( err>=limit ) {
        std::cout << "solution_on_edge() ERROR err= " << err << "\n";
        std::cout << " edge: " << vd->g[ vd->g.source(edge) ].index << " - " << vd->g[ vd->g.target(edge) ].index << "\n";
    }
    return (err<limit);
}

double VertexPositioner::edge_error(HEEdge e, Solution& s) {
    Point ep = vd->g[e].point( s.t );
    return (ep-s.p).norm();
}

// new vertices should lie within the far_radius
bool VertexPositioner::check_far_circle(const Point& p) {
    if (!(p.norm() < 5*vd->far_radius)) {
        std::cout << "WARNING check_far_circle() new vertex outside far_radius! \n";
        std::cout << p << " norm=" << p.norm() << " far_radius=" << vd->far_radius << "\n"; 
        return false;
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

// all vertices should be of degree three, i.e. three adjacent faces/sites
// distance to the three adjacent sites should be equal
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
