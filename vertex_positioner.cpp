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

#include <boost/array.hpp>

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
    assert( check_dist(edge, sl, s) );
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

    std::cout << " solver() done \n";
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
            std::cout << "WARNING: EDGE ERROR TOO LARGE\n";
            std::cout << " s1 = " << s1->str2() << " k1= " << k1 << "\n";
            std::cout << " s2 = " << s2->str2() << " k2= " << k2 << "\n";
            std::cout << " s3 = " << s3->str2() << "\n";
            std::cout << " sln=" << min_solution.p << " err=" << min_error << "\n";
            std::cout << " edge: " << vd->g[ vd->g.source(edge) ].position << " - " << vd->g[ vd->g.target(edge) ].position;
            std::cout << " edge-point(t="<<min_solution.t << ")= " << vd->g[edge].point(min_solution.t) << "\n";
            //assert(0);
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
    std::vector< Eq<qd_real> > quads,lins; // equation-parameters, in quad-precision
    boost::array<Site*,3> sites = {{s1,s2,s3}};
    boost::array<double,3> kvals = {{k1,k2,kk3}};
    for (unsigned int i=0;i<3;i++) {
        Eq<qd_real> eqn = sites[i]->eqp_qd( kvals[i] );
        if (sites[i]->is_linear() ) // store site-equations in lins or quads
            lins.push_back( eqn ); 
        else
            quads.push_back( eqn );
    }
    
    if (lins.size() == 3) // all lines.
        return lll_solver( lins, kk3, solns); // kk3 just passes through, has no effect in lll_solver
    
    assert( !quads.empty() );
    
    if ( lins.size()==1 || lins.size() == 0 ) {
        assert( quads.size() == 3 || quads.size() == 2 );
        for (unsigned int i=1;i<quads.size();i++) {
            quads[i] = quads[i] - quads[0]; // subtract the first quad from the other one(s), to create new linear eqn(s)
            lins.push_back(quads[i]);
        }
    }
    assert( lins.size() == 2);  // At this point, we should have exactly two linear equations.
   
    // TODO:  pick the solution appraoch with the best numerical stability.
    // index shuffling determines if we solve:
    // x and y in terms of t
    // y and t in terms of x
    // t and x in terms of y
    /*
    int scount = qll_solver( lins, 0, 1, 2, quads[0], kk3, solns);
    if (scount <= 0) { // negative scount when discriminant is zero, so shuffle around coord-indexes:
        scount = qll_solver(lins, 2, 0, 1, quads[0], kk3, solns);
        if (scount <= 0) {
            scount = qll_solver(lins, 1, 2, 0, quads[0], kk3, solns);
        }
    }*/
    
    // call all three permutations
    qll_solver( lins, 0, 1, 2, quads[0], kk3, solns);
    qll_solver( lins, 2, 0, 1, quads[0], kk3, solns);
    qll_solver( lins, 1, 2, 0, quads[0], kk3, solns);
    
    //std::cout << " solver() found " << scount << " roots\n";
    return solns.size();
}

// l0 first linear eqn
// l1 second linear eqn
// xi,yi,ti  indexes to shuffle around
// xk, yk, kk, rk = params of one ('last') quadratic site (point or arc)
// solns = output solution triplets (x,y,t) or (u,v,t)
// returns number of solutions found
int VertexPositioner::qll_solver( const std::vector< Eq<qd_real> >& lins, int xi, int yi, int ti, 
      const Eq<qd_real>& quad, qd_real kk3, std::vector<Solution>& solns) { 
    assert( lins.size() == 2 );
    qd_real ai = lins[0][xi]; // first linear 
    qd_real bi = lins[0][yi];
    qd_real ki = lins[0][ti];
    qd_real ci = lins[0].c;
    
    qd_real aj = lins[1][xi]; // second linear
    qd_real bj = lins[1][yi];
    qd_real kj = lins[1][ti];
    qd_real cj = lins[1].c;
    
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
    aargs[0][1] = quad.a;
    aargs[1][0] = 1.0;
    aargs[1][1] = quad.b;
    aargs[2][0] = -1.0;
    aargs[2][1] = quad.k;
    
    qd_real isolns[2][3];
    // this solves for w, and returns either 0, 1, or 2 triplets of (u,v,t) in isolns
    // NOTE: indexes of aargs shuffled depending on (xi,yi,ti) !
    int scount = qll_solve( aargs[xi][0], aargs[xi][1],
                            aargs[yi][0], aargs[yi][1],
                            aargs[ti][0], aargs[ti][1],
                            quad.c, // xk*xk + yk*yk - rk*rk,
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
    std::cout << "qll_solver()\n";
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

// all three sites are lines
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
    double limit = 9E-5;
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

// all vertices should be of degree three, i.e. three adjacent faces/sites
// distance to the three adjacent sites should be equal
bool VertexPositioner::check_dist(HEEdge e, const Solution& sl, Site* s3) {
    //HEVertex trg = vd->g.target(e);
    //HEVertex src = vd->g.source(e);
    HEFace face = vd->g[e].face;     
    HEEdge tw_edge = vd->g[e].twin;
    HEFace twin_face = vd->g[tw_edge].face;      
    
    Site* s1 = vd->g[face].site;
    Site* s2 = vd->g[twin_face].site;
    
    double d1 = (sl.p - s1->apex_point(sl.p) ).norm();
    double d2 = (sl.p - s2->apex_point(sl.p) ).norm();  
    double d3 = (sl.p - s3->apex_point(sl.p) ).norm(); 
    
    double maxd = std::max( std::max( fabs(sl.t-d1),fabs(sl.t-d2)) , fabs(sl.t-d3));
    errstat.push_back(maxd);

        
    if ( !equal(d1,d2) || !equal(d1,d3) || !equal(d2,d3) ||
         !equal(sl.t,d1) || !equal(sl.t,d2) || !equal(sl.t,d3) ) {
        std::cout << "WARNING check_dist() ! \n";
        std::cout << "  sl.t= " << sl.t << "\n";
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
