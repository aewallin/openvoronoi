/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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
#include "common/numeric.hpp"

#include "solvers/solver_ppp.hpp"
#include "solvers/solver_lll.hpp"
#include "solvers/solver_qll.hpp"

using namespace ovd::numeric; // sq() chop()

namespace ovd {

VertexPositioner::VertexPositioner(VoronoiDiagram* vodi): vd(vodi) {
    ppp_solver = new PPPSolver<qd_real>();
    //ppp_solver = new PPPSolver<double>(); // faster, but inaccurate
    lll_solver = new LLLSolver();
    qll_solver = new QLLSolver();
    errstat.clear();
}

VertexPositioner::~VertexPositioner() {
    delete ppp_solver;
    delete lll_solver;
    delete qll_solver;
}

// calculate the position of a new vertex on the given edge e
// the edge e holds information about which face it belongs to.
// each face holds information about which site created it
// so the three sites defining the position of the vertex are:
// - site to the left of HEEdge e
// - site to the right of HEEdge e
// - given new Site s
Solution VertexPositioner::position(HEEdge e, Site* s) {
    edge = e;
    HEFace face = vd->g[e].face;     
    HEEdge twin = vd->g[e].twin;
    HEFace twin_face = vd->g[twin].face;      
    //assert(  vd->g[face].status == INCIDENT);
    //assert( vd->g[twin_face].status == INCIDENT);
    
    HEVertex src = vd->g.source(e);
    HEVertex trg = vd->g.target(e);
    double t_src = vd->g[src].dist();
    double t_trg = vd->g[trg].dist();
    t_min = std::min( t_src, t_trg );
    t_max = std::max( t_src, t_trg );
/*
    std::cout << "  sites: " << vd->g[face].site->str() << "(k="<< vd->g[e].k<< ") ";
    std::cout << vd->g[twin_face].site->str() << "(k="<< vd->g[twin].k;
    std::cout << ") new= " << s->str() << "\n";
    std::cout << " t-vals t_min= " << t_min << " t_max= " << t_max << "\n";
  */      
    Solution sl = position( vd->g[face].site  , vd->g[e].k, vd->g[twin_face].site  , vd->g[twin].k, s );
    /*
    std::cout << " new vertex positioned at " << sl.p << " t=" << sl.t << " k3=" << sl.k3;
    std::cout << " err=" << edge_error(edge,sl) << "\n";
    */
    assert( solution_on_edge(sl) );
    assert( check_far_circle(sl) );
    assert( check_dist(edge, sl, s) );
    errstat.push_back( dist_error(edge, sl, s) );
    if ( dist_error(edge, sl, s) > 1e-9 ) {
        std::cout << " VertexPositioner::position() ERROR dist_error = " << dist_error(edge,  sl, s) << "\n";
    }
    
    return sl;
}

// find vertex that is equidistant from s1, s2, s3
// should lie on the k1 side of s1, k2 side of s2
// we try both k3=-1 and k3=+1 for s3
Solution VertexPositioner::position(Site* s1, double k1, Site* s2, double k2, Site* s3) {
    assert( (k1==1) || (k1 == -1) );
    assert( (k2==1) || (k2 == -1) );
    std::vector<Solution> solutions;
    solver_dispatch(s1,k1,s2,k2,s3,+1, solutions);
    if (!s3->isPoint()) // for points k3=+1 allways
        solver_dispatch(s1,k1,s2,k2,s3,-1, solutions); // for lineSite or ArcSite we try k3=-1 also    

    if ( solutions.size() == 1 && (t_min<=solutions[0].t) && (t_max>=solutions[0].t) && (s3->in_region( solutions[0].p)) )
        return solutions[0];
    
    // choose only in_region() solutions
    solutions.erase( std::remove_if(solutions.begin(),solutions.end(), in_region_filter(s3) ), solutions.end() );
    if (solutions.empty() ) 
        std::cout << "WARNING in_region_filter() results in empty solution set!!\n";

    // choose only t_min < t < t_max solutions 
    solutions.erase( std::remove_if(solutions.begin(),solutions.end(), t_filter(t_min,t_max) ), solutions.end() );
    if (solutions.empty() ) 
        std::cout << "WARNING t_filter() results in empty solution set!!\n";

    
    
    if ( solutions.size() == 1)
        return solutions[0];
    else if (solutions.size()>1) {
        // two or more points remain so we must further filter here!
        // filter further using edge_error
        double min_error=100;
        Solution min_solution(Point(0,0),0,0);
        //std::cout << " edge_error filter: \n";
        BOOST_FOREACH(Solution s, solutions) {
            double err = vd->g[edge].error(s);
            //std::cout << s.p << " k3=" << s.k3 << " t=" <<  s.t << " err=" << err << "\n";
            if ( err < min_error) {
                min_solution = s;
                min_error = err;
            }
        }
        if (min_error >= 1e-6) {
            /*
            std::cout << "WARNING: EDGE ERROR TOO LARGE\n";
            std::cout << " s1 = " << s1->str2() << " k1= " << k1 << "\n";
            std::cout << " s2 = " << s2->str2() << " k2= " << k2 << "\n";
            std::cout << " s3 = " << s3->str2() << "\n";
            std::cout << " sln=" << min_solution.p << " err=" << min_error << "\n";
            std::cout << " edge: " << vd->g[ vd->g.source(edge) ].position << " - " << vd->g[ vd->g.target(edge) ].position;
            std::cout << " edge-point(t="<<min_solution.t << ")= " << vd->g[edge].point(min_solution.t) << "\n";
            */
            //assert(0);
        }
        //assert( min_error < 1e-6 );
        return min_solution;
    } 

    // either 0, or >= 2 solutions found. error.
    // std::cout << " None, or too many solutions found! solutions.size()=" << solutions.size() << "\n";
    std::cout << " solution edge: " << vd->g[ vd->g.source(edge) ].position << "[" << vd->g[ vd->g.source(edge) ].type << "](t=" << vd->g[ vd->g.source(edge) ].dist() << ")";
    std::cout << " - " << vd->g[ vd->g.target(edge) ].position << "[" << vd->g[ vd->g.target(edge) ].type << "](t=" << vd->g[ vd->g.target(edge) ].dist() << ") \n";
    std::cout << " solution edge: " << vd->g[ vd->g.source(edge) ].index << "[" << vd->g[ vd->g.source(edge) ].type<<"]{" << vd->g[ vd->g.source(edge) ].status<<"}";
    std::cout << " -[" << vd->g[edge].type << "]- ";
    std::cout << vd->g[ vd->g.target(edge) ].index << "[" << vd->g[ vd->g.target(edge) ].type << "]{" << vd->g[ vd->g.target(edge) ].status<<"}\n";
    //std::cout << " t-vals t_min= " << t_min << " t_max= " << t_max << "\n";
    //std::cout << "  sites: " << s1->str() << "(k="<< k1<< ") " << s2->str() << "(k="<< k2 << ") new= " << s3->str() << "\n";
    std::cout << "s1= " << s1->str2() << "(k=" << k1<< ")\n";
    std::cout << "s2= " << s2->str2() << "(k=" << k2<< ")\n";
    std::cout << "s3= " << s3->str2() << "\n";

    // run the solver(s) one more time in order to print out un-filtered solution points for debugging
    std::vector<Solution> solutions2;
    solver_dispatch(s1,k1,s2,k2,s3,+1, solutions2);
    if (!s3->isPoint()) // for points k3=+1 allways
        solver_dispatch(s1,k1,s2,k2,s3,-1, solutions2); // for lineSite or ArcSite we try k3=-1 also    

    std::cout << " The failing solutions are: \n";
    BOOST_FOREACH(Solution s, solutions2 ) {
        std::cout << s.p << " t=" << s.t << " k3=" << s.k3  << " e_err=" << vd->g[edge].error(s) <<"\n";
        std::cout << " min<t<max=" << ((s.t>=t_min) && (s.t<=t_max));
        std::cout << " s3.in_region=" << s3->in_region(s.p);
        std::cout <<  " region-t=" << s3->in_region_t(s.p) << "\n";
        std::cout <<  " t - t_min= " << s.t - t_min << "\n";
        std::cout <<  " t_max - t= " << t_max - s.t << "\n";
        //std::cout << std::scientific;
    }

    assert(0);
    // try a desperate solution
    double t_mid = 0.5*(t_min+t_max);
    Point p_mid = vd->g[edge].point(t_mid);
    Solution desp( p_mid, t_mid, 1 );
    std::cout << " Returning desperate solution: \n";
    std::cout << desp.p << " t=" << desp.t << " k3=" << desp.k3  << " e_err=" << vd->g[edge].error(desp) <<"\n";
    return desp;
}

int VertexPositioner::solver_dispatch(Site* s1, double k1, Site* s2, double k2, Site* s3, double k3, std::vector<Solution>& solns) {    
    if ( s1->isLine() && s2->isLine() && s3->isLine() ) 
        return lll_solver->solve( s1,k1,s2,k2,s3,k3, solns ); // all lines.
    else if ( s1->isPoint() && s2->isPoint() && s3->isPoint() )
        return ppp_solver->solve( s1,s2,s3, solns ); // all points, no need to specify k1,k2,k3, they are all +1
    else
        return qll_solver->solve( s1,k1,s2,k2,s3,k3, solns ); // general case solver
 
}

bool VertexPositioner::solution_on_edge(Solution& s) {
    double err = vd->g[edge].error(s);
    double limit = 9E-4;
    if ( err>=limit ) {
        std::cout << "solution_on_edge() ERROR err= " << err << "\n";
        std::cout << " solution edge: " << vd->g[ vd->g.source(edge) ].index << "[" << vd->g[ vd->g.source(edge) ].type<<"]{" << vd->g[ vd->g.source(edge) ].status<<"}";
        std::cout << " -[" << vd->g[edge].type << "]- ";
        std::cout << vd->g[ vd->g.target(edge) ].index << "[" << vd->g[ vd->g.target(edge) ].type << "]{" << vd->g[ vd->g.target(edge) ].status<<"}\n";


        std::cout << " edge: " << vd->g[ vd->g.source(edge) ].index << "(t=" << vd->g[ vd->g.source(edge) ].dist() << ")"; 
        std::cout << " - " << vd->g[ vd->g.target(edge) ].index << "(t=" << vd->g[ vd->g.target(edge) ].dist() << ")\n";
        std::cout << " edge: " << vd->g[ vd->g.source(edge) ].position << " - " << vd->g[ vd->g.target(edge) ].position << "\n";
        std::cout << " solution: " << s.p << " t=" << s.t << "\n";
    }
    return (err<limit);
}

// calculate the distance from the solution-point to the corresponding point on the edge.
/*
double VertexPositioner::edge_error(HEEdge e, Solution& s) {
    Point ep = vd->g[e].point( s.t, s );
    return (ep-s.p).norm();
}*/

// new vertices should lie within the far_radius
bool VertexPositioner::check_far_circle(Solution& s) {
    if (!(s.p.norm() < 18*vd->far_radius)) {
        std::cout << "WARNING check_far_circle() new vertex outside far_radius! \n";
        std::cout << s.p << " norm=" << s.p.norm() << " far_radius=" << vd->far_radius << "\n"; 
        return false;
    }
    return true;
}

// all vertices should be of degree three, i.e. three adjacent faces/sites
// distance to the three adjacent sites should be equal
bool VertexPositioner::check_dist(HEEdge e, const Solution& sl, Site* s3) {
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
        std::cout << " solution edge: " << vd->g[ vd->g.source(edge) ].index << "[" << vd->g[ vd->g.source(edge) ].type<<"]{" << vd->g[ vd->g.source(edge) ].status<<"}";
        std::cout << " -[" << vd->g[edge].type << "]- ";
        std::cout << vd->g[ vd->g.target(edge) ].index << "[" << vd->g[ vd->g.target(edge) ].type << "]{" << vd->g[ vd->g.target(edge) ].status<<"}\n";
    
        return false;
    }
    return true;
}

// new vertices should be equidistant to the three adjacent sites that define the vertex
// we here calculate the distances d1, d2, d3 from the Solution to the three sites s1, s2, s3
// and return the max deviation from the solution t-value.
// this works as a sanity check for the solver.
// a high error value here is also an indication of numerical instability in the solver
double VertexPositioner::dist_error(HEEdge e, const Solution& sl, Site* s3) {
    HEFace face = vd->g[e].face;     
    HEEdge tw_edge = vd->g[e].twin;
    HEFace twin_face = vd->g[tw_edge].face;      
    
    Site* s1 = vd->g[face].site;
    Site* s2 = vd->g[twin_face].site;
    
    double d1 = (sl.p - s1->apex_point(sl.p) ).norm();
    double d2 = (sl.p - s2->apex_point(sl.p) ).norm();  
    double d3 = (sl.p - s3->apex_point(sl.p) ).norm(); 
    
    return std::max( std::max( fabs(sl.t-d1),fabs(sl.t-d2)) , fabs(sl.t-d3));

}

bool VertexPositioner::equal(double d1, double d2) {
    bool tol = 1e-3;
    if ( fabs(d1-d2) < 1e-15 )
        return true;
    if ( fabs(d1-d2) > tol*std::max(d1,d2) )
        return false;
    return true;
}
    
    
} // end namespace
