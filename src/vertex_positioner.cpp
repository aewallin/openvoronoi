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
#include <boost/math/tools/minima.hpp> // brent_find_minima

#include "vertex_positioner.hpp"
#include "voronoidiagram.hpp"
#include "common/numeric.hpp"

#include "solvers/solver_ppp.hpp"
#include "solvers/solver_lll.hpp"
#include "solvers/solver_lll_para.hpp"

#include "solvers/solver_qll.hpp"
#include "solvers/solver_sep.hpp"
#include "solvers/solver_alt_sep.hpp"

using namespace ovd::numeric; // sq() chop()

namespace ovd {

/// create positioner, set graph.
VertexPositioner::VertexPositioner(HEGraph& gi): g(gi) {
    //ppp_solver = new solvers::PPPSolver<double>(); // faster, but inaccurate
    ppp_solver =      new solvers::PPPSolver<qd_real>(); // slower, more accurate
    lll_solver =      new solvers::LLLSolver();
    qll_solver =      new solvers::QLLSolver();
    sep_solver =      new solvers::SEPSolver();
    alt_sep_solver =  new solvers::ALTSEPSolver();
    lll_para_solver = new solvers::LLLPARASolver();
    silent = false;
    solver_debug(false);
    errstat.clear();
}

/// delete all solvers
VertexPositioner::~VertexPositioner() {
    //std::cout << "~VertexPositioner()..";
    delete ppp_solver;
    delete lll_solver;
    delete qll_solver;
    delete sep_solver;
    delete alt_sep_solver;
    delete lll_para_solver;

    errstat.clear();
    //std::cout << "DONE.\n";
}

/// \brief position a new vertex on given HEEdge \a e when inserting the new Site \a s3
///
/// calculate the position of a new voronoi-vertex lying on the given edge.
/// The new vertex is equidistant to the two sites that defined the edge
/// and to the new site. 
// the edge e holds information about which face it belongs to.
// each face holds information about which site created it
// so the three sites defining the position of the vertex are:
// - site to the left of HEEdge e
// - site to the right of HEEdge e
// - given new Site s
solvers::Solution VertexPositioner::position(HEEdge e, Site* s3) {
    edge = e;
    HEFace face = g[e].face;     
    HEEdge twin = g[e].twin;
    HEFace twin_face = g[twin].face;

    HEVertex src = g.source(e);
    HEVertex trg = g.target(e);
    double t_src = g[src].dist();
    double t_trg = g[trg].dist();
    t_min = std::min( t_src, t_trg ); // the solution we seek must have t_min<t<t_max
    t_max = std::max( t_src, t_trg );

    Site* s1 =  g[face].site;
    Site* s2 = g[twin_face].site;

    solvers::Solution sl = position(  s1 , g[e].k, s2, g[twin].k, s3 );

    assert( solution_on_edge(sl) );
    //assert( check_far_circle(sl) );
    assert( check_dist(edge, sl, s3) );
    
    // error logging (FIXME: make optional, for max performance?)
    #ifndef NDEBUG
    {
        errstat.push_back( dist_error(edge, sl, s3) );
        if ( dist_error(edge, sl, s3) > 1e-6 ) {
            // 2012-02-04: 1e-9 passes 79/79 tests
            //             1e-10 passes 79/79
            //             1e-12 passes 79/79
            //             1e-13  17 FAILED out of 79
            //             1e-14  38 FAILED out of 79
            std::cout << " VertexPositioner::position() WARNING; large dist_error = " << dist_error(edge,  sl, s3) << "\n";
            double s1_dist = (sl.p - s1->apex_point(sl.p)).norm();
            double s2_dist = (sl.p - s2->apex_point(sl.p)).norm();
            double s3_dist = (sl.p - s3->apex_point(sl.p)).norm();
            std::cout << " s1 dist = " << s1_dist << "\n";
            std::cout << " s2 dist = " << s2_dist << "\n";
            std::cout << " s3 dist = " << s3_dist << "\n";
            std::cout << " t       = " << sl.t << "\n";
            exit(-1);
            //return fabs(t-s3_dist);
        }
    }
    #endif
    
    return sl;
}

/// position new vertex
// find vertex that is equidistant from s1, s2, s3
// should lie on the k1 side of s1, k2 side of s2
// we try both k3=-1 and k3=+1 for s3
solvers::Solution VertexPositioner::position(Site* s1, double k1, Site* s2, double k2, Site* s3) {
    assert( (k1==1) || (k1 == -1) );
    assert( (k2==1) || (k2 == -1) );
    std::vector<solvers::Solution> solutions;
        
    solver_dispatch(s1,k1,s2,k2,s3,+1, solutions); // a single k3=+1 call for s3->isPoint()
    
    if (!s3->isPoint()) 
        solver_dispatch(s1,k1,s2,k2,s3,-1, solutions); // for lineSite or ArcSite we try k3=-1 also    
    
    if ( solutions.size() == 1 && (t_min<=solutions[0].t) && (t_max>=solutions[0].t) && (s3->in_region( solutions[0].p)) )
        return solutions[0];
            
    if (solutions.empty() && !silent ) 
        std::cout << "WARNING empty solution set!!\n";
    
    // choose only in_region() solutions
    solutions.erase( std::remove_if(solutions.begin(),solutions.end(), in_region_filter(s3) ), solutions.end() );
    if (solutions.empty() && !silent ) 
        std::cout << "WARNING in_region_filter() results in empty solution set!!\n";
    
    
    // choose only t_min < t < t_max solutions 
    solutions.erase( std::remove_if(solutions.begin(),solutions.end(), t_filter(t_min,t_max) ), solutions.end() );
    if (solutions.empty() && !silent ) 
        std::cout << "WARNING t_filter() results in empty solution set!!\n";

    if ( solutions.size() == 1) // if only one solution is found, return that.
        return solutions[0];
    
    if (solutions.size()>1) {
        std::vector<solvers::Solution> equidistant_solutions;
        // If s3 is a linesite or arcsite, we look for a solution in both halfspaces delimited by s3
        // Usually, the t_min / t_max should take care of filtering out the invalid solution, but in some
        // cases, t remains nearly constant (e.g. if s1 and s2 are nearly parallel) and the out-of-region solution slips through. 
        // Therefore, we remove any solutions which deviate from the equidistance constraint by more than 1%.
        for (std::size_t i = 0; i < solutions.size(); i++) {
            const solvers::Solution& s = solutions[i];
            double d1 = (s.p - s1->apex_point(s.p)).norm();
            double d2 = (s.p - s2->apex_point(s.p)).norm();
            double d3 = (s.p - s3->apex_point(s.p)).norm();
            double err = std::max(std::abs(d1-d2), std::max(std::abs(d2-d3), std::abs(d3-d1)));
            double mindist = std::min(d1, std::min(d2, d3));
            if (err/mindist > 0.01) {
                std::cout << "Solution "<<i<<" violates equidistance constraint. Distances of solution were: "<<sqrt(d1)<<", "<<sqrt(d2)<<", "<<sqrt(d3)<<std::endl;
            }
            else {
                equidistant_solutions.push_back(s);
            }
        }
        solutions = equidistant_solutions;
    }
    if ( solutions.size() == 1) // if only one solution is found, return that.
        return solutions[0];

    if (solutions.size()>1) {
        // two or more points remain so we must further filter here!
        // filter further using edge_error
        double min_error=100;
        solvers::Solution min_solution(Point(0,0),0,0);
        //std::cout << " edge_error filter: \n";
        BOOST_FOREACH(solvers::Solution s, solutions) {
            double err = edge_error(s); //g[edge].error(s);
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
            std::cout << " edge: " << g[ g.source(edge) ].position << " - " << g[ g.target(edge) ].position;
            std::cout << " edge-point(t="<<min_solution.t << ")= " << g[edge].point(min_solution.t) << "\n";
            */
            //assert(0);
        }
        //assert( min_error < 1e-6 );
        return min_solution;
    }
    

    // either 0, or >= 2 solutions found. This is an error.
    // std::cout << " None, or too many solutions found! solutions.size()=" << solutions.size() << "\n";
     
    if ( !silent) {
        std::cout << " solution edge: " << g[ g.source(edge) ].position << "[" << g[ g.source(edge) ].type << "](t=" << g[ g.source(edge) ].dist() << ")";
        std::cout << " - " << g[ g.target(edge) ].position << "[" << g[ g.target(edge) ].type << "](t=" << g[ g.target(edge) ].dist() << ") \n";
        std::cout << " solution edge: " << g[ g.source(edge) ].index << "[" << g[ g.source(edge) ].type<<"]{" << g[ g.source(edge) ].status<<"}";
        std::cout << " -[" << g[edge].type << "]- ";
        std::cout << g[ g.target(edge) ].index << "[" << g[ g.target(edge) ].type << "]{" << g[ g.target(edge) ].status<<"}\n";
        //std::cout << " t-vals t_min= " << t_min << " t_max= " << t_max << "\n";
        //std::cout << "  sites: " << s1->str() << "(k="<< k1<< ") " << s2->str() << "(k="<< k2 << ") new= " << s3->str() << "\n";
        std::cout << " s1= " << s1->str2() << "(k=" << k1<< ")\n";
        std::cout << " s2= " << s2->str2() << "(k=" << k2<< ")\n";
        std::cout << " s3= " << s3->str2() << "\n";
        
        std::cout << "Running solvers again: \n";
    }
    solver_debug(true);
    // run the solver(s) one more time in order to print out un-filtered solution points for debugging
    std::vector<solvers::Solution> solutions2;
    solver_dispatch(s1,k1,s2,k2,s3,+1, solutions2);
    if (!s3->isPoint()) // for points k3=+1 allways
        solver_dispatch(s1,k1,s2,k2,s3,-1, solutions2); // for lineSite or ArcSite we try k3=-1 also    
    solver_debug(false);
    
    if ( !silent) { 
        if ( !solutions2.empty() ) {
            std::cout << "The failing " << solutions2.size() << " solutions are: \n";
            BOOST_FOREACH(solvers::Solution s, solutions2 ) {
                std::cout << s.p << " t=" << s.t << " k3=" << s.k3  << " e_err=" << edge_error(s) <<"\n";
                std::cout << " min<t<max=" << ((s.t>=t_min) && (s.t<=t_max));
                std::cout << " s3.in_region=" << s3->in_region(s.p);
                std::cout <<  " region-t=" << s3->in_region_t(s.p) << "\n";
                std::cout <<  " t - t_min= " << s.t - t_min << "\n";
                std::cout <<  " t_max - t= " << t_max - s.t << "\n";
                std::cout <<  " edge type : " << g[edge].type << "\n"; //std::scientific;
            }   
        } else {
            std::cout << "No solutions found by solvers!\n";
        }
    }

    //assert(0); // in Debug mode, stop here.
    
    solvers::Solution desp = desperate_solution(s3);  // ( p_mid, t_mid, desp_k3 ); 
    
    VertexError s1_err_functor(g, edge, s1);
    VertexError s2_err_functor(g, edge, s2);
    VertexError s3_err_functor(g, edge, s3);
    
    if ( !silent) { 
        std::cout << "WARNING: Returning desperate solution: \n";
        std::cout << desp.p << " t=" << desp.t << " k3=" << desp.k3  << " e_err=" << edge_error(desp) <<"\n";
        std::cout << "     s1_err= " << s1_err_functor(desp.t) << "\n";
        std::cout << "     s2_err= " << s2_err_functor(desp.t) << "\n";
        std::cout << "     s3_err= " << s3_err_functor(desp.t) << "\n";
    }
    //exit(-1);
    return desp;
}

/// search numerically for a desperate solution along the solution-edge
solvers::Solution VertexPositioner::desperate_solution(Site* s3) {
    VertexError err_functor(g, edge, s3);
    //HEFace face = g[edge].face;     
    //HEEdge twin = g[edge].twin;
    //HEFace twin_face = g[twin].face;
    //Site* s1 =  g[face].site;
    //Site* s2 = g[twin_face].site;
    HEVertex src = g.source(edge);
    HEVertex trg = g.target(edge);
    Point src_p = g[src].position;
    Point trg_p = g[trg].position;
    
    if (!silent) {
        std::cout << "VertexPositioner::desperate_solution() \n";
        std::cout << " edge: " << src_p << " - " << trg_p << "\n";
        std::cout << " dist(): " << g[src].dist() << " - " << g[trg].dist() << "\n";
    }
    
    /*
    if (s1->isLine() && s2->isLine() ) {
        std::cout << s1->str2() << "\n";
        std::cout << s2->str2() << "\n";
        std::cout << s3->str2() << "\n";
    }
    boost::array<double,8> x = g[edge].x;
    boost::array<double,8> y = g[edge].y;
    for (unsigned int n=0 ; n < 8 ; n++ ) {
        std::cout << n << "  " << x[n] << "  " << y[n] << "\n";
    }
    
    for (int n=0 ; n < 20 ; n++ ) {
        VertexError s1_err_functor(g, edge, s1);
        VertexError s2_err_functor(g, edge, s2);
        double ts = t_min + ((t_max-t_min)/(20-1))*n;
        std::cout << n << " ts=" << ts << " s3_error= " << err_functor(ts) << " p="<< g[edge].point(ts) << "\n";
        std::cout << "     s1_err=" << s1_err_functor(ts) << " s2_err=" << s2_err_functor(ts) << "\n";

    }
    */
    
    typedef std::pair<double, double> Result;
    Result r = boost::math::tools::brent_find_minima( err_functor, t_min, t_max, 64);
    double t_sln = r.first;
    //Point p_sln = g[edge].point(t_sln);
    Point p_sln = err_functor.edge_point(t_sln); //g[edge].point(t_sln);
    double desp_k3(0);
    if (s3->isPoint())
        desp_k3 = 1;
    else if ( s3->isLine() ) {
        // find out on which side the desperate solution lies
        Point src_se = s3->start();
        Point trg_se = s3->end();
        Point left = 0.5*(src_se+trg_se) + (trg_se-src_se).xy_perp(); 
        if (p_sln.is_right(src_se,trg_se)) {
            desp_k3 = (s3->k()==1) ? -1 : 1;
        } else {
            desp_k3 = (s3->k()==1) ? 1 : -1;
        }
    }
    solvers::Solution desp( p_sln, t_sln, desp_k3 ); 
    return desp;
}

/// set debug output true/false
void VertexPositioner::solver_debug(bool b) {
    ppp_solver->set_debug(b);
    lll_solver->set_debug(b);
    qll_solver->set_debug(b);
    sep_solver->set_debug(b);
    alt_sep_solver->set_debug(b);
    lll_para_solver->set_debug(b);
}

void VertexPositioner::set_silent(bool b) {
    silent=b;
    ppp_solver->set_silent(b);
    lll_solver->set_silent(b);
    lll_para_solver->set_silent(b);
    qll_solver->set_silent(b);
    sep_solver->set_silent(b);
    alt_sep_solver->set_silent(b);
}
    
/// dispatch to the correct solver based on the sites
int VertexPositioner::solver_dispatch(Site* s1, double k1, Site* s2, double k2, Site* s3, double k3, 
                                        std::vector<solvers::Solution>& solns) {


    if ( g[edge].type == SEPARATOR ) {
        // this is a SEPARATOR edge with two LineSites adjacent.
        // find the PointSite that defines the SEPARATOR, so that one LineSite and one PointSite
        // can be submitted to the Solver.
        if ( s1->isLine() && s2->isLine() ) {
            // the parallell lineseg case      v0 --s1 --> pt -- s2 --> v1
            // find t
            if ( g[edge].has_null_face ) {
                s2 = g[ g[edge].null_face ].site;
                assert( s2->isPoint() ); // the sites of null-faces are allwais PointSite
                k2 = +1;
            } else if ( g[ g[edge].twin ].has_null_face ) {
                s2 = g[ g[ g[edge].twin ].null_face ].site;
                assert( s2->isPoint() );
                k2 = +1;
            }
        } else if ( s1->isPoint() && s2->isLine() ) {
            // a normal SEPARATOR edge, defined by a PointSite and a LineSite 
            // swap sites, so SEPSolver can assume s1=line s2=point
            Site* tmp = s1;
            double k_tmp = k1;
            s1 = s2;
            s2 = tmp;
            k1 = k2;
            k2 = k_tmp;
            assert( s1->isLine() );
            assert( s2->isPoint() );
        }
        assert( s1->isLine() && s2->isPoint() ); // we have previously set s1(line) s2(point)
        return sep_solver->solve(s1,k1,s2,k2,s3,k3,solns); 
    } else if ( g[edge].type == PARA_LINELINE  && s3->isLine() ) { // an edge betwee parallel LineSites
        //std::cout << " para lineline! \n";
        return lll_para_solver->solve( s1,k1,s2,k2,s3,k3, solns );
    } else if ( s1->isLine() && s2->isLine() && s3->isLine() ) 
        return lll_solver->solve( s1,k1,s2,k2,s3,k3, solns ); // all lines.
    else if ( s1->isPoint() && s2->isPoint() && s3->isPoint() )
        return ppp_solver->solve( s1,1,s2,1,s3,1, solns ); // all points, no need to specify k1,k2,k3, they are all +1
    else if ( (s3->isLine() && s1->isPoint() ) || 
              (s1->isLine() && s3->isPoint() ) ||
              (s3->isLine() && s2->isPoint() ) ||
              (s2->isLine() && s3->isPoint() ) // bad coverage for this line?
            ) {
        // if s1/s2 form a SEPARATOR-edge, this is dispatched automatically to sep-solver
        // here we detect for a separator case between
        // s1/s3
        // s2/s3
        if (s3->isLine() && s1->isPoint() ) {
            if ( detect_sep_case(s3,s1) ) {
                alt_sep_solver->set_type(0);
                return alt_sep_solver->solve(s1, k1, s2, k2, s3, k3, solns );
            }
        }
        if (s3->isLine() && s2->isPoint() ) {
            if ( detect_sep_case(s3,s2) ) {
                alt_sep_solver->set_type(1);
                return alt_sep_solver->solve(s1, k1, s2, k2, s3, k3, solns );
            }
        }
    } 
    
    // if we didn't dispatch to a solver above, we try the general solver
    return qll_solver->solve( s1,k1,s2,k2,s3,k3, solns ); // general case solver
    
}

/// detect separator-case, so we can dispatch to the correct Solver
bool VertexPositioner::detect_sep_case(Site* lsite, Site* psite) {
    HEEdge le = lsite->edge();
    HEVertex src = g.source(le);
    HEVertex trg = g.target(le);
    //std::cout << " detect_sep_case() Linesite from " << g[src].index << " to " << g[trg].index << "\n";
    // now from segment end-points get the null-vertex
    HEEdge src_out;
    BOOST_FOREACH(HEEdge e, g.out_edge_itr(src) ) {
        if ( g[e].type == NULLEDGE )
            src_out = e;
    }
    HEEdge trg_out;
    BOOST_FOREACH(HEEdge e, g.out_edge_itr(trg) ) {
        if ( g[e].type == NULLEDGE )
            trg_out = e;
    }
    //std::cout << " detect_sep_case() src null-edge "; g.print_edge(src_out); // << g[src].index << " to " << g[trg].index << "\n";
    //std::cout << " detect_sep_case() trg null-edge "; g.print_edge(trg_out);
    
    HEFace src_null_face = g[src_out].face;
    if (g[src_null_face].null == false ) {
        // take twin face instead
        HEEdge src_out_twin = g[src_out].twin;
        src_null_face = g[src_out_twin].face;
    }
    
    HEFace trg_null_face = g[trg_out].face;
    if ( g[trg_null_face].null == false ) {
        HEEdge trg_out_twin = g[trg_out].twin;
        trg_null_face = g[trg_out_twin].face;
    }
    assert( g[src_null_face].null && g[trg_null_face].null );
        
    // do we want src_out face??
    // OR src_out_twin face??
    // we want the null-face !
        
    //std::cout << " detect_sep_case() src null-face " << src_null_face << "\n";
    //g.print_face(src_null_face);
    //std::cout << " detect_sep_case() trg null-face " << trg_null_face << "\n";
    //g.print_face(trg_null_face);
    
    Site* src_site = g[src_null_face].site;
    Site* trg_site = g[trg_null_face].site;
    if (src_site == NULL || trg_site == NULL ) {
        exit(-1);
    }
    if ( !src_site->isPoint() || !trg_site->isPoint() ) {
        exit(-1);
    }
    //std::cout << " detect_sep_case() src PointSite is "  << src_site->str() << "\n"; // << src_site->vertex();
    //std::cout << " detect_sep_case() trg PointSite is "  << trg_site->str() << "\n";
    HEVertex src_vertex = src_site->vertex();
    HEVertex trg_vertex = trg_site->vertex();
    //std::cout << " detect_sep_case() 1st end-point vertex is " << g[src_vertex].index << "\n";
    //std::cout << " detect_sep_case() 2nd end-point vertex is " << g[trg_vertex].index << "\n";
    //std::cout << " detect_sep_case() psite is " << psite->str() << "\n";
    //std::cout << " detect_sep_case() psite is " << g[psite->vertex()].index << "\n"; // g[ psite->vertex()].index << "\n";
    if ( src_vertex == psite->vertex() ) {
        //std::cout << " detect_sep_case(): src separator case!\n";
        //std::cout << " detect_sep_case(): line is " << g[src].index << " - " << g[trg].index << " with psites " << g[src_vertex].index << " - " << g[trg_vertex].index << "\n";
        //std::cout << " detect_sep_case(): psite vertex is " << g[ psite->vertex() ].index << "\n";
        return true;
    }
    if ( trg_vertex == psite->vertex() ) {
        //std::cout << " detect_sep_case(): trg separator case!\n";
        //std::cout << " detect_sep_case(): line is " << g[src].index << " - " << g[trg].index << " with psites " << g[src_vertex].index << " - " << g[trg_vertex].index << "\n";
        //std::cout << " detect_sep_case(): psite vertex is " << g[ psite->vertex() ].index << "\n";

        return true;
    }
    //std::cout << " detect_sep_case()   NOT a separator case.\n";
    return false;
}

/// error from solution to corresponding point on the edge
double VertexPositioner::edge_error(solvers::Solution& sl) {
    Point p;
    if (g[edge].type==PARA_LINELINE) {
        p = projection_point( sl );
    } else {
        p = g[edge].point( sl.t );
    }
    return (p-sl.p).norm();
}

/// when the edge is not parametrized by t-value as normal edges
/// so we need a projection of sl onto the edge instead
Point VertexPositioner::projection_point(solvers::Solution& sl) {
    assert( g[edge].type == PARA_LINELINE );
    // edge given by
    // p = p0 + t * (p1-p0)   with t in [0,1]
    Point p0( g[ g.source(edge) ].position );
    Point p1( g[ g.target(edge) ].position ); 
    //std::cout << " edge is  " << p0 << " - " << p1 << "\n";
    //std::cout << " edge direction: " << v << "\n";
    Point v = p1-p0;
    
    double t = (sl.p - p0).dot(v) / v.dot(v);
    // clamp to [0,1]
    if ( t>1)
        t=1;
    else if (t<0)
        t=0;
    //std::cout << " projection of solution " << sl.p << " is " << (p0+v*t) << "\n";
    return (p0+v*t);
}

/// check that the new solution lies on the edge
bool VertexPositioner::solution_on_edge(solvers::Solution& s) {
    double err = edge_error(s);
    double limit = 9E-4;
    if ( err>=limit ) {
        std::cout << "solution_on_edge() ERROR err= " << err << "\n";
        std::cout << " solution edge: " << g[ g.source(edge) ].index << "[" << g[ g.source(edge) ].type<<"]{" << g[ g.source(edge) ].status<<"}";
        std::cout << " -[" << g[edge].type << "]- ";
        std::cout << g[ g.target(edge) ].index << "[" << g[ g.target(edge) ].type << "]{" << g[ g.target(edge) ].status<<"}\n";


        std::cout << " edge: " << g[ g.source(edge) ].index << "(t=" << g[ g.source(edge) ].dist() << ")"; 
        std::cout << " - " << g[ g.target(edge) ].index << "(t=" << g[ g.target(edge) ].dist() << ")\n";
        std::cout << " edge: " << g[ g.source(edge) ].position << " - " << g[ g.target(edge) ].position << "\n";
        std::cout << " solution: " << s.p << " t=" << s.t << "\n";
    }
    return (err<limit);
}

// calculate the distance from the solution-point to the corresponding point on the edge.
/*
double VertexPositioner::edge_error(HEEdge e, Solution& s) {
    Point ep = g[e].point( s.t, s );
    return (ep-s.p).norm();
}*/

/// new vertices should lie within the far_radius
bool VertexPositioner::check_far_circle(solvers::Solution& s) {
    if (!(s.p.norm() < 18*1)) {
        std::cout << "WARNING check_far_circle() new vertex outside far_radius! \n";
        std::cout << s.p << " norm=" << s.p.norm() << " far_radius=" << 1 << "\n"; 
        return false;
    }
    return true;
}

/// distance sanity check
// all vertices should be of degree three, i.e. three adjacent faces/sites
// distance to the three adjacent sites should be equal
bool VertexPositioner::check_dist(HEEdge e, const solvers::Solution& sl, Site* s3) {
    HEFace face = g[e].face;     
    HEEdge tw_edge = g[e].twin;
    HEFace twin_face = g[tw_edge].face;      
    
    Site* s1 = g[face].site;
    Site* s2 = g[twin_face].site;
    
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
        std::cout << " solution edge: " << g[ g.source(edge) ].index << "[" << g[ g.source(edge) ].type<<"]{" << g[ g.source(edge) ].status<<"}";
        std::cout << " -[" << g[edge].type << "]- ";
        std::cout << g[ g.target(edge) ].index << "[" << g[ g.target(edge) ].type << "]{" << g[ g.target(edge) ].status<<"}\n";
    
        return false;
    }
    return true;
}

/// distance-error
// new vertices should be equidistant to the three adjacent sites that define the vertex
// we here calculate the distances d1, d2, d3 from the Solution to the three sites s1, s2, s3
// and return the max deviation from the solution t-value.
// this works as a sanity check for the solver.
// a high error value here is also an indication of numerical instability in the solver
double VertexPositioner::dist_error(HEEdge e, const solvers::Solution& sl, Site* s3) {
    HEFace face = g[e].face;     
    HEEdge tw_edge = g[e].twin;
    HEFace twin_face = g[tw_edge].face;      
    
    Site* s1 = g[face].site;
    Site* s2 = g[twin_face].site;
    
    double d1 = (sl.p - s1->apex_point(sl.p) ).norm();
    double d2 = (sl.p - s2->apex_point(sl.p) ).norm();  
    double d3 = (sl.p - s3->apex_point(sl.p) ).norm(); 
    
    return std::max( std::max( fabs(sl.t-d1),fabs(sl.t-d2)) , fabs(sl.t-d3));

}

/// are \a d1 and \a d2 roughly equal?
bool VertexPositioner::equal(double d1, double d2) {
    double tol = 1e-3;
    if ( fabs(d1-d2) < 1e-15 )
        return true;
    if ( fabs(d1-d2) > tol*std::max(d1,d2) )
        return false;
    return true;
}
    
    
} // end namespace
