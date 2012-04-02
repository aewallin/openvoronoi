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

#pragma once

#include <qd/qd_real.h> // http://crd.lbl.gov/~dhbailey/mpdist/

#include "graph.hpp"
#include "vertex.hpp"
#include "solvers/solution.hpp"

namespace ovd {

namespace solvers {
class Solver; // fwd decl
}

/// Calculates the (x,y) position of a VoronoiVertex in the VoronoiDiagram
class VertexPositioner {
public:
    VertexPositioner(HEGraph& gi);
    virtual ~VertexPositioner();
    solvers::Solution position( HEEdge e, Site* s);
    /// return vector of errors
    std::vector<double> get_stat() {return errstat;}
    double dist_error(HEEdge e, const solvers::Solution& sl, Site* s3);
    void solver_debug(bool b);
    void set_silent(bool b); ///< no warning messages when silent==true
private:

    /// predicate for rejecting out-of-region solutions
    struct in_region_filter {
        /// \param s Site for in_region check
        in_region_filter(Site* s): site_(s) {}
        /// is Solution \a s in_region of Site \a site_ ?
        bool operator()(solvers::Solution s) { 
            return !site_->in_region(s.p); 
        }
        /// the Site
        Site* site_;
    };
    /// predicate for filtering solutions based on t-value in [tmin,tmax] range
    struct t_filter {
        /// create filter for [tmin,tmax]
        t_filter(double tmin, double tmax): tmin_(tmin),tmax_(tmax) {}
        /// is the given Solution \a s in the offset-distance interval [tmin,tmax] ?
        bool operator()(solvers::Solution s) { 
            double eps=1e-9;
            double tround=s.t;
            if ( fabs(s.t-tmin_) < eps )
                tround=tmin_;
            else if (fabs(s.t-tmax_)<eps)
                tround=tmax_;
            return (tround<tmin_) || (tround>tmax_); // these points rejected!
        }
        /// minimum offset-distance value
        double tmin_;
        /// maximum offset-distance value
        double tmax_;
    };

    solvers::Solution position(Site* s1, double k1, Site* s2, double k2, Site* s3);
    int solver_dispatch(Site* s1, double k1, 
               Site* s2, double k2, 
               Site* s3, double k3, std::vector<solvers::Solution>& slns ); 
    bool detect_sep_case(Site* lsite, Site* psite);

// solution-filtering
    double edge_error(solvers::Solution& sl);
    Point projection_point(solvers::Solution& sl);
// geometry-checks
    bool solution_on_edge(solvers::Solution& s);
    bool check_far_circle(solvers::Solution& s);
    bool check_dist(HEEdge e, const solvers::Solution& s, Site* s3);
    bool equal(double d1, double d2);
    
    solvers::Solution desperate_solution(Site* s3);

// solvers, to which we dispatch, depending on the input sites
    
    solvers::Solver* ppp_solver; ///< point-point-point solver
    solvers::Solver* lll_solver; ///< line-line-line solver
    solvers::Solver* lll_para_solver; ///< solver
    solvers::Solver* qll_solver; ///< solver
    solvers::Solver* sep_solver; ///< separator solver
    solvers::Solver* alt_sep_solver; ///< alternative separator solver
// DATA
    HEGraph& g;  ///< reference to the VD graph.
    double t_min; ///< minimum offset-distance
    double t_max; ///< maximum offset-distance
    HEEdge edge;  ///< the edge on which we position a new vertex
    std::vector<double> errstat; ///< error-statistics
    bool silent; ///< silent mode (outputs no warnings to stdout)
};

/// \brief error functor for edge-based desperate solver
///
/// minimize error by searching for a point on the solution-edge
class VertexError {
public:
    /// \param gi vd-graph
    /// \param sln_edge solution edge
    /// \param si3 newly inserted Site
    VertexError(HEGraph& gi, HEEdge sln_edge, Site* si3) :
    g(gi),  edge(sln_edge), s3(si3)
    {}
    /// return the vertex-error t-d3 where
    /// t3 is the distance from edge-point(t) to s3, and
    /// t is the offset-distance of the solution
    double operator()(const double t) {
        Point p = edge_point(t);
        double s3_dist = (p - s3->apex_point(p)).norm();
        return fabs(t-s3_dist);
    }
    /// return a point on the edge at given offset-distance
    /// \param t offset-distance ( >= 0 )
    Point edge_point(const double t) {
        Point p;
        if ( g[edge].type == LINELINE ) { // this is a workaround because the LINELINE edge-parameters are wrong? at least in some cases?
            HEVertex src = g.source(edge);
            HEVertex trg = g.target(edge);
            Point src_p = g[src].position;
            Point trg_p = g[trg].position;
            double src_t = g[src].dist();
            double trg_t = g[trg].dist();
            // edge is src_p -> trg_p
            if ( trg_t > src_t ) {
                double frac = (t-src_t) / (trg_t-src_t);
                p = src_p + frac*(trg_p-src_p);
            } else {
                double frac = (t-trg_t) / (src_t-trg_t);
                p = trg_p + frac*(src_p-trg_p);
            }
            
        } else
            p = g[edge].point(t);
        return p;
    }
private:
    HEGraph& g; ///< vd-graph
    HEEdge edge; ///< existing edge on which we have positioned a new vertex
    Site* s3; ///< newly inserted Site
};

}

