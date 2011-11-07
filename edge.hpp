/* 
 *  Copyright 2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include <cassert>
#include <cmath>

#include <boost/graph/adjacency_list.hpp>

#include "point.hpp"
#include "site.hpp"

namespace ovd {

#define OUT_EDGE_CONTAINER boost::vecS 

// note: cannot use vecS since remove_vertex invalidates iterators/edge_descriptors (?)
#define VERTEX_CONTAINER boost::listS
#define EDGE_LIST_CONTAINER boost::listS

// type of edge-descriptors in the graph
typedef boost::adjacency_list_traits<OUT_EDGE_CONTAINER, 
                                     VERTEX_CONTAINER, 
                                     boost::bidirectionalS, 
                                     EDGE_LIST_CONTAINER >::edge_descriptor HEEdge;
                                     
typedef unsigned int HEFace;

enum VoronoiEdgeType {LINE, LINELINE, OUTEDGE, PARABOLA, ELLIPSE, HYPERBOLA, SEPARATOR, LINESITE};

/// properties of an edge in the voronoi diagram
/// each edge stores a pointer to the next HEEdge 
/// and the HEFace to which this HEEdge belongs
struct EdgeProps {
    EdgeProps();
    EdgeProps(HEEdge n, HEFace f): next(n), face(f) {}
    /// create edge with given next, twin, and face
    EdgeProps(HEEdge n, HEEdge t, HEFace f): next(n), twin(t), face(f) {}
    /// the next edge, counterclockwise, from this edge
    HEEdge next; 
    /// the twin edge
    HEEdge twin;
    /// the face to which this edge belongs
    HEFace face; // each face corresponds to an input Site/generator
    double k; // offset-direction from the adjacent site, either +1 or -1
    VoronoiEdgeType type;
    double x[8];
    double y[8];
    bool sign; // choose either +/- in front of sqrt()
    
    inline double sq(double x) const {return x*x;}
    inline double round(double x) const {
        double eps = 1e-8;
        if (fabs(x) < eps)
            return 0.0;
        else
            return x;
    }
    

    Point point(double t) const; 


    double minimum_t( Site* s1, Site* s2) {
        if (s1->isPoint() && s2->isPoint())        // PP
            return minimum_pp_t(s1,s2);
        else if (s1->isPoint() && s2->isLine())    // PL
            return minimum_pl_t(s1,s2);
        else if (s2->isPoint() && s1->isLine())    // LP
            return minimum_pl_t(s2,s1);
        else if (s1->isLine() && s2->isLine())     // LL
            return 0;
        else
            assert(0);
        
        return -1;
    }
    double minimum_pp_t(Site* s1, Site* s2) {
        assert( s1->isPoint() && s2->isPoint() );
        double p1p2 = (s1->position() - s2->position()).norm() ;
        assert( p1p2 >=0 );
        return p1p2/2;
    }
    double minimum_pl_t(Site* s1, Site* s2) {
        double mint = - x[6]/(2.0*x[7]);
        assert( mint >=0 );
        return mint;
    }
    void copy_parameters(EdgeProps& other);   
    void set_parameters(Site* s1, Site* s2, bool sig);
    void set_pp_parameters(Site* s1, Site* s2);
    void set_pl_parameters(Site* s1, Site* s2);
    void set_ll_parameters(Site* s1, Site* s2);
    void print_params() const;

};

} // end namespace
