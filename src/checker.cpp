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

#include "checker.hpp"
#include "voronoidiagram.hpp"

namespace ovd {

VoronoiDiagramChecker::VoronoiDiagramChecker(HEGraph& gi) : g(gi) {}

VoronoiDiagramChecker::~VoronoiDiagramChecker() {}

/// overall sanity-check for the diagram, calls other sanity-check functions
bool VoronoiDiagramChecker::is_valid() {
    return  (   all_faces_ok() && 
                vertex_degree_ok() &&
                face_count_equals_generator_count()
            );
}

/// check that number of faces equals the number of generators
/// \todo not implemented!
bool VoronoiDiagramChecker::face_count_equals_generator_count() {
    // Euler formula for planar graphs
    // v - e + f = 2
    // in a half-edge diagram all edges occur twice, so:
    // f = 2-v+e
    //int vertex_count = hedi::num_vertices(g);
    /*int vertex_count = 0;
    BOOST_FOREACH( HEVertex v, hedi::vertices( g ) ) {
        if ( g[v].type == NORMAL )
            vertex_count++;
    }
    int face_count = (vertex_count- 4)/2 + 3; // degree three graph
    //int face_count = hed.num_faces();
    if (face_count != gen_count) {
        std::cout << " face_count_equals_generator_count() ERROR:\n";
        std::cout << " num_vertices = " << vertex_count << "\n";
        std::cout << " gen_count = " << gen_count << "\n";
        std::cout << " face_count = " << face_count << "\n";
    }
    return ( face_count == gen_count );
    * */
    return true;
}

/// check that the diagram is of degree three.
/// however ::SPLIT and ::APEX vertices are of degree 2.
bool VoronoiDiagramChecker::vertex_degree_ok() {
    BOOST_FOREACH(HEVertex v, g.vertices() ) {
        if ( g.degree(v) != VoronoiVertex::expected_degree[ g[v].type ] ) {
            std::cout << " vertex_degree_ok() ERROR\n";
            std::cout << " vertex " << g[v].index << " type = " << g[v].type << "\n";
            std::cout << " vertex degree = " << g.degree(v) << "\n";
            std::cout << " expected degree = " << VoronoiVertex::expected_degree[ g[v].type ]  << "\n";
            return false;
        }
    }
    return true;
}

// traverse the incident faces and check next-pointers
/*
bool VoronoiDiagramChecker::allIncidentFacesOK() { // have this take incident_faces as a parameter?
    // all incident faces should pass the sanity-check
    BOOST_FOREACH( HEFace f, incident_faces ) {
        if ( !faceVerticesConnected(  f, IN ) )
            return false; // IN vertices should be connected
        if ( !faceVerticesConnected(  f, OUT ) )  // OUT vertices should be connected
            return false;
        if ( !noUndecidedInFace( f ) )            // no UNDECIDED vertices should remain
            return false;
    }
    return true;
}*/

/// check that all vertices in the input vector have status ::IN
bool VoronoiDiagramChecker::all_in( const VertexVector& q) {
    BOOST_FOREACH( HEVertex v, q) {
        if ( g[v].status != IN )
            return false;
    }
    return true;
}

/// check that no undecided vertices remain in the face
bool  VoronoiDiagramChecker::noUndecidedInFace( HEFace f ) { // is this true??
    VertexVector face_verts = g.face_vertices(f);
    BOOST_FOREACH( HEVertex v, face_verts ) {
        if ( g[v].status == UNDECIDED )
            return false;
    }
    return true;
}

/// check that for HEFace f the vertices TYPE are connected
bool VoronoiDiagramChecker::faceVerticesConnected(  HEFace f, VertexStatus Vtype ) {
    VertexVector face_verts = g.face_vertices(f);
    VertexVector type_verts;
    BOOST_FOREACH( HEVertex v, face_verts ) {
        if ( g[v].status == Vtype )
            type_verts.push_back(v); // build a vector of all Vtype vertices
    }
    assert( !type_verts.empty() );
    if (type_verts.size()==1) // set of 1 is allways connected
        return true;
    
    // check that type_verts are connected
    HEEdge currentEdge = g[f].edge;
    HEVertex endVertex = g.source(currentEdge); // stop when target here
    EdgeVector startEdges;
    bool done = false;
    while (!done) { 
        HEVertex src = g.source( currentEdge );
        HEVertex trg = g.target( currentEdge );
        if ( g[src].status != Vtype ) { // seach ?? - Vtype
            if ( g[trg].status == Vtype ) { // we have found ?? - Vtype
                startEdges.push_back( currentEdge );
            }
        }
        currentEdge = g[currentEdge].next;
        if ( trg == endVertex ) {
            done = true;
        }
    }
    assert( !startEdges.empty() );
    if ( startEdges.size() != 1 ) // when the Vtype vertices are connected, there is exactly one startEdge
        return false;
    else 
        return true;
}

// sanity check: IN-vertices on each INCIDENT face should be connected
/*
bool VoronoiDiagramChecker::incidentFaceVerticesConnected( VoronoiVertexStatus  ) {    
    BOOST_FOREACH( HEFace f1, incident_faces ) {
        if ( !faceVerticesConnected(  f1, IN ) ) {
            std::cout << " VoronoiDiagramChecker::incidentFaceVerticesConnected() ERROR, IN-vertices not connected.\n";
            std::cout << " printing all incident faces for debug: \n";
            BOOST_FOREACH( HEFace f2, incident_faces ) {
                g.print_face( f2 );
            } 
            return false;
        }
    }
    return true;
}*/


/// check that all faces are ok. calls face_ok()
bool VoronoiDiagramChecker::all_faces_ok() {
    for(HEFace f=0;f< g.num_faces() ; f++ ) {
        if (!face_ok(f)) {
            std::cout << "VoronoiDiagramChecker::all_faces_ok() ERROR: f= " << f << "\n";
            g.print_face(f);
            return false;
        }
    }
    return true;
}

/// check that the face is ok
bool VoronoiDiagramChecker::face_ok(HEFace f, bool debug) {
    HEEdge current_edge = g[f].edge;
    //std::cout << " face_ok( " << f << " )\n";
    
    HEEdge start_edge= current_edge;
    double k = g[current_edge].k;
    if ( !((k==1) || (k==-1)) ) {
        std::cout << " VoronoiDiagramChecker::face_ok() f=" << f << " ERROR:\n";
        std::cout << " illegal k-value for edge:";
        std::cout << g[ g.source(current_edge)].index << " - "; 
        std::cout <<  g[ g.target(current_edge)].index  ;
        std::cout << " k= " << k << "\n";
        return false;
    }
    if (g[f].site!=0) { // guard against null-faces that dont have Site
        if ( g[f].site->isPoint() ) {
            if ( !(k==1) ) {
                std::cout << " VoronoiDiagramChecker::face_ok() f=" << f << " ERROR:\n";
                std::cout << " f = " << f << " site is " << g[f].site->str() << " but k=" << k  << "\n";
                std::cout << " null? " << g[f].null << "\n";
                return false;
            }
        }
    }
    int n=0;
    if (debug) std::cout << " checking face " << f << "\n";
    do {
        if(debug) {
            std::cout << " edge: " << g[ g.source(current_edge)].index << " - "; 
            std::cout <<  g[ g.target(current_edge)].index  ;
            std::cout << " edge.face= " << g[ current_edge ].face;
            std::cout << " edge.k= " << g[ current_edge ].k;
            std::cout << "\n";
        }
        if (g[current_edge].k != k )  { // all edges should have the same k-value
            std::cout << " face_ok() g[current_edge].k != k ! \n";
            return false;
        }
        if ( !current_face_equals_next_face(current_edge) )  {// all edges should have the same face
            std::cout << " face_ok() !current_face_equals_next_face(current_edge) ! \n";
            return false;
        }
        
        if ( !check_edge(current_edge) ) {
            std::cout << " VoronoiDiagramChecker::face_ok() f= " << f << " check_edge ERROR\n";
            std::cout << " edge: " << g[ g.source(current_edge)].index << " t=" << g[ g.source(current_edge)].type; 
            std::cout << " - ";
            std::cout <<  g[ g.target(current_edge)].index << " t=" << g[ g.target(current_edge) ].type << "\n"; 
            HEEdge twin = g[current_edge].twin;
            std::cout << " twin: " << g[ g.source(twin)].index  << " t=" << g[ g.source(twin)].type; 
            std::cout << " - ";
            std::cout <<  g[ g.target(twin)].index << " t=" << g[ g.target(twin) ].type << "\n"; 
            //std::cout << " twin: " << g[ g.source(twin)].index << " - " <<  g[ g.target(twin)].index << "\n";
            
            std::cout << " edge-type= " << g[current_edge].type << "\n";

            return false;
        }
        
        current_edge = g[current_edge].next; 
        n++;
        assert( n < 10000 ); // reasonable max
    } while( current_edge != start_edge);
    //std::cout << " face ok, edges=" << n-1 << "\n";
    return true;
}

/// check that current edge and next-edge are on the same face
bool VoronoiDiagramChecker::current_face_equals_next_face( HEEdge e) {
    if ( g[e].face !=  g[ g[e].next ].face) {
        std::cout << " current_face_equals_next_face() ERROR.\n";
        std::cout << "   current.face = " << g[e].face << " IS NOT next_face = " << g[ g[e].next ].face << std::endl;
        HEVertex c_trg = g.target( e );
        HEVertex c_src = g.source( e );

        HEVertex n_trg = g.target( g[e].next );
        HEVertex n_src = g.source( g[e].next );
        
        std::cout << "   current_edge = " << g[c_src].index << " - " << g[c_trg].index << " type=" << g[e].type << " face=" << g[e].face  <<"\n";
        std::cout << "   next_edge = " << g[n_src].index << " - " << g[n_trg].index << " type=" << g[ g[e].next ].type << " face="<< g[ g[e].next ].face << "\n";
        
        g.print_face( g[e].face );
        g.print_face( g[ g[e].next ].face );
        
        //std::cout << " printing all incident faces for debug: \n";
        //BOOST_FOREACH( HEFace f, incident_faces ) {
        //    g.print_face( f );
        //} 
        return false;
    }
    return true;
}

/// sanity-check for edge
bool VoronoiDiagramChecker::check_edge(HEEdge e) const {
    HEVertex src = g.source(e);
    HEVertex trg = g.target(e);
    HEEdge twine = g[e].twin;

    if (twine == HEEdge() ) {
        return true;
    } else if ( !( e == g[twine].twin ) ) {
        std::cout << " VoronoiDiagramChecker::check_edge() twinning error!\n";
        return false;
    }
    
    HEVertex tw_src = g.source(twine);
    HEVertex tw_trg = g.target(twine);
    //std::cout << (src==tw_trg) << " && " << (trg==tw_src) << "\n";
    if ( !((src==tw_trg) && (trg==tw_src)) ) {
        std::cout << "VoronoiDiagramChecker::check_edge() ERROR: \n";
        std::cout << "      edge: " << g[src].index << "("<< g[src].type << ")";
        std::cout << " - " << g[trg].index << "("<< g[trg].type << ")" << "\n";
        //std::cout << "      twin: " << g[tw_src].index << " - " << g[tw_src].index << "\n";
        std::cout << "      edge: " << e << "\n";
        std::cout << "      twin: " << twine << "\n";
        std::cout << "      edge: " << src << " - " << trg << "\n";
        std::cout << "      twin: " << tw_src << " - " << tw_trg << "\n";

    }
    return ( (src==tw_trg) && (trg==tw_src) );
}

} // end namespace
