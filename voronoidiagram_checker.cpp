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

#include "voronoidiagram_checker.hpp"
#include "voronoidiagram.hpp"

namespace ovd {

    
/// sanity-check for the diagram, calls other sanity-check functions
bool VoronoiDiagramChecker::isValid() {
    if (!vertex_degree_ok() )
        return false;
    if (!face_count_equals_generator_count())
        return false;
    return true;
}
    
    
    
/// check that number of faces equals the number of generators
bool VoronoiDiagramChecker::face_count_equals_generator_count() {
    // Euler formula for planar graphs
    // v - e + f = 2
    // in a half-edge diagram all edges occur twice, so:
    // f = 2-v+e
    //int vertex_count = hedi::num_vertices(vd->g);
    /*int vertex_count = 0;
    BOOST_FOREACH( HEVertex v, hedi::vertices( vd->g ) ) {
        if ( vd->g[v].type == NORMAL )
            vertex_count++;
    }
    int face_count = (vertex_count- 4)/2 + 3; // degree three graph
    //int face_count = hed.num_faces();
    if (face_count != vd->gen_count) {
        std::cout << " face_count_equals_generator_count() ERROR:\n";
        std::cout << " num_vertices = " << vertex_count << "\n";
        std::cout << " gen_count = " << vd->gen_count << "\n";
        std::cout << " face_count = " << face_count << "\n";
    }
    return ( face_count == vd->gen_count );
    * */
    return true;
}

    
    
/// the diagram should be of degree three (at least with point generators)
bool VoronoiDiagramChecker::vertex_degree_ok() {
    // the outermost init() vertices have special degree, all others == 6
    BOOST_FOREACH(HEVertex v, vd->g.vertices() ) {
        
        if ( vd->g.degree(v) != VoronoiVertex::expected_degree[ vd->g[v].type ] )
            return false;
            
       /* 
        if ( vd->g[v].type == NORMAL ) {
            if ( vd->g.degree( v ) != 6 )
                return false;
        }
        if ( vd->g[v].type == OUTER ) {
            if ( vd->g.degree( v ) != 4 )
                return false;
        }
        if ( vd->g[v].type == VERTEXGEN ) {
            if ( vd->g.degree( v ) != 0 )
                return false;
        }*/
    }
    return true;
}
    
    
    /// traverse the incident faces and check next-pointers
    bool VoronoiDiagramChecker::allIncidentFacesOK() {
        // all incident faces should pass the sanity-check
        BOOST_FOREACH( HEFace f, vd->incident_faces ) {
            if ( !faceVerticesConnected(  f, IN ) )
                return false; // IN vertices should be connected
            if ( !faceVerticesConnected(  f, OUT ) )  // OUT vertices should be connected
                return false;
            if ( !noUndecidedInFace( f ) )            // no UNDECIDED vertices should remain
                return false;
        }
        return true;
    }
    
    
    /// check that all vertices in the input vector are of type IN
    bool VoronoiDiagramChecker::allIn( const VertexVector& q) {
        BOOST_FOREACH( HEVertex v, q) {
            if ( vd->g[v].status != IN )
                return false;
        }
        return true;
    }

    /// check that no undecided vertices remain in the face
    bool  VoronoiDiagramChecker::noUndecidedInFace(  HEFace f ) {
        VertexVector face_verts = vd->g.face_vertices(f);
        BOOST_FOREACH( HEVertex v, face_verts ) {
            if ( vd->g[v].status == UNDECIDED )
                return false;
        }
        return true;
    }
        
// check that for HEFace f the vertices TYPE are connected
bool VoronoiDiagramChecker::faceVerticesConnected(  HEFace f, VoronoiVertexStatus Vtype ) {
    VertexVector face_verts = vd->g.face_vertices(f);
    VertexVector type_verts;
    BOOST_FOREACH( HEVertex v, face_verts ) {
        if ( vd->g[v].status == Vtype )
            type_verts.push_back(v); // build a vector of all Vtype vertices
    }
    assert( !type_verts.empty() );
    if (type_verts.size()==1) // set of 1 is allways connected
        return true;
    
    // check that type_verts are connected
    HEEdge currentEdge = vd->g[f].edge;
    HEVertex endVertex = vd->g.source(currentEdge); // stop when target here
    EdgeVector startEdges;
    bool done = false;
    while (!done) { 
        HEVertex src = vd->g.source( currentEdge );
        HEVertex trg = vd->g.target( currentEdge );
        if ( vd->g[src].status != Vtype ) { // seach ?? - Vtype
            if ( vd->g[trg].status == Vtype ) {
                // we have found ?? - Vtype
                startEdges.push_back( currentEdge );
            }
        }
        currentEdge = vd->g[currentEdge].next;
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

bool VoronoiDiagramChecker::incidentFaceVerticesConnected( VoronoiVertexStatus Vtype ) {
    // sanity check: IN-vertices for each face should be connected
    BOOST_FOREACH( HEFace f, vd->incident_faces ) {
        if ( !faceVerticesConnected(  f, IN ) ) {
            std::cout << " VoronoiDiagramChecker::incidentFaceVerticesConnected() ERROR, IN-vertices not connected.\n";
            std::cout << " printing all incident faces for debug: \n";
            BOOST_FOREACH( HEFace f, vd->incident_faces ) {
                vd->print_face( f );
            } 
            return false;
        }
    }
    return true;
}

bool VoronoiDiagramChecker::detH_is_negative(  const Point& p, HEFace f, HEVertex minimalVertex ) {
    double minimumH = vd->g[minimalVertex].detH(p);
    if (!(minimumH < 0) ) {
        std::cout << " VD find_seed_vertex() WARNING\n";
        std::cout << " WARNING: searching for seed when inserting " << p << " \n";
        std::cout << " WARNING: closest face is " << f << " with generator " << vd->g[f].generator << " \n";
        std::cout << " WARNING: minimal vd-vertex " << vd->g[minimalVertex].index << " has detH= " << minimumH  << "\n";
    }
    
    return (minimumH < 0 );
}
    
bool VoronoiDiagramChecker::current_face_equals_next_face( HEEdge e) {
    if ( vd->g[e].face !=  vd->g[ vd->g[e].next ].face) {
        std::cout << " VD remove_vertex_set() ERROR.\n";
        std::cout << "current.face = " << vd->g[e].face << " IS NOT next_face = " << vd->g[ vd->g[e].next ].face << std::endl;
        HEVertex c_trg = vd->g.target( e );
        HEVertex c_src = vd->g.source( e );
        HEVertex n_trg = vd->g.target( vd->g[e].next );
        HEVertex n_src = vd->g.source( vd->g[e].next );
        
        std::cout << "current_edge = " << vd->g[c_src].index << " - " << vd->g[c_trg].index << "\n";
        std::cout << "next_edge = " << vd->g[n_src].index << " - " << vd->g[n_trg].index << "\n";
        
        vd->print_face( vd->g[e].face );
        vd->print_face( vd->g[ vd->g[e].next ].face );
        
        std::cout << " printing all incident faces for debug: \n";
        BOOST_FOREACH( HEFace f, vd->incident_faces ) {
            vd->print_face( f );
        } 
        return false;
    }
    return true;
}
                    
/* OLD CODE NOT USED ANYMORE
int VoronoiDiagram::outVertexCount(HEFace f) {
    int outCount = 0;
    VertexVector face_verts = hedi::face_vertices(f, g);
    BOOST_FOREACH( HEVertex v, face_verts ) {
        if (g[v].status == OUT )
            ++outCount;
    }
    return outCount;
}*/


} // end namespace
