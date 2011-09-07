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

//#define NDEBUG  // this turns off assertions
#include <cassert>

#include <boost/foreach.hpp>

#include "voronoidiagram.hpp"
#include "facegrid.hpp"
#include "voronoidiagram_checker.hpp"

namespace ovd {

VoronoiDiagram::VoronoiDiagram(double far, unsigned int n_bins) {
    fgrid = new FaceGrid(far, n_bins);
    vd_checker = new VoronoiDiagramChecker(this);
    far_radius=far;
    gen_count=3;
    initialize();
}

VoronoiDiagram::~VoronoiDiagram() { 
    delete fgrid; 
    delete vd_checker;
}

void VoronoiDiagram::run() {
    //int n_bins = (int)( sqrt(2)*sqrt( vertex_sites.size() ) );
    //if (n_bins < 10)
    //    n_bins = 10;
        
    //fgrid = new FaceGrid(far_radius, n_bins);
    BOOST_FOREACH( Point p, vertex_sites ) {
        add_vertex_site(p);
    }
}

void VoronoiDiagram::push_vertex_site(const Point& p) {
    vertex_sites.push_back(p);
}


// add one vertex at origo and three vertices at 'infinity' and their associated edges
void VoronoiDiagram::initialize() {
    // add init vertices
    HEVertex v0  = g.add_vertex();
    HEVertex v01 = g.add_vertex();
    HEVertex v02 = g.add_vertex();
    HEVertex v03 = g.add_vertex();
    double far_multiplier = 6;
    g[v01] = VoronoiVertex(Point(             0                 , -3.0*far_radius*far_multiplier    )               , OUT, OUTER);
    g[v02] = VoronoiVertex(Point(  +3.0*sqrt(3.0)*far_radius*far_multiplier/2.0, +3.0*far_radius*far_multiplier/2.0), OUT, OUTER);
    g[v03] = VoronoiVertex(Point(  -3.0*sqrt(3.0)*far_radius*far_multiplier/2.0, +3.0*far_radius*far_multiplier/2.0), OUT, OUTER);
    

        
    out_verts[0]=v01; out_verts[1]=v02; out_verts[2]=v03;
    // the locations of the initial generators:
    Point gen1 = Point( 0, 3.0*far_radius);
    Point gen2 = Point( -3.0*sqrt(3.0)*far_radius/2.0, -3.0*far_radius/2.0 );
    Point gen3 = Point( +3.0*sqrt(3.0)*far_radius/2.0, -3.0*far_radius/2.0 );
    g[v0].position = VertexPositioner::position( gen1, gen2, gen3 );
    
     g[v0].init_dist( gen1 );
    g[v01].init_dist( gen3 );
    g[v02].init_dist( gen1 );
    g[v03].init_dist( gen2 );
    
    
    HEVertex g1 = g.add_vertex();
    HEVertex g2 = g.add_vertex();
    HEVertex g3 = g.add_vertex();
    g[g1] = VoronoiVertex( gen1 , OUT, VERTEXGEN);
    g[g2] = VoronoiVertex( gen2 , OUT, VERTEXGEN);
    g[g3] = VoronoiVertex( gen3 , OUT, VERTEXGEN);
    // add face 1: v0-v1-v2 which encloses gen3
    HEEdge e1 =  g.add_edge( v0 , v01 );   
    HEEdge e2 =  g.add_edge( v01, v02 );
    HEEdge e3 =  g.add_edge( v02, v0  ); 
    HEFace f1 =  g.add_face( ); 
    g[f1].edge = e2;
    g[f1].generator = gen3;
    g[f1].status = NONINCIDENT;
    
    fgrid->add_face( g[f1] );
    g[e1].face = f1;
    g[e2].face = f1;
    g[e3].face = f1;
    g[e1].next = e2;
    g[e2].next = e3;
    g[e3].next = e1;
    g[g3].face = f1;
    
    
    // add face 2: v0-v02-v03 which encloses gen1
    HEEdge e4 = g.add_edge( v0, v02   );   
    HEEdge e5 = g.add_edge( v02, v03  );
    HEEdge e6 = g.add_edge( v03, v0   ); 
    HEFace f2 =  g.add_face();
    g[f2].edge = e5;
    g[f2].generator = gen1;
    g[f2].status = NONINCIDENT;
    
    fgrid->add_face( g[f2] );
    g[e4].face = f2;
    g[e5].face = f2;
    g[e6].face = f2;
    g[e4].next = e5;
    g[e5].next = e6;
    g[e6].next = e4;
    g[g1].face = f2;
    
    // add face 3: v0-v3-v1 which encloses gen2
    HEEdge e7 = g.add_edge( v0 , v03 );   
    HEEdge e8 = g.add_edge( v03, v01 );
    HEEdge e9 = g.add_edge( v01, v0  ); 
    HEFace f3 =  g.add_face();
    g[f3].edge = e8;
    g[f3].generator = gen2;
    g[f3].status = NONINCIDENT;
    fgrid->add_face( g[f3] );
    g[e7].face = f3;
    g[e8].face = f3;
    g[e9].face = f3;
    g[e7].next = e8;
    g[e8].next = e9;
    g[e9].next = e7;
    g[g2].face = f3;
    
    // twin edges
    g[e1].twin = e9;
    g[e9].twin = e1;
    g[e2].twin = HEEdge(); // the outermost edges have invalid twins
    g[e5].twin = HEEdge();
    g[e8].twin = HEEdge();
    g[e3].twin = e4;
    g[e4].twin = e3;
    g[e6].twin = e7;
    g[e7].twin = e6;

    assert( vd_checker->isValid() );
}

// comments relate to Sugihara-Iri 1994 paper
// this is roughly "algorithm A" from the paper, page 15/50
int VoronoiDiagram::add_vertex_site(const Point& p) {
    HEVertex new_vert = g.add_vertex( );
    g[new_vert].position=p; 
    g[new_vert].type=VERTEXGEN; 
    g[new_vert].status=OUT; 
    assert( p.norm() < far_radius );     // only add vertices within the far_radius circle
    gen_count++;
    // 1)
    HEFace closest_face = fgrid->grid_find_closest_face( p ); // faster grid-search
    // 2)
    HEVertex v_seed = find_seed_vertex(closest_face, p) ;
    // 3)
    augment_vertex_set(v_seed, p); 
    // 4)
    add_new_voronoi_vertices( p );    
    // 5)
    HEFace newface = split_faces( p );
    // 6)
    remove_vertex_set( newface );
    g[new_vert].face = newface;
    // 7
    reset_status();
    assert( vd_checker->isValid() );
    return g[new_vert].index;
}

// evaluate H on all face vertices and return vertex with the lowest H
/*
HEVertex VoronoiDiagram::find_seed_vertex(HEFace f, const Point& p) {
    VertexVector face_verts = g.face_vertices(f);
    assert( face_verts.size() >= 3 );
    double minimumH(0.0); 
    HEVertex minimalVertex=  HEVertex() ;
    bool first = true;
    BOOST_FOREACH( HEVertex q, face_verts) { // go thorugh all the vertices and find the one with smallest detH
        if ( g[q].status != OUT ) {
            double h = g[q].detH( p ); 
            std::cout << "    H: " << g[q].index << " detH= " << g[q].detH(p) << " r=" << g[q].inCircle(p) << "\n";
            if ( first || (h<minimumH) ) {
                minimumH = h;
                minimalVertex = q;
                first = false;
            }
        }
    }
    std::cout << "RETURN H: " << g[minimalVertex].index << " detH= " << g[minimalVertex].detH(p) << " r=" << g[minimalVertex].inCircle(p) << "\n";
    assert( vd_checker->detH_is_negative( p, f, minimalVertex ) );
    return minimalVertex;
}*/


// evaluate H on all face vertices and return vertex with the lowest H
HEVertex VoronoiDiagram::find_seed_vertex(HEFace f, const Point& p) {
    VertexVector face_verts = g.face_vertices(f);
    assert( face_verts.size() >= 3 );
    double minPred(0.0); 
    HEVertex minimalVertex=  HEVertex() ;
    bool first = true;
    BOOST_FOREACH( HEVertex q, face_verts) { // go thorugh all the vertices and find the one with smallest detH
        if ( g[q].status != OUT ) {
            double h = g[q].in_circle( p ); 
            if ( first || (h<minPred) ) {
                minPred = h;
                minimalVertex = q;
                first = false;
            }
        }
    }
    assert( vd_checker->inCircle_is_negative( p, f, minimalVertex ) );
    return minimalVertex;
}

// growing the v0/delete-tree of "IN" vertices by "weighted breadth-first search"
// we start at the seed and add vertices with detH<0 provided that:
// (C4) v should not be adjacent to two or more IN vertices (this would result in a loop/cycle!)
// (C5) for an incident face containing v: v is adjacent to an IN vertex on this face
//  we process UNDECIDED vertices adjacent to known IN-vertices in a "weighted breadth-first-search" manner
//  where vertices with a large fabs(detH) are processed first, since we assume the detH to be more reliable the larger fabs(detH) is.
void VoronoiDiagram::augment_vertex_set(HEVertex& v_seed, const Point& p) {
    mark_vertex( v_seed, p );
    modified_vertices.push_back( v_seed );
    while( !Q.empty() ) {
        HEVertex v;
        double h;
        boost::tie( v, h ) = Q.top();      assert( g.g[v].status == UNDECIDED );
        Q.pop(); 
        if ( h < 0.0 ) { // mark IN if detH<0 and passes (C4) and (C5) tests. otherwise mark OUT
            if ( (adjacent_in_count(v) >= 2) || (!incidentFacesHaveAdjacentInVertex(v)) ) 
                g[v].status = OUT; // C4 or C5 violated, so mark OUT
            else
                mark_vertex( v,  p); // h<0 and no violations, so mark IN. push adjacent UNDECIDED vertices onto Q.
        } else { 
            g[v].status = OUT; // detH was positive (or zero), so mark OUT
        }
        modified_vertices.push_back( v );
    }
    // sanity-check: for all incident_faces the IN-vertices should be connected
    assert( vd_checker->incidentFaceVerticesConnected(  IN ) );
    assert( vd_checker->incidentFaceVerticesConnected(  OUT ) );
}

void VoronoiDiagram::mark_vertex(HEVertex& v,  const Point& p) {
    g[v].status = IN;
    v0.push_back( v );
    mark_adjacent_faces( v );
    push_adjacent_vertices( v ,  p);
}

void VoronoiDiagram::push_adjacent_vertices( HEVertex v, const Point& p) {
    BOOST_FOREACH( HEVertex w, g.adjacent_vertices(v) ) {
        if ( (g[w].status == UNDECIDED) && (!g[w].in_queue) ) {
                Q.push( VertexDetPair(w , g[w].in_circle(p) ) ); // push adjacent undecided verts for testing.
                g[w].in_queue=true;
        }
    }
}


// IN-Vertex v has three adjacent faces, mark nonincident faces incident
// and push them to incident_faces
void VoronoiDiagram::mark_adjacent_faces( HEVertex v) {
    assert( g[v].status == IN );
    FaceVector new_adjacent_faces = g.adjacent_faces( v ); 
    assert( new_adjacent_faces.size()==3 );
    BOOST_FOREACH( HEFace adj_face, new_adjacent_faces ) {
        if ( g[adj_face].status  != INCIDENT ) {
            g[adj_face].status = INCIDENT; 
            incident_faces.push_back(adj_face);
        }
    }
}


// the set v0 are IN vertices that should be removed
// generate new voronoi-vertices on all edges connecting v0 to OUT-vertices
void VoronoiDiagram::add_new_voronoi_vertices( const Point& p ) {
    assert( !v0.empty() );
    EdgeVector q_edges = find_in_out_edges(); //, OUT); // new vertices generated on these IN-OUT edges
    
    for( unsigned int m=0; m<q_edges.size(); ++m )  {  // create new vertices on all edges q_edges[]
        HEVertex q = g.add_vertex();
        g[q].status = NEW;
        modified_vertices.push_back(q);
        HEFace face = g[q_edges[m]].face;     assert(  g[face].status == INCIDENT);
        HEEdge twin = g[q_edges[m]].twin;
        HEFace twin_face = g[twin].face;      assert( g[twin_face].status == INCIDENT);

        g[q].position = VertexPositioner::position( g[face].generator  , g[twin_face].generator  , p );
        g[q].init_dist(p);
        check_vertex_on_edge(q, q_edges[m]);
        g.insert_vertex_in_edge( q, q_edges[m] );
    }
}

// check that vertex q is positioned on the edge e

void VoronoiDiagram::check_vertex_on_edge(HEVertex q, HEEdge e) {
    // sanity check on new vertex
    if (!(g[q].position.norm() < 18*far_radius)) {
        std::cout << "WARNING check_vertex_on_edge() new vertex outside far_radius! \n";
        std::cout << g[q].position << " norm=" << g[q].position.norm() << " 6.1*far_radius=" << 6.1*far_radius << "\n"; //"WARNING check_vertex_on_edge() new vertex outside far_radius! \n";
    }
    assert( g[q].position.norm() < 18*far_radius);
    
    HEVertex trg = g.target(e);
    HEVertex src = g.source(e);
    Point trgP = g[trg].position;
    Point srcP = g[src].position;
    Point newP = g[q].position;
    
    {
        double dtl_orig = g[q].position.xyDistanceToLine(srcP, trgP);
        //double dtl(dtl_orig);
        if (dtl_orig > 1e-3* ( trgP - srcP ).norm() ) {
            double t = ( g[q].position - srcP).dot( trgP - srcP ) / ( trgP - srcP ).dot( trgP - srcP ) ;
            //g[q].position = srcP + t*( trgP-srcP);
            //dtl = g[q].position.xyDistanceToLine(srcP, trgP);
            std::cout << "WARNING!! check_vertex_on_edge()  dtl= " << dtl_orig << " t= " << t  << " edgelength= " << ( trgP - srcP ).norm()   <<"\n";
        }
    }
    
    
    if (( trgP - srcP ).norm() <= 0 ) {
        std::cout << "WARNING check_vertex_on_edge() zero-length edge! \n";
        g[q].position = srcP;
    } else {
        assert( ( trgP - srcP ).norm() > 0.0 ); // edge has finite length
        assert( ( trgP - srcP ).dot( trgP - srcP ) > 0.0 ); // length squared
        double torig = ((newP - srcP).dot( trgP - srcP )) / ( trgP - srcP ).dot( trgP - srcP ) ;
        bool warn = false;
        double t(torig);
        if (torig < 0.0) { // clamp the t-parameter to [0,1]
            warn = true;
            t=0.0;
        } else if (torig> 1.0) {
            warn = true;
            t=1.0;
        }
        if ( warn ) {
            std::cout << "WARNING: check_vertex_on_edge() t_old= " << torig << " CORRECTED t_new= " << t << "\n";
            std::cout << "src= " << srcP << " new= " << newP << " trg= " << trgP << "\n";
            std::cout << " (src-trg).norm()= " << (srcP-trgP).norm() << "\n";
            g[q].position = srcP + t*( trgP-srcP);
            t = ( g[q].position - srcP).dot( trgP - srcP ) / ( trgP - srcP ).dot( trgP - srcP ) ;
            
        }
        // now we are clamped:
        assert( t >= 0.0 );
        assert( t <= 1.0 );        
        double dtl_orig = g[q].position.xyDistanceToLine(srcP, trgP);
        double dtl(dtl_orig);
        if (dtl_orig > 1e-3* ( trgP - srcP ).norm() ) {
            t = ( g[q].position - srcP).dot( trgP - srcP ) / ( trgP - srcP ).dot( trgP - srcP ) ;
            g[q].position = srcP + t*( trgP-srcP);
            dtl = g[q].position.xyDistanceToLine(srcP, trgP);
            std::cout << "WARNING check_vertex_on_edge()  old_dtl= " << dtl_orig << " new_dtl= " << dtl  << " edgelength= " << ( trgP - srcP ).norm()   <<"\n";
        }
        assert( dtl < 1e-3* ( trgP - srcP ).norm() );
    }
}

HEFace VoronoiDiagram::split_faces(const Point& p) {
    HEFace newface =  g.add_face(); 
    g[newface].generator = p;
    g[newface].status = NONINCIDENT;
    fgrid->add_face( g[newface] );
    BOOST_FOREACH( HEFace f, incident_faces ) {
        split_face(newface, f); // each INCIDENT face is split into two parts: newface and f
    }
    return newface;
}

// 1) repair the next-pointers for newface that are broken.
// 2) remove IN vertices in the set v0
void VoronoiDiagram::remove_vertex_set( HEFace newface ) {
    HEEdge current_edge = g[newface].edge; 
    HEEdge start_edge = current_edge;
    do {
        HEVertex current_target = g.target( current_edge ); // an edge on the new face
        HEVertex current_source = g.source( current_edge );
        BOOST_FOREACH( HEEdge edge, g.out_edges( current_target ) ) { // loop through potential "next" candidates
            HEVertex out_target = g.target( edge );
            if ( g[out_target].status == NEW ) { // the next vertex along the face should be "NEW"
                if ( out_target != current_source ) { // but not where we came from
                    g[current_edge].next = edge; // this is the edge we want to take
                    assert( vd_checker->current_face_equals_next_face( current_edge ) );
                }
            }
        }
        current_edge = g[current_edge].next; // jump to the next edge
    } while (g[current_edge].next != start_edge);
    
    BOOST_FOREACH( HEVertex v, v0 ) {      // it should now be safe to delete all IN vertices
        assert( g[v].status == IN );
        g.delete_vertex(v); // this also removes edges connecting to v
    }
}

// reset status of modified_vertices and incident_faces
void VoronoiDiagram::reset_status() {
    BOOST_FOREACH( HEVertex v, modified_vertices ) {
        g[v].reset();
    }
    modified_vertices.clear();
    g[out_verts[0]].status = OUT; // the outer vertices are special.
    g[out_verts[1]].status = OUT;
    g[out_verts[2]].status = OUT;
    BOOST_FOREACH(HEFace f, incident_faces ) { 
        g[f].status = NONINCIDENT; 
    }
    incident_faces.clear();
    v0.clear();
}


// among the edges of HEFace f
// find the NEW vertex with status-signature NEW->s
boost::tuple<HEEdge, HEVertex, HEEdge> VoronoiDiagram::find_new_vertex(HEFace f, VoronoiVertexStatus s) {
    HEVertex v;
    HEEdge prev;
    HEEdge twin_next;
    HEEdge current_edge = g[f].edge;
    bool found = false;                             
    while (!found) {
        HEVertex current_vertex = g.target( current_edge );
        HEEdge next_edge = g[current_edge].next;
        HEVertex next_vertex = g.target( next_edge );
        if ( g[current_vertex].status == s ) {
            if ( g[next_vertex].status == NEW) {
                v = next_vertex;
                prev = next_edge;
                twin_next = g[next_edge].next;
                found = true;
            }
        }
        current_edge = g[current_edge].next;   
    }
    return boost::tuple<HEEdge, HEVertex, HEEdge>( prev, v, twin_next );
}


// split the face f into one part which is newface, and the other part is the old f
void VoronoiDiagram::split_face(HEFace newface, HEFace f) {
    HEVertex new_source; // this Vertex is found as OUT-NEW-IN
    HEVertex new_target; // this Vertex is found as IN-NEW-OUT
    HEEdge new_previous, new_next, twin_next, twin_previous;
    boost::tie( new_previous, new_source, twin_next) = find_new_vertex(f, OUT);
    boost::tie( twin_previous, new_target, new_next) = find_new_vertex(f, IN);
    // now connect:   new_previous -> new_source -> new_target -> new_next
    // and:              twin_next <- new_source <- new_target <- twin_previous    
    HEEdge e_new = g.add_edge( new_source, new_target );
    g[e_new].next = new_next;
    g[e_new].face = f;
    g[new_previous].next = e_new;
    g[f].edge = e_new; 
    // the twin edge that bounds the new face
    HEEdge e_twin = g.add_edge( new_target, new_source );
    g[e_twin].next = twin_next; 
    g[e_twin].face = newface;
    g[twin_previous].next = e_twin;
    g[newface].edge = e_twin; 
    g.twin_edges(e_new,e_twin);
}

// given a list inVertices of "IN" vertices, find the adjacent IN-OUT edges 
EdgeVector VoronoiDiagram::find_in_out_edges() { 
    assert( !v0.empty() );
    EdgeVector output; // new vertices generated on these edges
    BOOST_FOREACH( HEVertex v, v0 ) {                                   
        assert( g[v].status == IN ); // all verts in v0 are IN
        BOOST_FOREACH( HEEdge edge, g.out_edges( v ) ) {
            HEVertex adj_vertex = g.target( edge );
            if ( g[adj_vertex].status == OUT ) 
                output.push_back(edge); // this is an IN-OUT edge
        }
    }
    assert( !output.empty() );
    return output;
}


int VoronoiDiagram::adjacent_in_count(HEVertex v) {
    int in_count=0;
    BOOST_FOREACH( HEVertex w, g.adjacent_vertices(v) ) {
        if ( g[w].status == IN )
            in_count++;
    }
    return in_count;
}

// any voronoi-vertex v has three adjacent faces
// return those that are INCIDENT
FaceVector VoronoiDiagram::adjacent_incident_faces(HEVertex v) {
    FaceVector adj_faces = g.adjacent_faces(v);   assert( adj_faces.size() == 3 );
    FaceVector inc_faces;
    BOOST_FOREACH( HEFace f, adj_faces ) {
        if ( g[f].status == INCIDENT )
            inc_faces.push_back( f );
    }
    assert( !inc_faces.empty() );
    return inc_faces;
}

bool VoronoiDiagram::incidentFacesHaveAdjacentInVertex(HEVertex v) {
    bool all_found = true;
    BOOST_FOREACH( HEFace f, adjacent_incident_faces(v) ) { // check each face f
        // v should be adjacent to an IN vertex on the face
        bool face_found=false;
        BOOST_FOREACH( HEVertex w, g.face_vertices(f) ) {
            if ( w != v && g[w].status == IN && g.has_edge(w,v) ) 
                face_found = true;
        }
        if (!face_found)
            all_found=false;
    }
    return all_found;
}



void VoronoiDiagram::print_face(HEFace f) {
    std::cout << " Face " << f << ": ";
    VertexVector face_verts = g.face_vertices(f);    
    unsigned count=1;
    BOOST_FOREACH( HEVertex v, face_verts ) {
        std::cout << g[v].index  << "(" << g[v].status  << ")";
        if (count != face_verts.size() )
            std::cout << "-";
        count++;
    }
    std::cout << "\n";
}

void VoronoiDiagram::print_vertices(VertexVector& q) {
    BOOST_FOREACH( HEVertex v, q) {
        std::cout << g[v].index << " ";
    }
    std::cout << std::endl;
}

std::string VoronoiDiagram::str() const {
    std::ostringstream o;
    o << "VoronoiDiagram (nVerts="<< g.num_vertices() << " , nEdges="<< g.num_edges() <<"\n";
    return o.str();
}

} // end namespace
// end file voronoidiagram.cpp
