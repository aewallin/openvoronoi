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

#include <cassert>

#include <boost/foreach.hpp>

#include "voronoidiagram.hpp"
#include "facegrid.hpp"
#include "voronoidiagram_checker.hpp"

namespace ovd {

VoronoiDiagram::VoronoiDiagram(double far, unsigned int n_bins) {
    fgrid = new FaceGrid(far, n_bins);
    vd_checker = new VoronoiDiagramChecker(this); // helper-class that checks topology/geometry
    vpos = new VertexPositioner(this); // helper-class that positions vertices
    far_radius=far;
    gen_count=3;
    initialize();
}

VoronoiDiagram::~VoronoiDiagram() { 
    delete fgrid; 
    delete vd_checker;
    delete vpos;
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
    g[v0].position = Point(0,0);    // OR we could run vpos->position( gen1, gen2, gen3 );
    
    g[ v0].init_dist( gen1 );
    g[v01].init_dist( gen3 );
    g[v02].init_dist( gen1 );
    g[v03].init_dist( gen2 );
    // add initial sites to graph
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
    g[f1].site = new PointSite(gen3);
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
    g[f2].site = new PointSite(gen1);
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
    g[f3].site = new PointSite(gen2);
    g[f3].status = NONINCIDENT;
    fgrid->add_face( g[f3] );
    g[e7].face = f3;
    g[e8].face = f3;
    g[e9].face = f3;
    g[e7].next = e8;
    g[e8].next = e9;
    g[e9].next = e7;
    g[g2].face = f3;
    
    // set type. (note that edge-params x[8] and y[8] are not set!
    g[e1].type = LINE;
    g[e2].type = LINE; 
    g[e3].type = LINE;
    g[e4].type = LINE;
    g[e5].type = LINE;
    g[e6].type = LINE;
    g[e7].type = LINE;
    g[e8].type = LINE;
    g[e9].type = LINE;
    
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
    
    // k-values all positive for PointSite generators
    g[e1].k = 1.0;
    g[e2].k = 1.0;
    g[e3].k = 1.0;
    g[e4].k = 1.0;
    g[e5].k = 1.0;
    g[e6].k = 1.0;
    g[e7].k = 1.0;
    g[e8].k = 1.0;
    g[e9].k = 1.0;
    
    assert( vd_checker->isValid() );
}

// comments relate to Sugihara-Iri 1994 paper
// this is roughly "algorithm A" from the paper, page 15/50
//
// 1) find the face that is closest to the new site
// 2) among the vertices on the closest face, find the seed vertex
// 3) grow the tree of IN-vertices
// 4) add new voronoi-vertices on all IN-OUT edges so they becone IN-NEW-OUT
// 5) add new face by splitting each INCIDENT face into two parts by inserting a NEW-NEW edge. 
// 6) remove IN-IN edges and IN-NEW edges
// 7) reset vertex/face status to be ready for next incremental operation
int VoronoiDiagram::insert_point_site(const Point& p) {
    assert( p.norm() < far_radius );     // only add vertices within the far_radius circle
    
    HEVertex new_vert = g.add_vertex();
    g[new_vert].position=p; 
    g[new_vert].type=VERTEXGEN; 
    g[new_vert].status=OUT; 
    g[new_vert].site = new PointSite(p);
    gen_count++;
    HEFace closest_face = fgrid->grid_find_closest_face( p ); 
    HEVertex v_seed = find_seed_vertex(closest_face, g[new_vert].site ) ;
    augment_vertex_set(v_seed, g[new_vert].site );
    add_new_vertices( g[new_vert].site );    
    HEFace newface = add_new_face( g[new_vert].site );
    remove_vertex_set( newface );
    g[new_vert].face = newface;
    
    reset_status();
    assert( vd_checker->isValid() );
    return g[new_vert].index;
}
void VoronoiDiagram::insert_line_site1(int idx1, int idx2) {
    // find the vertices corresponding to idx1 and idx2
    HEVertex start, end;
    bool start_found=false;
    bool end_found=false;
    BOOST_FOREACH( HEVertex v, g.vertices() ) {
        if ( g[v].index == idx1 ) {
            start = v;
            start_found = true;
        }
        if (g[v].index == idx2) {
            end = v;
            end_found = true;
        }
    }
    assert(start_found);
    assert(end_found);
    std::cout << " found startvert = " << start << " " << g[start].position <<"\n";
    std::cout << "   found endvert = " << end << " " << g[end].position << "\n";
    
    g[start].type=ENDPOINT; 
    g[start].status=OUT; 
    g[end].type=ENDPOINT; 
    g[end].status=OUT; 
    
    LineSite* line_site = new LineSite( g[start].position, g[end].position) ;
    // seed-face is face of start-point
    HEFace closest_face = g[start].face; 
    // seed 
    HEVertex v_seed = find_seed_vertex(closest_face, line_site ) ;
    g[v_seed].status = IN;
    
    std::cout << "   seed  = " << v_seed << " " << g[v_seed].position << "\n";
}

void VoronoiDiagram::insert_line_site2(int idx1, int idx2) {
    // find the vertices corresponding to idx1 and idx2
    HEVertex start, end;
    bool start_found=false;
    bool end_found=false;
    BOOST_FOREACH( HEVertex v, g.vertices() ) {
        if ( g[v].index == idx1 ) {
            start = v;
            start_found = true;
        }
        if (g[v].index == idx2) {
            end = v;
            end_found = true;
        }
    }
    assert(start_found);
    assert(end_found);
    std::cout << " found startvert = " << start << " " << g[start].position <<"\n";
    std::cout << "   found endvert = " << end << " " << g[end].position << "\n";
    
    g[start].type=ENDPOINT; 
    g[start].status=OUT; 
    g[end].type=ENDPOINT; 
    g[end].status=OUT; 
    
    LineSite* line_site = new LineSite( g[start].position, g[end].position) ;
    // seed-face is face of start-point
    HEFace closest_face = g[start].face; 
    // seed 
    HEVertex v_seed = find_seed_vertex(closest_face, line_site ) ;
    std::cout << "   seed  = " << v_seed << " " << g[v_seed].position << "\n";
    augment_vertex_set(v_seed, line_site );
    std::cout << "   after augment: v0.size() = " << v0.size() << "\n";
    
    
}

void VoronoiDiagram::insert_line_site3(int idx1, int idx2) {
    // find the vertices corresponding to idx1 and idx2
    HEVertex start, end;
    bool start_found=false;
    bool end_found=false;
    BOOST_FOREACH( HEVertex v, g.vertices() ) {
        if ( g[v].index == idx1 ) {
            start = v;
            start_found = true;
        }
        if (g[v].index == idx2) {
            end = v;
            end_found = true;
        }
    }
    assert(start_found);
    assert(end_found);
    std::cout << " found startvert = " << start << " " << g[start].position <<"\n";
    std::cout << "   found endvert = " << end << " " << g[end].position << "\n";
    
    g[start].type=ENDPOINT; 
    g[start].status=OUT; 
    g[end].type=ENDPOINT; 
    g[end].status=OUT; 
    
    LineSite* line_site = new LineSite( g[start].position, g[end].position) ;
    // seed-face is face of start-point
    HEFace closest_face = g[start].face; 
    // seed 
    HEVertex v_seed = find_seed_vertex(closest_face, line_site ) ;
    augment_vertex_set(v_seed, line_site );
    std::cout << "   after augment: v0.size() = " << v0.size() << "\n";
    add_new_vertices( line_site ); 
    
}

void VoronoiDiagram::insert_line_site(int idx1, int idx2) {
    // find the vertices corresponding to idx1 and idx2
    HEVertex start, end;
    bool start_found=false;
    bool end_found=false;
    BOOST_FOREACH( HEVertex v, g.vertices() ) {
        if ( g[v].index == idx1 ) {
            start = v;
            start_found = true;
        }
        if (g[v].index == idx2) {
            end = v;
            end_found = true;
        }
    }
    assert(start_found);
    assert(end_found);
    std::cout << " found startvert = " << start << " " << g[start].position <<"\n";
    std::cout << "   found endvert = " << end << " " << g[end].position << "\n";
    
    g[start].type=ENDPOINT; 
    g[start].status=OUT; 
    g[end].type=ENDPOINT; 
    g[end].status=OUT; 
    
    LineSite* line_site = new LineSite( g[start].position, g[end].position) ;
    
    HEEdge s_e = g.add_edge(start,end);
    HEEdge e_s = g.add_edge(end,start);
    g[s_e].type = LINESITE;
    g[e_s].type = LINESITE;
    g.twin_edges(s_e,e_s);
    
    // seed-face is face of start-point
    HEFace start_face = g[start].face; 
    HEFace end_face = g[end].face; 
    // seed 
    HEVertex v_seed = find_seed_vertex(start_face, line_site ) ;
    std::cout << "   seed  = " << v_seed << " " << g[v_seed].position << "\n";
    augment_vertex_set(v_seed, line_site );
    std::cout << "   after augment: v0.size() = " << v0.size() << "\n";
    add_new_vertices( line_site );  
    
    add_separator(start_face, line_site, start);
    add_separator(end_face, line_site, end);
    
    add_new_face( line_site );
    
    //HEFace newface = add_new_face( line_site );
    //remove_vertex_set( newface );
    
    //g[new_vert].face = newface;
    
    //reset_status();
    //assert( vd_checker->isValid() );
    return; // g[new_vert].index;
}

// add separator on the face of the endpoint of s
void VoronoiDiagram::add_separator(HEFace f, Site* s, HEVertex endp) {
    EdgeVector out_edges = g.out_edges(endp);
    assert( out_edges.size() == 1);
    HEEdge segment_e = out_edges[0];
    HEEdge segment_tw = g[segment_e].twin;
    /*
    VertexVector verts = g.face_vertices(f);
    VertexVector new_verts;
    BOOST_FOREACH(HEVertex v, verts) {
        if (g[v].status==NEW)
            new_verts.push_back(v);
    }
    assert(new_verts.size()==2);
    HEVertex v1=new_verts[0];
    HEVertex v2=new_verts[1];
    */
    HEVertex v1; // this Vertex is found as OUT-NEW-IN
    HEVertex v2; // this Vertex is found as IN-NEW-OUT
    HEEdge v2_previous, v1_next;
    HEEdge v2_next, v1_previous;
    //Site* f_site = g[f].site;
    //Site* new_site = g[newface].site;
                 // in-new    new   new-out
    boost::tie( v2_previous, v2,   v2_next) = find_new_vertex(f, OUT); // NEW->OUT vertex
    
                 // out-new    new   new-in
    boost::tie( v1_previous, v1, v1_next) = find_new_vertex(f, IN);  // NEW->IN  vertex
    
    // figure out which one is left and which is right.
    assert( s->in_region( g[v1].position ) );
    assert( s->in_region( g[v2].position ) );
    std::cout << " separator endpoint k3=" << g[v1].k3 << "\n";
    std::cout << " separator endpoint k3=" << g[v2].k3 << "\n";
    assert( g[v1].k3 !=   g[v2].k3 ); // should be on opposite sides
    /*
    HEVertex vpos,vneg;
    if ( g[v1].k3 == 1 ) {
        assert(g[v2].k3 == -1);
        vpos = v1;
        vneg = v2;
    } else {
        assert(g[v1].k3 == -1);
        assert(g[v2].k3 == +1);
        vpos = v2;
        vneg = v1;
    }*/
    
    HEEdge e2 = g.add_edge( endp, v2 );
    HEEdge e2_tw = g.add_edge( v2, endp );
    HEEdge e1 = g.add_edge( endp, v1 );
    HEEdge e1_tw = g.add_edge( v1, endp );
    // set these for new edges: next, face, twin, k
    
    // twin
    g.twin_edges(e1,e1_tw);
    g.twin_edges(e2,e2_tw);
    
    g[e2].type = SEPARATOR;
    g[e2_tw].type = SEPARATOR;
    g[e1].type = SEPARATOR;
    g[e1_tw].type = SEPARATOR;
    
    // face
    g[e2].face = f;
    g[e1_tw].face = f;
    g[f].edge = e2;
    
    // k
    g[e2].k    = g[v2].k3;
    g[e2_tw].k = g[v2].k3;
    g[e1].k    = g[v1].k3;
    g[e1_tw].k = g[v1].k3;
    
    // next-pointers for endpoint-face
    g[v1_previous].next = e1_tw;
    g[e1_tw].next = e2;
    g[e2].next = v2_next;
    
    // next-pointers for segment face(s)
    g[v2_previous].next = e2_tw;
    g[e2_tw].next = segment_e;
    g[segment_tw].next = e1;
    g[e1].next = v1_next;
    
    g[f].status = NONINCIDENT; // face is "done"
    assert( vd_checker->face_ok( f ) );
}

// find amount of clearance-disk violation on all face vertices and return vertex with the largest violation
HEVertex VoronoiDiagram::find_seed_vertex(HEFace f, Site* site) const { //const Point& p) {
    VertexVector face_verts = g.face_vertices(f);
    assert( face_verts.size() >= 3 );
    double minPred(0.0); 
    HEVertex minimalVertex =  HEVertex() ;
    bool first = true;
    BOOST_FOREACH( HEVertex q, face_verts) { // go thorugh all the vertices and find the one with smallest detH
        if ( g[q].status != OUT ) {
            double h = g[q].in_circle( site->apex_point( g[q].position ) ); 
            if ( first || ( (h<minPred) && (site->in_region(g[q].position) ) ) ) {
                minPred = h;
                minimalVertex = q;
                first = false;
            }
        }
    }
    assert( minPred < 0 );
    // FIXME not using Point p anymore: assert( vd_checker->inCircle_is_negative( p, f, minimalVertex ) );
    return minimalVertex;
}

// growing the v0/delete-tree of "IN" vertices by "weighted breadth-first search"
// we start at the seed and add vertices with detH<0 provided that:
// (C4) v should not be adjacent to two or more IN vertices (this would result in a loop/cycle!)
// (C5) for an incident face containing v: v is adjacent to an IN vertex on this face
// C4 and C5 refer to the Sugihara&Iri 1992 "one million" paper 
//  we process UNDECIDED vertices adjacent to known IN-vertices in a "weighted breadth-first-search" manner
//  where vertices with a large fabs(detH) are processed first, since we assume the detH to be more reliable the larger fabs(detH) is.
void VoronoiDiagram::augment_vertex_set( HEVertex& v_seed, Site* site ) {
    mark_vertex( v_seed, site );
    modified_vertices.push_back( v_seed );
    while( !vertexQueue.empty() ) {
        HEVertex v;
        double h;
        boost::tie( v, h ) = vertexQueue.top();      assert( g.g[v].status == UNDECIDED );
        vertexQueue.pop(); 
        if ( h < 0.0 ) { // mark IN if detH<0 and passes (C4) and (C5) tests and in_region(). otherwise mark OUT
            if ( predicate_c4(v) || !predicate_c5(v) || !site->in_region(g[v].position) ) 
                g[v].status = OUT; // C4 or C5 violated, so mark OUT
            else {
                mark_vertex( v,  site); // h<0 and no violations, so mark IN. push adjacent UNDECIDED vertices onto Q.
            }
        } else { 
            g[v].status = OUT; // detH was positive (or zero), so mark OUT
        }
        modified_vertices.push_back( v );
    }
    // sanity-check: for all incident_faces the IN-vertices should be connected
    assert( vd_checker->incidentFaceVerticesConnected(  IN ) );
    assert( vd_checker->incidentFaceVerticesConnected(  OUT ) );
}

// mark vertex IN. mark adjacent faces INCIDENT
// push adjacent UNDECIDED vertices onto queue 
void VoronoiDiagram::mark_vertex(HEVertex& v,  Site* site) {
    g[v].status = IN;
    v0.push_back( v );
    mark_adjacent_faces( v );
    
    // also push the v-adjacent vertices onto the queue
    BOOST_FOREACH( HEVertex w, g.adjacent_vertices(v) ) {
        if ( (g[w].status == UNDECIDED) && (!g[w].in_queue) ) {
                // when pushing onto queue we also evaluate in_circle predicate so that we process vertices in the correct order
                vertexQueue.push( VertexDetPair(w , g[w].in_circle(site->apex_point(g[w].position)) ) ); 
                g[w].in_queue=true;
        }
    }
}

// IN-Vertex v has three adjacent faces, mark nonincident faces incident
// and push them to the incident_faces queue
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

// generate new voronoi-vertices on all IN-OUT edges 
void VoronoiDiagram::add_new_vertices( Site* new_site ) {
    assert( !v0.empty() );
    EdgeVector q_edges = find_in_out_edges();       // new vertices generated on these IN-OUT edges
    for( unsigned int m=0; m<q_edges.size(); ++m )  {   
        HEVertex q = g.add_vertex();
        g[q].status = NEW;
        modified_vertices.push_back(q);
        g[q].position = vpos->position( q_edges[m], new_site ); // set position
        g[q].k3 = vpos->get_k3();
        g[q].init_dist( new_site->apex_point( g[q].position ) ); // set initial clearance-disk
        add_vertex_in_edge(q, q_edges[m] );
    }
}

void VoronoiDiagram::add_vertex_in_edge( HEVertex v, HEEdge e) {
    // the vertex v is inserted into the middle of edge e
    // edge e and its twin are replaced by four new edges: e1,e2 and their twins te2,te1
    // 
    //                    face
    //                    e1   e2
    // previous-> source  -> v -> target -> next
    //            tw_trg  <- v <- tw_src <- tw_previous
    //                    te2  te1
    //                    twin_face
    HEEdge twin = g[e].twin;
    HEVertex source = g.source(e); //boost::source( e , g );
    HEVertex target = g.target(e); //boost::target( e , g);
    HEVertex twin_source = g.source(twin); //boost::source( twin , g);
    HEVertex twin_target = g.target(twin); //boost::target( twin , g );
    assert( source == twin_target );    
    assert( target == twin_source );
    HEFace face = g[e].face;
    HEFace twin_face = g[twin].face;
    HEEdge previous = g.previous_edge(e);
    assert( g[previous].face == g[e].face );
    HEEdge twin_previous = g.previous_edge(twin);
    assert( g[twin_previous].face == g[twin].face );
    
    HEEdge e1 = g.add_edge( source, v ); // these replace e
    HEEdge e2 = g.add_edge( v, target );    
    // preserve the left/right face link
    g[e1].face = face;
    g[e2].face = face;
    // next-pointers
    g[previous].next = e1;
    g[e1].next = e2;
    g[e2].next = g[e].next;
    // k-vals
    g[e1].k = g[e].k;
    g[e2].k = g[e].k;
    // type
    g[e1].type = g[e].type;
    g[e2].type = g[e].type;
    
    HEEdge te1 = g.add_edge( twin_source, v  ); // these replace twin
    HEEdge te2 = g.add_edge( v, twin_target  );
    
    g[te1].face = twin_face;
    g[te2].face = twin_face;
    
    g[twin_previous].next = te1;
    g[te1].next = te2;
    g[te2].next = g[twin].next;
    
    // TWINNING (note indices 'cross', see ASCII art above)
    g[e1].twin = te2;
    g[te2].twin = e1;
    g[e2].twin = te1;
    g[te1].twin = e2;
    // k-vals
    g[te1].k = g[twin].k;
    g[te2].k = g[twin].k;
    // type
    g[te1].type = g[twin].type;
    g[te2].type = g[twin].type;
    
    // update the faces (required here?)
    g[face].edge = e1;
    g[twin_face].edge = te1;
    
    // finally, remove the old edge
    g.remove_edge(e);
    g.remove_edge(twin);
}

// add a new face corresponding to the new Site
// call add_new_edge() on all the incident_faces that should be split
HEFace VoronoiDiagram::add_new_face(Site* s) { 
    HEFace newface =  g.add_face(); 
    g[newface].site = s;
    g[newface].status = NONINCIDENT;
    if (s->isPoint() )
        fgrid->add_face( g[newface] ); 
    BOOST_FOREACH( HEFace f, incident_faces ) {
        if ( g[f].status == INCIDENT ) // end-point faces already dealt with in add_separator()
            add_new_edge(newface, f); // each INCIDENT face is split into two parts: newface and f
    }
    return newface;
}

// by adding a NEW-NEW edge, split the face f into one part which is newface, and the other part is the old f
void VoronoiDiagram::add_new_edge(HEFace newface, HEFace f) {
    HEVertex new_source; // this Vertex is found as OUT-NEW-IN
    HEVertex new_target; // this Vertex is found as IN-NEW-OUT
    HEEdge new_previous, new_next, twin_next, twin_previous;
    Site* f_site = g[f].site;
    Site* new_site = g[newface].site;
    
    boost::tie( new_previous, new_source, twin_next) = find_new_vertex(f, OUT); // NEW->OUT vertex
    boost::tie( twin_previous, new_target, new_next) = find_new_vertex(f, IN);  // NEW->IN  vertex
    // now connect:   new_previous -> new_source -> new_target -> new_next
    // and:              twin_next <- new_source <- new_target <- twin_previous    
    HEEdge e_new = g.add_edge( new_source, new_target );
    g[e_new].set_parameters( f_site, new_site );
    g[e_new].next = new_next;
    g[e_new].k = g[new_next].k; // the next edge is on the same face, so has the correct k-value
    g[e_new].face = f;
    g[new_previous].next = e_new;
    g[f].edge = e_new; 
    // the twin edge that bounds the new face
    HEEdge e_twin = g.add_edge( new_target, new_source );
    
    g[e_twin].set_parameters( new_site, f_site );
    std::cout << g[e_new].type << "\n";
    
    g[e_twin].next = twin_next; 
    
    // the NEW vertice were generated with vpos, which has set the offset direction k3
    /*
    if ( g[new_target].k3 != g[new_source].k3 ) {
        std::cout << " Error: g[new_target].k3= " << g[new_target].k3 << " g[new_source].k3 " << g[new_source].k3 << "\n";
        std::cout << " Error: g[new_target].position= " << g[new_target].position << " g[new_source].position " << g[new_source].position << "\n";
    }
    
    */
    assert( g[new_target].k3 == g[new_source].k3 );
    g[e_twin].k = g[new_source].k3; //g[twin_next].k;
    g[e_twin].face = newface;
    g[twin_previous].next = e_twin;
    g[newface].edge = e_twin; 
    g.twin_edges(e_new,e_twin);
}

// among the edges of HEFace f, which has had NEW vertices inserted into two (?exactly two?) of its edges,
// find the NEW vertex with status-signature NEW->s.
// return a tuple with the previous_edge, found_vertex, and the twin_next edge
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


// this is called after split_face()
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

// at the end after an incremental insertion of a new site,
// reset status of modified_vertices to UNDECIDED and incident_faces to NONINCIDENT,
// so that we are ready for the next insertion.
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

// given the set v0 of "IN" vertices, find and return the adjacent IN-OUT edges 
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

// number of IN vertices adjacent to given vertex v
// predicate C4 i.e. "adjacent in-count" from Sugihara&Iri 1992 "one million" paper
bool VoronoiDiagram::predicate_c4(HEVertex v) {
    int in_count=0;
    BOOST_FOREACH( HEVertex w, g.adjacent_vertices(v) ) {
        if ( g[w].status == IN )
            in_count++;
    }
    return (in_count >= 2);
}

// do any of the three faces that are adjacent to the given IN-vertex v have an IN-vertex ?
// predicate C5 i.e. "connectedness"  from Sugihara&Iri 1992 "one million" paper
bool VoronoiDiagram::predicate_c5(HEVertex v) {
    FaceVector adj_faces = g.adjacent_faces(v);   
    assert( adj_faces.size() == 3 );
    FaceVector adjacent_incident_faces; // find the ajacent incident faces
    BOOST_FOREACH( HEFace f, adj_faces ) {
        if ( g[f].status == INCIDENT )
            adjacent_incident_faces.push_back( f );
    }
    assert( !adjacent_incident_faces.empty() );
    
    bool all_found = true;
    BOOST_FOREACH( HEFace f, adjacent_incident_faces ) { // check each adjacent face f for an IN-vertex
        bool face_ok=false;
        BOOST_FOREACH( HEVertex w, g.face_vertices(f) ) { 
            if ( w != v && g[w].status == IN && g.has_edge(w,v) )  // v should be adjacent to an IN vertex on the face
                face_ok = true;
        }
        if (!face_ok)
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

std::string VoronoiDiagram::print() const {
    std::ostringstream o;
    o << "VoronoiDiagram (nVerts="<< g.num_vertices() << " , nEdges="<< g.num_edges() <<"\n";
    return o.str();
}

} // end namespace
// end file voronoidiagram.cpp
