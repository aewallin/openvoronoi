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

#include <cassert>

#include <boost/foreach.hpp>
#include <boost/math/tools/roots.hpp> // for toms748



#include "voronoidiagram.hpp"
#include "facegrid.hpp"
#include "checker.hpp"

namespace ovd {

VoronoiDiagram::VoronoiDiagram(double far, unsigned int n_bins) {
    fgrid = new FaceGrid(far, n_bins);
    vd_checker = new VoronoiDiagramChecker(this); // helper-class that checks topology/geometry
    vpos = new VertexPositioner(this); // helper-class that positions vertices
    far_radius=far;
    initialize();
    num_psites=3;
    num_lsites=0;
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
    g[g1] = VoronoiVertex( gen1 , OUT, POINTSITE);
    g[g2] = VoronoiVertex( gen2 , OUT, POINTSITE);
    g[g3] = VoronoiVertex( gen3 , OUT, POINTSITE);

    // APEX 1
    Point apex1 = 0.5*(gen2+gen3);
    HEVertex a1 = g.add_vertex();
    g[a1] = VoronoiVertex( apex1, UNDECIDED, APEX );
    g[a1].init_dist(gen2);
    
    // APEX 2
    Point apex2 = 0.5*(gen1+gen3);
    HEVertex a2 = g.add_vertex();
    g[a2] = VoronoiVertex( apex2, UNDECIDED, APEX );
    g[a2].init_dist(gen3);   
    
    // APEX 3
    Point apex3 = 0.5*(gen1+gen2);
    HEVertex a3 = g.add_vertex();
    g[a3] = VoronoiVertex( apex3, UNDECIDED, APEX );
    g[a3].init_dist(gen1);   

    // add face 1: v0-v1-v2 which encloses gen3
    HEEdge e1_1 =  g.add_edge( v0 , a1 ); //v01 );   
    HEEdge e1_2 =  g.add_edge( a1 , v01); //v01 );  
    HEEdge e2 =  g.add_edge( v01, v02 );
    HEEdge e3_1 =  g.add_edge( v02, a2 ); 
    HEEdge e3_2 =  g.add_edge( a2 , v0 ); 
    HEFace f1 =  g.add_face();
    g[f1].edge = e2;
    g[f1].site = new PointSite(gen3,f1);
    g[f1].status = NONINCIDENT;
    
    fgrid->add_face( g[f1] );
    g[e1_1].face = f1;
    g[e1_2].face = f1;
    g[e2].face = f1;
    g[e3_1].face = f1;
    g[e3_2].face = f1;
    g[e1_1].next = e1_2;
    g[e1_2].next = e2;
    g[e2].next = e3_1;
    g[e3_1].next = e3_2;
    g[e3_2].next = e1_1;

    // add face 2: v0-v02-v03 which encloses gen1
    HEEdge e4_1 = g.add_edge( v0, a2  );
    HEEdge e4_2 = g.add_edge( a2, v02 );

    HEEdge e5 = g.add_edge( v02, v03  );
    HEEdge e6_1 = g.add_edge( v03, a3 );
    HEEdge e6_2 = g.add_edge( a3, v0 ); 
    HEFace f2 =  g.add_face();
    g[f2].edge = e5;
    g[f2].site = new PointSite(gen1,f2);
    g[f2].status = NONINCIDENT;
    
    fgrid->add_face( g[f2] );
    g[e4_1].face = f2;
    g[e4_2].face = f2;
    g[e5].face = f2;
    g[e6_1].face = f2;
    g[e6_2].face = f2;
    g[e4_1].next = e4_2;
    g[e4_2].next = e5;
    g[e5].next = e6_1;
    g[e6_1].next = e6_2;
    g[e6_2].next = e4_1;
    
    // add face 3: v0-v3-v1 which encloses gen2
    HEEdge e7_1 = g.add_edge( v0 , a3 );  
    HEEdge e7_2 = g.add_edge( a3 , v03 );   
    HEEdge e8 = g.add_edge( v03, v01 );
    HEEdge e9_1 = g.add_edge( v01, a1  ); 
    HEEdge e9_2 = g.add_edge( a1 , v0  ); 

    HEFace f3 =  g.add_face();
    g[f3].edge = e8;
    g[f3].site = new PointSite(gen2,f3);
    g[f3].status = NONINCIDENT;
    fgrid->add_face( g[f3] );
    g[e7_1].face = f3;
    g[e7_2].face = f3;
    g[e8].face = f3;
    g[e9_1].face = f3;
    g[e9_2].face = f3;
    g[e7_1].next = e7_2;
    g[e7_2].next = e8;
    g[e8].next = e9_1;
    g[e9_1].next = e9_2;
    g[e9_2].next = e7_1;

    // set type. (note that edge-params x[8] and y[8] are not set!
    g[e1_1].type = LINE;  g[e1_1].set_parameters(g[f1].site, g[f3].site, false);
    g[e1_2].type = LINE;  g[e1_2].set_parameters(g[f1].site, g[f3].site, true);
    g[e2].type = OUTEDGE; 
    g[e3_1].type = LINE; g[e3_1].set_parameters(g[f2].site, g[f1].site, true);
    g[e3_2].type = LINE; g[e3_2].set_parameters(g[f2].site, g[f1].site, false);
    g[e4_1].type = LINE; g[e4_1].set_parameters(g[f2].site, g[f1].site, false);
    g[e4_2].type = LINE; g[e4_2].set_parameters(g[f2].site, g[f1].site, true);
    g[e5].type = OUTEDGE;
    g[e6_1].type = LINE; g[e6_1].set_parameters(g[f2].site, g[f3].site, false);
    g[e6_2].type = LINE; g[e6_2].set_parameters(g[f2].site, g[f3].site, true);
    g[e7_1].type = LINE; g[e7_1].set_parameters(g[f2].site, g[f3].site, true);
    g[e7_2].type = LINE; g[e7_2].set_parameters(g[f2].site, g[f3].site, false);
    g[e8].type = OUTEDGE;
    g[e9_1].type = LINE; g[e9_1].set_parameters(g[f1].site, g[f3].site, true);
    g[e9_2].type = LINE; g[e9_2].set_parameters(g[f1].site, g[f3].site, false);
    
    // twin edges
    g[e1_1].twin = e9_2;
    g[e1_2].twin = e9_1;
    g[e9_1].twin = e1_2;
    g[e9_2].twin = e1_1;

    g[e2].twin = HEEdge(); // the outermost edges have invalid twins
    g[e5].twin = HEEdge();
    g[e8].twin = HEEdge();
    g[e3_1].twin = e4_2;
    g[e3_2].twin = e4_1;
    g[e4_1].twin = e3_2;
    g[e4_2].twin = e3_1;
    g[e6_1].twin = e7_2;
    g[e6_2].twin = e7_1;
    g[e7_1].twin = e6_2;
    g[e7_2].twin = e6_1;

    // k-values all positive for PointSite generators
    g[e1_1].k = 1.0;
    g[e1_2].k = 1.0;
    g[e2].k = 1.0;
    g[e3_1].k = 1.0;
    g[e3_2].k = 1.0;
    g[e4_1].k = 1.0;
    g[e4_2].k = 1.0;
    g[e5].k = 1.0;
    g[e6_1].k = 1.0;
    g[e6_2].k = 1.0;
    g[e7_1].k = 1.0;
    g[e7_2].k = 1.0;
    g[e8].k = 1.0;
    g[e9_1].k = 1.0;
    g[e9_2].k = 1.0;

    assert( vd_checker->is_valid() );
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
    num_psites++;
    assert( p.norm() < far_radius );     // only add vertices within the far_radius circle
    
    HEVertex new_vert = g.add_vertex();
    g[new_vert].position=p; 
    g[new_vert].type=POINTSITE; 
    g[new_vert].status=OUT; 
    g[new_vert].site = new PointSite(p);
    
    HEFace closest_face = fgrid->grid_find_closest_face( p ); 
    HEVertex v_seed = find_seed_vertex(closest_face, g[new_vert].site ) ;
    augment_vertex_set(v_seed, g[new_vert].site );
    add_vertices( g[new_vert].site );    
    HEFace newface = add_face( g[new_vert].site );
      
    BOOST_FOREACH( HEFace f, incident_faces ) {
        add_edge(newface, f); // no newface2 parameter given!
    }
    repair_face( newface );
    remove_vertex_set();
    reset_status();
    
    assert( vd_checker->face_ok( newface ) );
    assert( vd_checker->is_valid() );
    print_faces();
    return g[new_vert].index;
}

bool VoronoiDiagram::insert_line_site_step(int idx1, int idx2, int step) {
    num_lsites++;
    // find the vertices corresponding to idx1 and idx2
    HEVertex start, end;
    bool start_found=false;
    bool end_found=false;
    // FIXME: search for vertex-descriptor from an index/descriptor table ( O(1) ?)
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
    std::cout << " found startvert = " << g[start].index << " " << g[start].position <<"\n";
    std::cout << "   found endvert = " << g[end].index << " " << g[end].position << "\n";
    
    g[start].type=ENDPOINT; 
    g[start].status=OUT; 
    g[end].type=ENDPOINT; 
    g[end].status=OUT; 
    
    // create a point which is left of src->trg
    // determine k (offset-dir) for this point
    // the we know which site/face is the k==+1 and which is k==-1
    Point src_se = g[start].position;
    Point trg_se = g[end  ].position;
    Point left = 0.5*(src_se+trg_se) + (trg_se-src_se).xyPerp();
    LineSite* pos_site;
    LineSite* neg_site;
    bool se_sign = left.is_right(src_se,trg_se);
    HEEdge pos_edge, neg_edge;
    if ( se_sign ) {
        pos_site = new LineSite( g[start].position, g[end  ].position, +1 );
        neg_site = new LineSite( g[end  ].position, g[start].position, -1 );
        pos_edge = g.add_edge(start,end);
        neg_edge = g.add_edge(end,start);
    } else {
        pos_site = new LineSite( g[end  ].position, g[start].position, +1 );
        neg_site = new LineSite( g[start].position, g[end  ].position, -1 );
        pos_edge = g.add_edge(end,start);
        neg_edge = g.add_edge(start,end);
    }
    g[pos_edge].type = LINESITE;
    g[neg_edge].type = LINESITE;
    g[pos_edge].k = +1;
    g[neg_edge].k = -1;
    g.twin_edges(pos_edge,neg_edge);
    
    // seed-face is face of start-point
    HEFace start_face = g[start].site->face; 
    HEFace   end_face = g[end  ].site->face;
    
    // seed 
    HEVertex v_seed = find_seed_vertex(start_face, pos_site ) ;
    std::cout << " start face seed  = " << g[v_seed].index << " " << g[v_seed].position << "\n";
    g[v_seed].status = IN;
    
    if (step==1) return false;
        
    augment_vertex_set(v_seed, pos_site ); // should not matter if we use pos_site or neg_site here
    std::cout << "   after augment: v0.size() = " << v0.size() << "\n";
    std::cout << "   delete-set is: "; print_vertices(v0);    
    
    if (step==2) return false;

    // check that end_face is INCIDENT?
    // check that tree includes end_face_seed ?
    
    int n_verts=0;
    assert( !v0.empty() );
    EdgeVector q_edges = find_in_out_edges();       // new vertices generated on these IN-OUT edges
    print_edges(q_edges);
    
    for( unsigned int m=0; m<q_edges.size(); ++m )  {   
        HEVertex q = g.add_vertex();
        g[q].status = NEW;
        modified_vertices.push_back(q);
        std::cout << "position new vertex " << g[q].index << " on " <<  g[ g.source(q_edges[m])].index << "(t=" << g[ g.source(q_edges[m])].dist() << ")-"; 
        std::cout << g[ g.target(q_edges[m])].index << "(t=" << g[ g.target(q_edges[m])].dist() <<  " edge, type=" << g[q_edges[m]].type << "\n";

        g[q].position = vpos->position( q_edges[m], pos_site ); // set position
        g[q].k3 = vpos->get_k3();
        g[q].init_dist( pos_site->apex_point( g[q].position ) ); // set initial clearance-disk
        add_vertex_in_edge(q, q_edges[m] );
        n_verts++;
        if (step == (2+n_verts)) return false;
    }
    
    
    HEFace pos_face = add_face( pos_site ); //  this face to the left of start->end edge    
    HEFace neg_face = add_face( neg_site ); //  this face is to the left of end->start edge
    g[pos_face].edge = pos_edge;
    g[neg_face].edge = neg_edge;
    g[pos_edge].face = pos_face;
    g[neg_edge].face = neg_face;
    
    std::cout << " add_separator(" << g[start].index << ") \n";
    add_separator(start_face, start, pos_site, neg_site);
    if (step == (2+n_verts+1)) return false;
    
    std::cout << " add_separator(" << g[end].index << ") \n";
    add_separator(end_face  , end  , pos_site, neg_site); 
    if (step == (2+n_verts+2)) return false;

    assert( vd_checker->face_ok( start_face ) );
    assert( vd_checker->face_ok( end_face ) );
    //print_faces();
    std::cout << "start face: ";
    print_face(start_face);
    std::cout << "end face: ";
    print_face(end_face);

    BOOST_FOREACH( HEFace f, incident_faces ) {
        if ( g[f].status == INCIDENT )  {// end-point faces already dealt with in add_separator()
            add_edge(pos_site->face, f, neg_site->face); // each INCIDENT face is split into two parts: newface and f
        }
    }
    if (step == (2+n_verts+3)) return false;

    std::cout << "new edges added \n";
    remove_vertex_set();
    std::cout << " k=+1 face is " << pos_site->face << "\n";
    std::cout << " k=-1 face is " << neg_site->face << "\n";
    
    repair_face( pos_face );
    repair_face( neg_face );
    
    std::cout << " pos face: "; print_face( pos_face );
    std::cout << " neg face: "; print_face( neg_face );
    std::cout << "faces " << start_face << " " << end_face << " " << pos_face << " " << neg_face << " repaired \n";

    reset_status();
    std::cout << "insert_line_site(" << g[start].index << "-"<< g[end].index << ") done.\n";
    assert( vd_checker->face_ok( pos_face ) );
    assert( vd_checker->face_ok( neg_face ) );
    assert( vd_checker->is_valid() );
    print_faces();
    return true;
}
    
void VoronoiDiagram::insert_line_site(int idx1, int idx2) {
    num_lsites++;
    // find the vertices corresponding to idx1 and idx2
    HEVertex start, end;
    bool start_found=false;
    bool end_found=false;
    // FIXME: search for vertex-descriptor from an index/descriptor table ( O(1) ?)
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
    std::cout << " found startvert = " << g[start].index << " " << g[start].position <<"\n";
    std::cout << "   found endvert = " << g[end].index << " " << g[end].position << "\n";
    
    g[start].type=ENDPOINT; 
    g[start].status=OUT; 
    g[end].type=ENDPOINT; 
    g[end].status=OUT; 
    
    // create a point which is left of src->trg
    // determine k (offset-dir) for this point
    // the we know which site/face is the k==+1 and which is k==-1
    Point src_se = g[start].position;
    Point trg_se = g[end  ].position;
    Point left = 0.5*(src_se+trg_se) + (trg_se-src_se).xyPerp();
    LineSite* pos_site;
    LineSite* neg_site;
    bool se_sign = left.is_right(src_se,trg_se);
    HEEdge pos_edge, neg_edge;
    if ( se_sign ) {
        pos_site = new LineSite( g[start].position, g[end  ].position , +1);
        neg_site = new LineSite( g[end  ].position, g[start].position , -1);
        pos_edge = g.add_edge(start,end);
        neg_edge = g.add_edge(end,start);
    } else {
        pos_site = new LineSite( g[end  ].position, g[start].position , +1);
        neg_site = new LineSite( g[start].position, g[end  ].position , -1);
        pos_edge = g.add_edge(end,start);
        neg_edge = g.add_edge(start,end);
    }
    g[pos_edge].type = LINESITE;
    g[neg_edge].type = LINESITE;
    g[pos_edge].k = +1;
    g[neg_edge].k = -1;
    g.twin_edges(pos_edge,neg_edge);
    
    // seed-face is face of start-point
    HEFace start_face = g[start].site->face; 
    HEFace   end_face = g[end  ].site->face;
    
    // seed 
    HEVertex v_seed = find_seed_vertex(start_face, pos_site ) ;
    std::cout << " start face seed  = " << g[v_seed].index << " " << g[v_seed].position << "\n";
    g[v_seed].status = IN;
        
    augment_vertex_set(v_seed, pos_site ); // should not matter if we use pos_site or neg_site here
    std::cout << "   delete-set is("<< v0.size() <<"): "; print_vertices(v0);
    
    // check that end_face is INCIDENT?
    // check that tree includes end_face_seed ?
    
    add_vertices( pos_site );  
    
    HEFace pos_face = add_face( pos_site ); //  this face to the left of start->end edge    
    HEFace neg_face = add_face( neg_site ); //  this face is to the left of end->start edge
    g[pos_face].edge = pos_edge;
    g[neg_face].edge = neg_edge;
    g[pos_edge].face = pos_face;
    g[neg_edge].face = neg_face;
    
    add_separator(start_face, start, pos_site, neg_site);
    add_separator(end_face  , end  , pos_site, neg_site); 

    assert( vd_checker->face_ok( start_face ) );
    assert( vd_checker->face_ok( end_face ) );
    //print_faces();
    std::cout << "start face: ";
    print_face(start_face);
    std::cout << "end face: ";
    print_face(end_face);

    BOOST_FOREACH( HEFace f, incident_faces ) {
        if ( g[f].status == INCIDENT )  {// end-point faces already dealt with in add_separator()
            add_edge(pos_site->face, f, neg_site->face); // each INCIDENT face is split into two parts: newface and f
        }
    }

    std::cout << "new edges added \n";
    remove_vertex_set();
    std::cout << " k=+1 face is " << pos_site->face << "\n";
    std::cout << " k=-1 face is " << neg_site->face << "\n";
    
    repair_face( pos_face );
    repair_face( neg_face );
    
    std::cout << " pos face: "; print_face( pos_face );
    std::cout << " neg face: "; print_face( neg_face );
    std::cout << "faces " << start_face << " " << end_face << " " << pos_face << " " << neg_face << " repaired \n";

    reset_status();
    std::cout << "insert_line_site(" << g[start].index << "-"<< g[end].index << ") done.\n";
    assert( vd_checker->face_ok( pos_face ) );
    assert( vd_checker->face_ok( neg_face ) );
    assert( vd_checker->is_valid() );
    print_faces();
    return; 
}

// add separator on the face f, which contains the endpoint
// 
void VoronoiDiagram::add_separator(HEFace f, HEVertex endp, Site* s1, Site* s2) {
    //std::cout << "add_separator()\n";
    EdgeVector out_edges = g.out_edges(endp);
    assert( out_edges.size() == 1);
    HEEdge segment_e = out_edges[0];
    HEEdge segment_tw = g[segment_e].twin;
    
    EdgeData ed = find_edge_data(f);
    HEVertex v1 = ed.v1; // this Vertex is found as OUT->NEW->IN
    HEVertex v2 = ed.v2; // this Vertex is found as IN->NEW->OUT
    HEEdge v2_previous = ed.v2_prv;
    HEEdge v1_next = ed.v1_nxt;
    HEEdge v2_next = ed.v2_nxt;
    HEEdge v1_previous = ed.v1_prv;
    
    if (g[v1].k3 ==  g[v2].k3) {
        std::cout << " g[" << g[v1].index << "].k3=" << g[v1].k3 << "  !=  g[" << g[v2].index << "].k3=" << g[v2].k3 << "\n";
        std::cout << "WARNING in add_separator(): NEW-vertices should be on opposite sides of new site!\n";
        std::cout << " g[" << g[v1].index << "] " << g[v1].position << "\n";
        std::cout << " g[" << g[v2].index << "] " << g[v2].position << "\n";
        std::cout << " identical? " << (g[v1].position==g[v2].position) << "\n";
        
        // since the vertices are positioned at the same position, force a k3 value.
        if ( (g[v1].position==g[v2].position) ) {
            g[v1].k3 = +1;
            g[v2].k3 = -1; // FIXME, check k-value of vertex adjacent to v1/v2
        } else {
            assert(0);
        }
    }
    assert( g[v1].k3 !=  g[v2].k3 ); // v1 and v2 should be on opposite sides
    
    assert( s1->in_region( g[v1].position ) ); // v1 and v2 should be in the region of the line-site
    assert( s1->in_region( g[v2].position ) );
    assert( s2->in_region( g[v1].position ) );
    assert( s2->in_region( g[v1].position ) );

    HEEdge e2 = g.add_edge( endp, v2 );
    HEEdge e2_tw = g.add_edge( v2, endp );
    HEEdge e1 = g.add_edge( endp, v1 );
    HEEdge e1_tw = g.add_edge( v1, endp );
    // set these for new edges: next, face, twin, k
    
    // twin
    g.twin_edges(e1,e1_tw);
    g.twin_edges(e2,e2_tw);
    // type
    g[e2].type = SEPARATOR;
    g[e2_tw].type = SEPARATOR;
    g[e1].type = SEPARATOR;
    g[e1_tw].type = SEPARATOR;

    // k, offset direction
    g[e2].k    = +1; // e2 is on the ENDPOINT side
    g[e2_tw].k = g[v2].k3; // e2 is on the segment side
    g[e1].k    = g[v1].k3; // e1 is on the segment side
    g[e1_tw].k = +1; // e1_tw is on the ENDPOINT side

    // old endpoint face
    g[e2].face = f;
    g[e1_tw].face = f;
    g[f].edge = e2;

    // new faces:
    if (g[e1].k == 1) {
        g[e1].face = s1->face; // s1 is the k==1 face
        g[s1->face].edge=e1;
        g[e2_tw].face = s2->face;
        g[s2->face].edge=e2_tw;
    } else {
        g[e1].face = s2->face; // s2 is the k==-1 face
        g[s2->face].edge=e1;
        g[e2_tw].face = s1->face;
        g[s1->face].edge=e2_tw;
    }
    
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
    assert( vd_checker->face_ok( f ) ); // check that the old face is OK
}

// find amount of clearance-disk violation on all face vertices and return vertex with the largest violation
HEVertex VoronoiDiagram::find_seed_vertex(HEFace f, Site* site) const { //const Point& p) {
    VertexVector face_verts = g.face_vertices(f);
    assert( face_verts.size() >= 3 );
    double minPred(0.0); 
    HEVertex minimalVertex =  HEVertex() ;
    bool first = true;
    BOOST_FOREACH( HEVertex q, face_verts) { // go thorugh all the vertices and find the one with smallest detH
        if ( (g[q].status != OUT) && (g[q].type == NORMAL) ) {
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
//  where vertices with a large fabs(detH) are processed first, since we assume the in-circle predicate
//  to be more reliable the larger fabs(in_circle()) is.
void VoronoiDiagram::augment_vertex_set( HEVertex& v_seed, Site* site ) {
    mark_vertex( v_seed, site );
    modified_vertices.push_back( v_seed );
    while( !vertexQueue.empty() ) {
        HEVertex v;
        double h;
        boost::tie( v, h ) = vertexQueue.top();      assert( g.g[v].status == UNDECIDED );
        vertexQueue.pop(); 
        if ( h < 0.0 ) { // mark IN if h<0 and passes (C4) and (C5) tests and in_region(). otherwise mark OUT
            if ( predicate_c4(v) || !predicate_c5(v) || !site->in_region(g[v].position) ) {
                g[v].status = OUT; // C4 or C5 violated, so mark OUT
                std::cout << g[v].index << " marked OUT (topo): c4="<< predicate_c4(v) << " c5=" << !predicate_c5(v) << " r=" << !site->in_region(g[v].position) << " h=" << h << "\n";
            } else {
                mark_vertex( v,  site); // h<0 and no violations, so mark IN. push adjacent UNDECIDED vertices onto Q.
                std::cout << g[v].index << " marked IN (in_circle) ( " << h << " )\n";
            }
        } else { 
            g[v].status = OUT; // detH was positive (or zero), so mark OUT
            std::cout << g[v].index << " marked OUT (in_circle) ( " << h << " )\n";
        }
        modified_vertices.push_back( v );
    }
    
    // sanity-check: for all incident_faces the IN-vertices should be connected
    //assert( vd_checker->incidentFaceVerticesConnected(  IN ) );
    //assert( vd_checker->incidentFaceVerticesConnected(  OUT ) );
}

// mark vertex IN. mark adjacent faces INCIDENT
// push adjacent UNDECIDED vertices onto queue 
void VoronoiDiagram::mark_vertex(HEVertex& v,  Site* site) {
    g[v].status = IN;
    v0.push_back( v );
    mark_adjacent_faces( v, site );
    
    // push the v-adjacent vertices onto the queue
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
void VoronoiDiagram::mark_adjacent_faces( HEVertex v, Site* site) {
    assert( g[v].status == IN );
    FaceVector new_adjacent_faces = g.adjacent_faces( v ); 
        
    if (g[v].type == APEX)
        assert( new_adjacent_faces.size()==2 );
    else
        assert( new_adjacent_faces.size()==3 );
    
    BOOST_FOREACH( HEFace adj_face, new_adjacent_faces ) {
        if ( g[adj_face].status  != INCIDENT ) {
            if ( site->isLine() )
                add_split_vertex(adj_face, site);
                
            g[adj_face].status = INCIDENT; 
            incident_faces.push_back(adj_face);
        }
    }
}

// walk around the face f
// return edges whose endpoints are on separate sides of pt1-pt2 line
EdgeVector VoronoiDiagram::find_split_edges(HEFace f, Point pt1, Point pt2) {
    EdgeVector out;
    //HEVertex split_src, split_trg;
    HEEdge current_edge = g[f].edge;
    HEEdge start_edge = current_edge;
    //HEEdge split_edge;
    bool done = false;
    int count=0;                             
    while (!done) { // FIND ALL! not just one.
        HEVertex trg = g.target( current_edge );
        HEVertex src = g.source( current_edge );
        bool src_sign = g[src].position.is_right(pt1,pt2);
        bool trg_sign = g[trg].position.is_right(pt1,pt2);
        //std::cout << g[src].index << ":" << src_sign << " " << g[trg].index << ":" << trg_sign << "\n";
        if ( g[src].type == NORMAL || g[src].type == APEX ) {
            if ( src_sign != trg_sign ) {
                    //split_src = src;
                    //split_trg = trg;
                    //found = true;
                    //split_edge = current_edge;    
                    out.push_back(current_edge);             
            }
        }
        current_edge = g[current_edge].next;   
        count++;
        assert(count<100); // some reasonable max number of edges in face, to avoid infinite loop
        if ( current_edge == start_edge )
            done = true;
    }
    return out;
}

// add one or many split-vertices to the edges of the give face
//
// these are projections/mirrors of the site of f with the new Site s acting as the mirror
//
// split edges are inserted to avoid deleting loops during augment_vertex_set
void VoronoiDiagram::add_split_vertex(HEFace f, Site* s) {
    if (s->isPoint())
        return; // no split-vertices when inserting point-sites
    Site* fs = g[f].site;
    
    // don't search for split-vertex on the start or end face
    if (fs->isPoint() && s->isLine()) {
        if ( fs->position() == s->start() || fs->position() == s->end() )
            return;
    }
    
    
    if ( fs->isPoint() && s->isLine() && s->in_region( fs->position() ) ) {
        
        // 1) find the correct edge
        Point pt1 = fs->position();
        Point pt2 = s->apex_point(pt1);       
        assert( (pt1-pt2).norm() > 0 ); 
        
        EdgeVector split_edges = find_split_edges(f, pt1, pt2);
        // the sought edge should have src on one side of pt1-pt2
        // and trg on the other side of pt1-pt2
        
        BOOST_FOREACH(HEEdge split_edge, split_edges) {
            if ( (g[split_edge].type == SEPARATOR) || (g[split_edge].type == LINESITE) )
                return; // don't place split points on linesites or separators(?)
            HEVertex split_src = g.source(split_edge);
            HEVertex split_trg = g.target(split_edge);
            
            std::cout << " split src=" << g[split_src].index << "("<< g[split_src].dist() << ")";
            std::cout << " trg=" << g[split_trg].index << "("<< g[split_trg].dist() << ") \n";
            std::cout << "is_right src=" << g[split_src].position.is_right(pt1,pt2) << "  trg="<< g[split_trg].position.is_right(pt1,pt2) << "\n";
            // find a point = src + u*(trg-src)
            // with min_t < u < max_t
            // and minimum distance to the pt1-pt2 line
            
            SplitPointError errFunctr(this, split_edge, pt1, pt2); // error functor
            typedef std::pair<double, double> Result;
            boost::uintmax_t max_iter=500;
            boost::math::tools::eps_tolerance<double> tol(30); // bits of tolerance?
            double min_t = std::min( g[split_src].dist() , g[split_trg].dist() );
            double max_t = std::max( g[split_src].dist() , g[split_trg].dist() );
            
            std::cout << " min_t err=" << errFunctr(min_t) << " ";
            std::cout << " max_t err=" << errFunctr(max_t) << "\n";
            Result r1 = boost::math::tools::toms748_solve(errFunctr, min_t, max_t, tol, max_iter);
            
            Point split_pt = g[split_edge].point( r1.first ); 
            HEVertex v = g.add_vertex();
            g[v].type = APEX;
            g[v].status = UNDECIDED;
            g[v].position = split_pt;
            g[v].init_dist( fs->position() );
            std::cout << " new split-vertex " << g[v].index << " t=" << r1.first;
            std::cout << " inserted into edge " << g[split_src].index << "-" << g[split_trg].index  << "\n";
            // 3) insert new SPLIT vertex into the edge
            add_vertex_in_edge(v, split_edge);
        }
    }
}

// generate new voronoi-vertices on all IN-OUT edges 
void VoronoiDiagram::add_vertices( Site* new_site ) {
    assert( !v0.empty() );
    EdgeVector q_edges = find_in_out_edges();       // new vertices generated on these IN-OUT edges
    for( unsigned int m=0; m<q_edges.size(); ++m )  {   
        HEVertex q = g.add_vertex();
        g[q].status = NEW;
        modified_vertices.push_back(q);
        std::cout << "position new vertex " << g[q].index << "\n";
        g[q].position = vpos->position( q_edges[m], new_site ); // set position
        g[q].k3 = vpos->get_k3();
        g[q].init_dist( new_site->apex_point( g[q].position ) ); // set initial clearance-disk
        add_vertex_in_edge(q, q_edges[m] );
    }
}

/// when a new vertex has been added and positioned it is inserted
/// into the edge where it belongs here.
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
    //
    HEEdge twin = g[e].twin;
    HEVertex source = g.source(e); 
    HEVertex target = g.target(e); 
    HEVertex twin_source = g.source(twin); 
    HEVertex twin_target = g.target(twin); 
    assert( source == twin_target );    
    assert( target == twin_source );
    HEFace face = g[e].face;
    HEFace twin_face = g[twin].face;
    HEEdge previous = g.previous_edge(e);
    assert( g[previous].face == g[e].face );
    HEEdge twin_previous = g.previous_edge(twin);
    assert( g[twin_previous].face == g[twin].face );
    
    HEEdge e1 = g.add_edge( source, v ); // e1 and e1 replace e
    HEEdge e2 = g.add_edge( v, target );    
    // preserve the left/right face link
    g[e1].face = face;
    g[e2].face = face;
    // next-pointers
    g[previous].next = e1;
    g[e1].next = e2;
    g[e2].next = g[e].next;
    // k-values
    g[e1].k = g[e].k;
    g[e2].k = g[e].k;
    // type
    g[e1].type = g[e].type;
    g[e2].type = g[e].type;
    
    // edge-parameters
    g[e1].copy_parameters( g[e] );
    g[e2].copy_parameters( g[e] );
    
    HEEdge te1 = g.add_edge( twin_source, v  ); // te1 and te2 replace twin
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
    
    // edge parameters
    g[te1].copy_parameters( g[twin] );
    g[te2].copy_parameters( g[twin] );

    // update the faces (required here?)
    g[face].edge = e1;
    g[twin_face].edge = te1;
    
    // finally, remove the old edge
    g.remove_edge(e);
    g.remove_edge(twin);
}

// add a new face corresponding to the new Site
// call add_new_edge() on all the incident_faces that should be split
HEFace VoronoiDiagram::add_face(Site* s) { 
    HEFace newface =  g.add_face(); 
    g[newface].site = s;
    g[newface].status = NONINCIDENT;
    if (s->isPoint() )
        fgrid->add_face( g[newface] ); 
    s->face = newface;
    return newface;
}

// by adding a NEW-NEW edge, split the face f into one part which is newface, and the other part is the old f
// for linesegment or arc sites we pass in both the k=+1 face newface and the k=-1 face newface2
void VoronoiDiagram::add_edge(HEFace newface, HEFace f, HEFace newface2) {
    if (new_vertex_count(f) == 2) {        // add one new edge src-target to fix the face
        EdgeData ed = find_edge_data(f);
        add_edge(ed, newface, newface2);
    } else if ( new_vertex_count(f) == 4 ) {        // need to add TWO new edges to fix the face
        std::cout << " 4 NEW case!!\n";
        EdgeData ed = find_edge_data(f);
        EdgeData ed2 = find_edge_data(f,ed.v1);
        std::cout << "1st edge is " << g[ed.v1].index << " - " << g[ed.v2].index << "\n";
        std::cout << "2nd edge is " << g[ed2.v1].index << " - " << g[ed2.v2].index << "\n";
        add_edge(ed, newface, newface2);
        add_edge(ed2, newface, newface2);
    } else { 
        assert(0); // unhandled case! panic.
    }
}

// newface = the k=+1 positive offset face
// newface2 = the k=-1 negative offset face
void VoronoiDiagram::add_edge(EdgeData ed, HEFace newface, HEFace newface2) {
    HEEdge new_previous = ed.v1_prv;
    HEVertex new_source = ed.v1;
    HEEdge twin_next = ed.v1_nxt;
    HEEdge twin_previous = ed.v2_prv;
    HEVertex new_target = ed.v2;
    HEEdge new_next = ed.v2_nxt;
    HEFace f = ed.f;
    Site* f_site = g[f].site;
    Site* new_site = g[newface].site;
    
    // both trg and src should be on same side of new site 
    if (g[new_target].k3 != g[new_source].k3) {
        std::cout << " g[" << g[new_target].index << "].k3=" << g[new_target].k3 << " != ";
        std::cout << "g[" << g[new_source].index << "].k3=" << g[new_source].k3<< "\n";
    }
    assert( g[new_target].k3 == g[new_source].k3 );

    //                                           f
    // now connect:   new_previous -> new_source -> new_target -> new_next
    // and:              twin_next <- new_source <- new_target <- twin_previous 
    //                                           new_face   

    // check for potential apex-split
    bool src_sign=true, trg_sign=true;
    if (f_site->isPoint()  && new_site->isLine() ) {
        Point pt1 = f_site->position();
        Point pt2 = new_site->apex_point(pt1);
        src_sign = g[new_source].position.is_right( pt1, pt2 );
        trg_sign = g[new_target].position.is_right( pt1, pt2 );
    } else if (f_site->isPoint() && new_site->isPoint() ) {
        src_sign = g[new_source].position.is_right( f_site->position(), new_site->position() );
        trg_sign = g[new_target].position.is_right( f_site->position(), new_site->position() );
    } else if (f_site->isLine() && new_site->isLine() )  {
        //  a line-line bisector, sign should not matter because there is no sqrt()
        
        std::cout << " LL-edge " << g[new_source].index << " - " << g[new_target].index ;
        std::cout << " f_site->k()= " << f_site->k() << " new_site->k()= "<< new_site->k() << "\n";
        std::cout << " f_site <-> src("<< g[new_source].index << ") = " << g[new_source].position.is_right( f_site->start(), f_site->end() ) << "\n";
        std::cout << " f_site <-> trg("<< g[new_target].index << ") = " << g[new_target].position.is_right( f_site->start(), f_site->end() ) << "\n";
        std::cout << " n_site <-> src("<< g[new_source].index << ") = " << g[new_source].position.is_right( new_site->start(), new_site->end() ) << "\n";
        std::cout << " n_site <-> trg("<< g[new_target].index << ") = " << g[new_target].position.is_right( new_site->start(), new_site->end() ) << "\n";

// this is essentially an in-region test (?)

        assert( !g[new_source].position.is_right( f_site->start(), f_site->end() ) );
        assert( !g[new_target].position.is_right( f_site->start(), f_site->end() ) );
        
// find out if newface or newface2 should be used
        if ( g[new_source].k3 == 1 )
            new_site = g[newface].site; 
        else
            new_site = g[newface2].site; 

        assert( !g[new_source].position.is_right( new_site->start(), new_site->end() ) );
        assert( !g[new_target].position.is_right( new_site->start(), new_site->end() ) );
    } else {
        assert(0);
    }
    
    // both src and trg are on the same side of the new site.
    // so no apex-split is required, just add a single edge.
    if ( src_sign == trg_sign ) {  // add a single src-trg edge
        HEEdge e_new = g.add_edge( new_source, new_target );
        
        g[e_new].next = new_next;
        assert( g[new_next].k == g[new_previous].k );
        g[e_new].k = g[new_next].k; // the next edge is on the same face, so has the correct k-value
        g[e_new].face = f; // src-trg edge has f on its left
        g[new_previous].next = e_new;
        g[f].edge = e_new; 
        g[e_new].set_parameters( f_site, new_site, !src_sign ); 
        // the twin edge that bounds the new face
        HEEdge e_twin = g.add_edge( new_target, new_source );
        g[twin_previous].next = e_twin;
        g[e_twin].next = twin_next;         
        g[e_twin].k = g[new_source].k3; 
        g[e_twin].set_parameters( new_site, f_site, src_sign );

        if (g[e_twin].k == 1) {
            g[e_twin].face = newface; // assumes newface is k==+1 face!
            g[newface].edge = e_twin;
        } else {
            g[e_twin].face = newface2; // assumes newface2 is k==-1 face
            g[newface2].edge = e_twin;
        }
        g.twin_edges(e_new,e_twin);
        std::cout << " added edge " << g[new_target].index << "(" << g[new_target].dist() <<")";
        std::cout << " - " << g[new_source].index << "(" << g[new_source].dist() << ")";
        std::cout << " f=" << g[e_new].face << " k=" << g[e_twin].k;
        std::cout << " twf=" << g[e_twin].face << " twk=" << g[e_new].k;
        std::cout << " type=" << g[e_twin].type <<   " \n";
        //std::cout << " k3 target-source: "<<  g[new_target].k3 << " - " << g[new_source].k3 << "\n";
    } else {
        // need to do apex-split
        //                         f               f  
        //   new_prv -> new_src -- e1 ----> APEX --e2 ---> new_trg -> new_nxt
        //   twi_nxt <- new_src <- e1_tw -- APEX <-e2_tw-- new_trg <- twn_prv    
        //                       new1/new2         new1/new2
        //   
        HEVertex apex = g.add_vertex();
        g[apex].type = APEX;
        g[apex].status = NEW;
        HEEdge e1 = g.add_edge( new_source, apex);
        HEEdge e2 = g.add_edge( apex, new_target);
        g[e1].set_parameters(f_site,new_site,!src_sign);
        g[e2].set_parameters(f_site,new_site,!trg_sign);
        g[new_previous].next = e1;
        g[e1].next=e2;
        g[e2].next=new_next;
        g[e1].face = f;
        g[e2].face = f;
        g[f].edge = e1;
        
        assert( g[new_next].k == g[new_previous].k );

        g[e1].k = g[new_next].k;
        g[e2].k = g[new_next].k;
    // twin edges
        HEEdge e1_tw = g.add_edge( apex, new_source );
        HEEdge e2_tw = g.add_edge( new_target, apex );
        g[e1_tw].set_parameters(new_site, f_site, src_sign);
        g[e2_tw].set_parameters(new_site, f_site, trg_sign);
        g[twin_previous].next = e2_tw;
        g[e2_tw].next = e1_tw;
        g[e1_tw].next = twin_next;
        g[e1_tw].k = g[new_source].k3;
        g[e2_tw].k = g[new_source].k3;
        
        if (g[e1_tw].k == 1) {
            g[e1_tw].face = newface;
            g[e2_tw].face = newface;
            g[newface].edge = e1_tw;
        } else {
            g[e1_tw].face = newface2;
            g[e2_tw].face = newface2;
            g[newface2].edge = e1_tw;
        }        
        g.twin_edges(e1,e1_tw);
        g.twin_edges(e2,e2_tw);

    // position the apex
        double min_t = g[e1].minimum_t(f_site,new_site);
        g[apex].position = g[e1].point(min_t);
        std::cout << "   apex= " << g[apex].index ; // << "\n";
        std::cout << " t: src=" << g[new_source].dist() << " tmin= " << min_t << " trg= " << g[new_target].dist() << "\n";
        
        std::cout << " added edge " << g[new_target].index << "-" << g[apex].index << " f=" << g[e1_tw].face << " k= " << g[e1_tw].k << " twin_f= " <<  g[e1].face << " \n";
        std::cout << " added edge " << g[apex].index << "-" << g[new_source].index << " f=" << g[e2_tw].face << " k= " << g[e2_tw].k << "\n";
        
        //std::cout << " k3 target-source: "<<  g[new_target].k3 << " - " << g[new_source].k3 << "\n";

        
        g[apex].init_dist(f_site->apex_point(g[apex].position));
        modified_vertices.push_back( apex );
        
    }
}

// count number of NEW vertices on the given face
int VoronoiDiagram::new_vertex_count(HEFace f) {
    VertexVector face_verts = g.face_vertices(f);
    int count=0;
    BOOST_FOREACH(HEVertex v, face_verts) {
        if (g[v].status == NEW)
            count++;
    }    
    return count;
}

// on a face which has IN and OUT-vertices, find the sequence
// OUT-OUT-OUT-..-OUT-NEW(v1)-IN-...-IN-NEW(v2)-OUT-OUT
// and return v1/v2 together with their previous and next edges
EdgeData VoronoiDiagram::find_edge_data(HEFace f, HEVertex v)  {
    EdgeData ed;
    ed.f = f;
    HEEdge current_edge = g[f].edge; // start on some edge of the face
    bool found = false;
    int count=0;                             
    while (!found) {
        HEVertex current_vertex = g.target( current_edge );
        HEEdge next_edge = g[current_edge].next;
        HEVertex next_vertex = g.target( next_edge );
        if ( g[current_vertex].status == OUT ) {
            if ( g[next_vertex].status == NEW && next_vertex!=v) {
                    ed.v1 = next_vertex;
                    ed.v1_prv = next_edge;
                    ed.v1_nxt = g[next_edge].next;
                    found = true;                 
            }  
        }
        current_edge = g[current_edge].next;   
        count++;
        assert(count<100); // some reasonable max number of edges in face, to avoid infinite loop
    }
    // now search for v2
    count=0; found=false;
    while (!found) {
        HEVertex current_vertex = g.target( current_edge );
        HEEdge next_edge = g[current_edge].next;
        HEVertex next_vertex = g.target( next_edge );
        if ( g[current_vertex].status == IN ) {
            if ( g[next_vertex].status == NEW) { // IN-IN-NEW(v2)
                    ed.v2 = next_vertex;
                    ed.v2_prv = next_edge;
                    ed.v2_nxt = g[next_edge].next;
                    found = true;                 
            }
        }
        current_edge = g[current_edge].next;   
        count++;
        assert(count<100); // some reasonable max number of edges in face, to avoid infinite loop
    }
    return ed;
}

// start on g[newface].edge, walk around the face and repair the next-pointers
// this is called on the newly created face after all NEW-NEW edges have been added
void VoronoiDiagram::repair_face( HEFace f ) {
    //std::cout << "repair_face( " << f << " )\n";
    HEEdge current_edge = g[f].edge; 
    HEEdge start_edge = current_edge;
    do {
        HEVertex current_target = g.target( current_edge ); // an edge on the new face
        HEVertex current_source = g.source( current_edge );
        bool found_next_edge= false;
        BOOST_FOREACH( HEEdge edge, g.out_edges( current_target ) ) { // loop through potential "next" candidates
            HEVertex out_target = g.target( edge );
            //std::cout << g[current_source].index << " - " << g[current_target].index << " f= "<< g[current_edge].face;
            //std::cout << " target= " << g[out_target].index << " f=" << g[edge].face << " "; 
            if ( ( (g[out_target].status == NEW) || (g[out_target].type == ENDPOINT) ) && (g[edge].face == f) ) { // the next vertex along the face should be "NEW"
                if ( out_target != current_source ) { // but not where we came from
                    g[current_edge].next = edge; // this is the edge we want to take
                    //std::cout << " VALID!.\n";
                    found_next_edge = true;
                    assert( vd_checker->current_face_equals_next_face( current_edge ) );
                }
            } 
        }
        assert(found_next_edge); // must find a next-edge!
        current_edge = g[current_edge].next; // jump to the next edge
    } while (g[current_edge].next != start_edge);
    //std::cout << " repair done:\n";
    //print_face(f);
}

void VoronoiDiagram::remove_vertex_set() {
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
    return (in_count >= 2); // two or more adjacent IN-vertices might create a loop
}

// do any of the three faces that are adjacent to the given IN-vertex v have an IN-vertex ?
// predicate C5 i.e. "connectedness"  from Sugihara&Iri 1992 "one million" paper
bool VoronoiDiagram::predicate_c5(HEVertex v) {
    if (g[v].type == APEX ) { return true; } // ?
    
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
            else if ( w!=v && ( g[w].type == ENDPOINT || g[w].type == APEX) ) // if we are next to an ENDPOINT, then ok(?)
                face_ok=true;
        }
        if (!face_ok)
            all_found=false;
    }
    return all_found; // if this returns false, we mark a vertex OUT, on topology grounds.
}

void VoronoiDiagram::print_faces() {
    for( HEFace f=0;f<g.num_faces();f++) {
        print_face(f);
    }
}
void VoronoiDiagram::print_face(HEFace f) {
    std::cout << " Face " << f << ": vert_id(status) ";
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

void VoronoiDiagram::print_edges(EdgeVector& q) {
    BOOST_FOREACH( HEEdge e, q ) {
        HEVertex src = g.source(e);
        HEVertex trg = g.target(e);
        std::cout << g[src].index << "-" << g[trg].index << "\n";
    }
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
