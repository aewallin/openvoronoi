/* 
 *  Copyright 2010-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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
#include <boost/tuple/tuple.hpp>
#include <boost/assign/list_of.hpp>

#include "voronoidiagram.hpp"

#include "checker.hpp"
#include "common/numeric.hpp" // for diangle

namespace ovd {

/// \brief create a VoronoiDiagram
/// \param far is the radius of a circle within which all sites must be located. use far==1.0
/// \param n_bins is the number of bins for FaceGrid, the bucket-search for nearest-neighbors 
///        used in insert_point_site(). Use roughly sqrt(N) for a voronoi-diagram with N sites.
VoronoiDiagram::VoronoiDiagram(double far, unsigned int n_bins) {
    //fgrid = new FaceGrid(far, n_bins); // helper-class for nearest-neighbor search 
    kd_tree = new kdtree::KDTree<kd_point>(2);
    
    vd_checker = new VoronoiDiagramChecker( g ); // helper-class that checks topology/geometry
     
    vpos = new VertexPositioner( g ); // helper-class that positions vertices
    far_radius=far;
    initialize();
    num_psites=3;
    num_lsites=0;
    num_asites=0;
    reset_vertex_count();
    debug = false;
}

/// \brief delete allocated resources: FaceGrid, checker, positioner
VoronoiDiagram::~VoronoiDiagram() { 
    //std::cout << "~VoronoiDiagram()\n";
    //delete fgrid; 
    delete kd_tree;
    delete vpos;
    delete vd_checker;
    
    //std::cout << "~VoronoiDiagram() DONE.\n";
}

/// \brief initialize the diagram with three generators
///
/// add one vertex at origo and three vertices at 'infinity' and their associated edges
void VoronoiDiagram::initialize() {
    using namespace boost::assign;
    double far_multiplier = 6;
    // initial generators/sites:
    Point gen1 = Point( 0, 3.0*far_radius );
    Point gen2 = Point( -3.0*sqrt(3.0)*far_radius/2.0, -3.0*far_radius/2.0 );
    Point gen3 = Point( +3.0*sqrt(3.0)*far_radius/2.0, -3.0*far_radius/2.0 );

    // initial vd-vertices
    Point vd1 = Point(            0                 , -3.0*far_radius*far_multiplier    );
    Point vd2 = Point( +3.0*sqrt(3.0)*far_radius*far_multiplier/2.0, +3.0*far_radius*far_multiplier/2.0);
    Point vd3 = Point( -3.0*sqrt(3.0)*far_radius*far_multiplier/2.0, +3.0*far_radius*far_multiplier/2.0);
    // add init vertices
    HEVertex v00 = g.add_vertex( VoronoiVertex( Point(0,0), UNDECIDED, NORMAL, gen1 ) );
    HEVertex v01 = g.add_vertex( VoronoiVertex( vd1, OUT, OUTER, gen3) );
    HEVertex v02 = g.add_vertex( VoronoiVertex( vd2, OUT, OUTER, gen1) );
    HEVertex v03 = g.add_vertex( VoronoiVertex( vd3, OUT, OUTER, gen2) );
    // add initial sites to graph 
    HEVertex vert1 = g.add_vertex( VoronoiVertex( gen1 , OUT, POINTSITE) );
    HEVertex vert2 = g.add_vertex( VoronoiVertex( gen2 , OUT, POINTSITE) );
    HEVertex vert3 = g.add_vertex( VoronoiVertex( gen3 , OUT, POINTSITE) );

    // apex-points on the three edges: 
    HEVertex a1 = g.add_vertex( VoronoiVertex( 0.5*(gen2+gen3), UNDECIDED, APEX, gen2 ) );
    HEVertex a2 = g.add_vertex( VoronoiVertex( 0.5*(gen1+gen3), UNDECIDED, APEX, gen3 ) );
    HEVertex a3 = g.add_vertex( VoronoiVertex( 0.5*(gen1+gen2), UNDECIDED, APEX, gen1 ) );

    // add face 1: v0-v1-v2 which encloses gen3
    HEEdge e1_1 =  g.add_edge( v00 , a1 );    
    HEEdge e1_2 =  g.add_edge( a1 , v01 );   
    HEEdge e2   =  g.add_edge( v01, v02 );
    HEEdge e3_1 =  g.add_edge( v02, a2  ); 
    HEEdge e3_2 =  g.add_edge( a2 , v00 ); 
    HEFace f1   =  g.add_face(); 
    g[f1].site  = new PointSite(gen3,f1, vert3);
    g[f1].status = NONINCIDENT;
    //fgrid->add_face( f1, gen3 ); // for grid search
    kd_tree->insert( kd_point(gen3,f1) );
    g.set_next_cycle( list_of(e1_1)(e1_2)(e2)(e3_1)(e3_2) , f1 ,1);

    // add face 2: v0-v02-v03 which encloses gen1
    HEEdge e4_1 = g.add_edge( v00, a2  );
    HEEdge e4_2 = g.add_edge( a2, v02 );
    HEEdge e5   = g.add_edge( v02, v03  );
    HEEdge e6_1 = g.add_edge( v03, a3 );
    HEEdge e6_2 = g.add_edge( a3, v00 ); 
    HEFace f2   =  g.add_face();
    g[f2].site  = new PointSite(gen1,f2, vert1);
    g[f2].status = NONINCIDENT;    
    //fgrid->add_face( f2, gen1 );
    kd_tree->insert( kd_point(gen1,f2) );
    g.set_next_cycle( list_of(e4_1)(e4_2)(e5)(e6_1)(e6_2) , f2 ,1);

    // add face 3: v0-v3-v1 which encloses gen2
    HEEdge e7_1 = g.add_edge( v00, a3 );  
    HEEdge e7_2 = g.add_edge( a3 , v03 );   
    HEEdge e8   = g.add_edge( v03, v01 );
    HEEdge e9_1 = g.add_edge( v01, a1  ); 
    HEEdge e9_2 = g.add_edge( a1 , v00 ); 
    HEFace f3   =  g.add_face();
    g[f3].site  = new PointSite(gen2,f3, vert2); // this constructor needs f3...
    g[f3].status = NONINCIDENT;    
    //fgrid->add_face( f3, gen2 );
    kd_tree->insert( kd_point(gen2,f3) );
    g.set_next_cycle( list_of(e7_1)(e7_2)(e8)(e9_1)(e9_2) , f3 , 1);    

    // set type. 
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
    g.twin_edges(e1_1,e9_2);
    g.twin_edges(e1_2,e9_1);
    g[e2].twin = HEEdge(); // the outermost edges have invalid twins
    g[e5].twin = HEEdge();
    g[e8].twin = HEEdge();
    g.twin_edges(e3_1, e4_2);
    g.twin_edges(e3_2, e4_1);
    g.twin_edges(e6_1, e7_2);
    g.twin_edges(e6_2, e7_1);
    
    assert( vd_checker->is_valid() );
}


/// \brief insert a PointSite into the diagram 
///
/// \param p position of site
/// \param step (optional, for debugging) stop at this step
/// \return integer handle to the inserted point. use this integer when inserting lines/arcs with insert_line_site
///
/// \details
/// \attention All PointSite:s must be inserted before any LineSite:s or ArcSite:s are inserted. 
/// \attention It is an error to insert duplicate PointSite:s (i.e. points with the same x,y coordinates)
///
/// All PointSite:s must be inserted before any LineSite:s or ArcSite:s are inserted.
/// This is roughly "algorithm A" from the Sugihara-Iri 1994 paper, page 15/50
///
/// -# find the face that is closest to the new site, see FaceGrid
/// -# among the vertices on the closest face, find the seed vertex, see find_seed_vertex()
/// -# grow the tree of IN-vertices, see augment_vertex_set()
/// -# add new voronoi-vertices on all IN-OUT edges so they becone IN-NEW-OUT, see add_vertices()
/// -# add new face by splitting each INCIDENT face into two parts by inserting a NEW-NEW edge. see add_edges()
/// -# repair the next-pointers of faces that have been modified. see repair_face()
/// -# remove IN-IN edges and IN-NEW edges, see remove_vertex_set()
/// -# reset vertex/face status to be ready for next incremental operation, see reset_status()
int VoronoiDiagram::insert_point_site(const Point& p, int step) {
    num_psites++;
    //int current_step=1;
    if (p.norm() >= far_radius ) {
        std::cout << "openvoronoi error. All points must lie within unit-circle. You are trying to add p= " << p 
        << " with p.norm()= " << p.norm() << "\n";
    } 
    assert( p.norm() < far_radius );     // only add vertices within the far_radius circle
    
    HEVertex new_vert = g.add_vertex( VoronoiVertex(p,OUT,POINTSITE) );
    PointSite* new_site =  new PointSite(p);
    new_site->v = new_vert;
    vertex_map.insert( VertexMapPair(g[new_vert].index,new_vert) ); // so that we can find the descriptor later based on its index
    std::pair<kd_point,bool> nearest = kd_tree->nearest( kd_point(p) );
    assert( nearest.second );
    HEVertex v_seed = find_seed_vertex( nearest.first.face , new_site);
    mark_vertex( v_seed, new_site );
//if (step==current_step) return -1; current_step++;
    augment_vertex_set( new_site ); // grow the tree to maximum size
//if (step==current_step) return -1; current_step++;
    add_vertices( new_site );  // insert new vertices on IN-OUT edges
//if (step==current_step) return -1; current_step++;
    HEFace newface = add_face( new_site );
    g[new_vert].face = newface; // Vertices that correspond to point-sites have their .face property set!
    BOOST_FOREACH( HEFace f, incident_faces ) { // add NEW-NEW edges on all INCIDENT faces
        add_edges(newface, f);
    }
//if (step==current_step) return -1; current_step++;
    repair_face( newface  );
    if (debug) { std::cout << " new face: "; g.print_face( newface ); }
    remove_vertex_set(); // remove all IN vertices and adjacent edges
//if (step==current_step) return -1; current_step++;
    reset_status(); // reset all vertices to UNDECIDED
 
    assert( vd_checker->face_ok( newface ) );
    assert( vd_checker->is_valid() );
 
    return g[new_vert].index;
}

/// \brief insert a LineSite into the diagram
///
/// \param idx1 int handle to startpoint of line-segment
/// \param idx2 int handle to endpoint of line-segment
/// \param step (optional, for debug) stop at step
///
/// \details
/// \attention All PointSite:s must be inserted before any LineSite:s are inserted. 
///   All LineSite:s should be inserted before any ArcSite:s are inserted.
/// \attention It is an error to insert a LineSite that intersects an existing LineSite in the diagram!
///
/// The basic idea of the incremental diagram update is similar to that in insert_point_site().
/// The major differences are:
/// - the handling of null-faces at the endpoints of the LineSite.
/// - addition of ::SEPARATOR edges
/// - addition of ::SPLIT vertices during augment_vertex_set()
/// - removal of ::SPLIT vertices at the end
///
/// The steps of the algorithm are:
/// -# based on \a idx1 and \a idx2, look up the corresponding vertex descriptors. It is an error if these are not found.
/// -# find a seed-vertex
/// -# grow the delete-tree of ::IN vertices.
/// -# create or modify the null-faces at the startpoint and endpoint of the LineSite
/// -# create and add ::LINESITE pseudo-edges
/// -# add ::NEW vertices on all ::IN-::OUT edges.
/// -# add up to four ::SEPARATOR edges, where applicable
/// -# add non-separator edges by calling add_edges() on all ::INCIDENT faces
/// -# repair the next-pointers of faces that have been modified. see repair_face()
/// -# remove IN-IN edges and IN-NEW edges, see remove_vertex_set()
/// -# remove ::SPLIT vertices
/// -# reset vertex/face status to be ready for next incremental operation, see reset_status()
bool VoronoiDiagram::insert_line_site(int idx1, int idx2, int step) {
    num_lsites++;
    int current_step=1;
    // find the vertices corresponding to idx1 and idx2
    HEVertex start=HEVertex(), end=HEVertex();
    boost::tie(start,end) = find_endpoints(idx1,idx2);
    g[start].status=OUT;
    g[end].status=OUT;   
    g[start].zero_dist();
    g[end].zero_dist();
    
    if (debug) std::cout << "insert_line_site( " << g[start].index << " - " << g[end].index << " )\n";
    
    // create a point which is left of src->trg
    // determine k (offset-dir) for this point
    // then we know which site/face is the k==+1 and which is k==-1
    Point src_se = g[start].position;
    Point trg_se = g[end  ].position;
    Point left = 0.5*(src_se+trg_se) + (trg_se-src_se).xy_perp(); // this is used below and in find_null_face()
    #ifndef NDEBUG
    bool linesite_k_sign = left.is_right(src_se,trg_se); 
    #endif
    
if (step==current_step) return false; current_step++;

    LineSite* pos_site;
    LineSite* neg_site;
    assert(!linesite_k_sign);
    // 2012-03-20: coverage testing shows this never happens:
    //if ( linesite_k_sign ) {
    //    pos_site = new LineSite( g[start].position, g[end  ].position , +1);
    //    neg_site = new LineSite( g[end  ].position, g[start].position , -1);
    //} else {
    pos_site = new LineSite( g[end  ].position, g[start].position , +1);
    neg_site = new LineSite( g[start].position, g[end  ].position , -1);
    //}

if (step==current_step) return false; current_step++;

    HEFace seed_face = g[start].face; // assumes this point-site has a face!
    
    // on the face of start-point, find the seed vertex
    HEVertex v_seed = find_seed_vertex(seed_face, pos_site ) ;
    if (debug) std::cout << " start face seed  = " << g[v_seed].index << "\n";
    mark_vertex( v_seed, pos_site  );

if (step==current_step) return false; current_step++;

    augment_vertex_set( pos_site  ); // it should not matter if we use pos_site or neg_site here
    // todo(?) sanity checks:
    // check that end_face is INCIDENT? 
    // check that tree (i.e. v0) includes end_face_seed ?
    if (debug) { std::cout << " delete-set |v0|="<< v0.size() <<" : ";  g.print_vertices(v0); }

if (step==current_step) return false; current_step++;

    // process the null-faces here
    HEVertex seg_start, seg_end; // new segment end-point vertices. these are created here.
    HEFace start_null_face, end_null_face; // either existing or new null-faces at endpoints
    HEVertex pos_sep_start = HEVertex(); // optional positive separator vertex at start
    HEVertex neg_sep_start = HEVertex(); // optional negative separator vertex at start 
    HEVertex pos_sep_end = HEVertex();   // optional positive separator vertex at end
    HEVertex neg_sep_end = HEVertex();   // optional negative separator vertex at end
    HEFace start_to_null = g.HFace(); // when a point-site completely disappears, we "null" the face belonging to this pointsite
    HEFace   end_to_null = g.HFace();
    // returns new seg_start/end vertices, new or existing null-faces, and separator endpoints (if separators should be added)
    Point dir2 = g[start].position - g[end].position;
    Point dir1 = g[end].position - g[start].position;
    boost::tie(seg_start, start_null_face, pos_sep_start, neg_sep_start, start_to_null) = find_null_face(start, end  , left, dir1, pos_site);
    boost::tie(seg_end  , end_null_face  , pos_sep_end  , neg_sep_end  , end_to_null  ) = find_null_face(end  , start, left, dir2, pos_site);

    // now safe to set the zero-face edge
    // in the collinear case, set the edge for the face that "disappears" to a null edge
    if (start_to_null!=g.HFace())
        g[start_to_null].edge = g[start_null_face].edge; 
    if (end_to_null!=g.HFace())
        g[end_to_null].edge = g[end_null_face].edge; 

if (step==current_step) return false; current_step++;
    // create LINESITE pseudo edges and faces
    HEFace pos_face, neg_face; 
    HEEdge pos_edge, neg_edge;
    {
        // 2012-03-20: coverage-testing shows this never happens!
        //if ( linesite_k_sign ) {
        //    boost::tie(pos_edge, neg_edge) = g.add_twin_edges( seg_start,   seg_end );
        //    g[pos_edge].inserted_direction = true;
        //    g[neg_edge].inserted_direction = false;
        //} else {
        assert(!linesite_k_sign);
            boost::tie( pos_edge, neg_edge) = g.add_twin_edges( seg_end  ,seg_start );
            g[pos_edge].inserted_direction = false;
            g[neg_edge].inserted_direction = true;
        //}
        g[pos_edge].type = LINESITE;
        g[neg_edge].type = LINESITE;
        g[pos_edge].k = +1;
        g[neg_edge].k = -1;
        pos_face = add_face( pos_site ); //  this face to the left of start->end edge    
        neg_face = add_face( neg_site ); //  this face is to the left of end->start edge
        g[pos_face].edge = pos_edge;
        g[neg_face].edge = neg_edge;
        g[pos_edge].face = pos_face;
        g[neg_edge].face = neg_face;
        
        // associate sites with LINESITE edges
        pos_site->e = pos_edge;
        neg_site->e = neg_edge;
        
        if (debug) std::cout << "created faces: pos_face=" << pos_face << " neg_face=" << neg_face << "\n";   
    }

if (step==current_step) return false; current_step++;

    add_vertices( pos_site );  // add NEW vertices on all IN-OUT edges.

if (step==current_step) return false; current_step++;

    { // add SEPARATORS
        // find SEPARATOR targets first
        typedef boost::tuple<HEEdge, HEVertex, HEEdge,bool> SepTarget;
        SepTarget pos_start_target, neg_start_target;
        pos_start_target = find_separator_target( g[start].face , pos_sep_start);
        neg_start_target = find_separator_target( g[start].face , neg_sep_start);
        
        // add positive separator edge at start
        add_separator( g[start].face , start_null_face, pos_start_target, pos_sep_start, g[pos_face].site , g[neg_face].site );
        
if (step==current_step) return false; current_step++;
        
        // add negative separator edge at start
        add_separator( g[start].face , start_null_face, neg_start_target, neg_sep_start, g[pos_face].site , g[neg_face].site );
        g[ g[start].face ].status = NONINCIDENT; // face is now done.
        assert( vd_checker->face_ok( g[start].face ) );
        
if (step==current_step) return false; current_step++;

        SepTarget pos_end_target, neg_end_target;
        pos_end_target = find_separator_target( g[end].face ,  pos_sep_end);
        neg_end_target = find_separator_target( g[end].face , neg_sep_end);
        // add positive separator edge at end
        add_separator( g[end].face , end_null_face, pos_end_target, pos_sep_end, g[pos_face].site , g[neg_face].site );
        
if (step==current_step) return false; current_step++;
        
        // add negative separator edge at end
        add_separator( g[end].face , end_null_face, neg_end_target, neg_sep_end, g[pos_face].site , g[neg_face].site );
        g[ g[end].face ].status = NONINCIDENT;
        assert( vd_checker->face_ok( g[end].face ) );

        if(debug) std::cout << "all separators  done.\n";
    } // end SEPARATORS

if (step==current_step) return false; current_step++;

// add non-separator edges by calling add_edges on all INCIDENT faces
    {
        if(debug) std::cout << "adding edges.\n";
        BOOST_FOREACH( HEFace f, incident_faces ) {
            if ( g[f].status == INCIDENT )  {// end-point faces already dealt with in add_separator()
                if(debug) { 
                    std::cout << " add_edges f= " << f << "\n";
                    g.print_face(f);
                }
                add_edges( pos_face, f, neg_face, std::make_pair(seg_start,seg_end)); // each INCIDENT face is split into two parts: newface and f
            }
        }
        if(debug) std::cout << "adding edges. DONE.\n";
    }

if (step==current_step) return false; current_step++;

// new vertices and edges inserted. remove the delete-set, repair faces.

    remove_vertex_set();

    if (debug) { std::cout << "will now repair pos/neg faces: " << pos_face << " " << neg_face << "\n"; }

    repair_face( pos_face, std::make_pair(seg_start,seg_end), 
                           std::make_pair(start_to_null,end_to_null),
                           std::make_pair(start_null_face,end_null_face) );
    assert( vd_checker->face_ok( pos_face ) );

    repair_face( neg_face, std::make_pair(seg_start,seg_end), 
                           std::make_pair(start_to_null,end_to_null),
                           std::make_pair(start_null_face,end_null_face) );
    assert( vd_checker->face_ok( neg_face ) );
    if (debug) { std::cout << "faces repaired.\n"; }

if (step==current_step) return false; current_step++;

    // we are done and can remove split-vertices
    BOOST_FOREACH(HEFace f, incident_faces) {
        remove_split_vertex(f);
    }
    reset_status();


    assert( vd_checker->face_ok( start_null_face ) );
    assert( vd_checker->face_ok( end_null_face ) );
    assert( vd_checker->face_ok( pos_face ) );
    assert( vd_checker->face_ok( neg_face ) );    
    assert( vd_checker->is_valid() );

    return true; 
}

/// \brief insert a circular arc Site into the diagram
/// \param idx1 index of start vertex
/// \param idx2 index of end vertex
/// \param center center Point of arc
/// \param cw bool flag true=CW arc, false=CCW arc
/// \param step for debug, stop algorithm at this sub-step
void VoronoiDiagram::insert_arc_site(int idx1, int idx2, const Point& center, bool cw, int step) {
    num_asites++;
    int current_step=1;
    // find the vertices corresponding to idx1 and idx2
    HEVertex start=HEVertex(), end=HEVertex();
    boost::tie(start,end) = find_endpoints(idx1,idx2);
    g[start].status=OUT;
    g[end].status=OUT;   
    g[start].zero_dist();
    g[end].zero_dist();
    double radius = (g[start].position - center).norm();
    if (debug) {
        std::cout << "insert_arc_site( " << g[start].index << " - " << g[end].index << "," ;
        std::cout << " c= " << center << ", cw= " << cw;
        std::cout << " )\n";
        std::cout << " radius= " << radius << "\n";
    }
if (step==current_step) return; current_step++;
    
    ArcSite* pos_site;
    ArcSite* neg_site;
    if (cw) {
        pos_site = new ArcSite( g[end  ].position, g[start].position , center, cw);
        neg_site = new ArcSite( g[start].position, g[end  ].position , center, !cw);
    } else {
        pos_site = new ArcSite( g[start].position, g[end].position , center, !cw);
        neg_site = new ArcSite( g[end].position, g[start].position , center, cw);
    }
    
    if (debug) {
        std::cout << " pos site =  " << pos_site->str2() << "\n";
        std::cout << " neg site =  " << neg_site->str2() << "\n";
    }
if (step==current_step) return; current_step++;

    HEFace seed_face = g[start].face; // assumes this point-site has a face!
    // on the face of start-point, find the seed vertex
    HEVertex v_seed = find_seed_vertex(seed_face, pos_site ) ;
    if (debug) std::cout << " start face seed  = " << g[v_seed].index << "\n";
    mark_vertex( v_seed, pos_site  );

if (step==current_step) return; current_step++;
    
    augment_vertex_set( pos_site  ); // it should not matter if we use pos_site or neg_site here
    // todo(?) sanity checks:
    // check that end_face is INCIDENT? 
    // check that tree (i.e. v0) includes end_face_seed ?
    if (debug) { std::cout << " delete-set |v0|="<< v0.size() <<" : ";  g.print_vertices(v0); }

if (step==current_step) return; current_step++;
    // process the null-faces here
    HEVertex seg_start, seg_end; // new segment end-point vertices. these are created here.
    HEFace start_null_face, end_null_face; // either existing or new null-faces at endpoints
    HEVertex pos_sep_start = HEVertex(); // optional positive separator vertex at start
    HEVertex neg_sep_start = HEVertex(); // optional negative separator vertex at start 
    HEVertex pos_sep_end = HEVertex();   // optional positive separator vertex at end
    HEVertex neg_sep_end = HEVertex();   // optional negative separator vertex at end
    HEFace start_to_null = g.HFace(); // when a point-site completely disappears, we "null" the face belonging to this pointsite
    HEFace   end_to_null = g.HFace();
    Point src_se = g[start].position;
    Point trg_se = g[end  ].position;

    Point left = 0.5*(src_se+trg_se) + (trg_se-src_se).xy_perp(); // this is used below and in find_null_face()
    // returns new seg_start/end vertices, new or existing null-faces, and separator endpoints (if separators should be added)
    Point dir1,dir2;
    if (cw) {
        dir1 = (g[start].position - center).xy_perp();
        dir2 = -1*(g[end].position - center).xy_perp();
    } else {
        dir1 = -1*(g[start].position - center).xy_perp();
        dir2 = (g[end].position - center).xy_perp();
    }
    if (debug) std::cout << "find_null_face( " << g[start].index << " )\n";
    
    boost::tie(seg_start, start_null_face, pos_sep_start, neg_sep_start, start_to_null) = find_null_face(start, end  , left, dir1, pos_site);
    if (debug) std::cout << "find_null_face( " << g[end].index << " )\n";
    boost::tie(seg_end  , end_null_face  , pos_sep_end  , neg_sep_end  , end_to_null  ) = find_null_face(end  , start, left, dir2, pos_site);

    // now safe to set the zero-face edge
    // in the collinear case, set the edge for the face that "disappears" to a null edge
    if (start_to_null!=g.HFace())
        g[start_to_null].edge = g[start_null_face].edge; 
    if (end_to_null!=g.HFace())
        g[end_to_null].edge = g[end_null_face].edge; 

if (step==current_step) return; current_step++;
    // create pseudo edges (sites) and faces
    HEFace pos_face, neg_face; 
    HEEdge pos_edge, neg_edge;
    {
        if (cw)
            boost::tie( pos_edge, neg_edge) = g.add_twin_edges( seg_end  ,seg_start );
        else
            boost::tie( neg_edge, pos_edge) = g.add_twin_edges( seg_end  ,seg_start );

        g[pos_edge].inserted_direction = false;
        g[neg_edge].inserted_direction = true;
        g[pos_edge].type = ARCSITE;
        g[neg_edge].type = ARCSITE;
        g[pos_edge].k = +1;
        g[neg_edge].k = -1;
        pos_face = add_face( pos_site ); //  this face to the left of start->end edge    
        neg_face = add_face( neg_site ); //  this face is to the left of end->start edge
        g[pos_face].edge = pos_edge;
        g[neg_face].edge = neg_edge;
        g[pos_edge].face = pos_face;
        g[neg_edge].face = neg_face;
        
        //if (debug) std::cout << "created ARCSITE edge " << g[seg_end].index << " -[f" <<  << "]- << g[seg_start].index << "\n";  
        
        // associate sites with ARCSITE edges
        pos_site->e = pos_edge;
        neg_site->e = neg_edge;
        
        if (debug) std::cout << "created faces: pos_face=" << pos_face << " neg_face=" << neg_face << "\n";   
    }

if (step==current_step) return; current_step++;

    add_vertices( pos_site );  // add NEW vertices on all IN-OUT edges.
    
if (step==current_step) return; current_step++;

    { // add SEPARATORS
    
        // find SEPARATOR targets first
        typedef boost::tuple<HEEdge, HEVertex, HEEdge,bool> SepTarget;
        SepTarget pos_start_target, neg_start_target;
        pos_start_target = find_separator_target( g[start].face , pos_sep_start);
        neg_start_target = find_separator_target( g[start].face , neg_sep_start);
        
        // add positive separator edge at start
        add_separator( g[start].face , start_null_face, pos_start_target, pos_sep_start, g[pos_face].site , g[neg_face].site );
        
if (step==current_step) return; current_step++;
        
        // add negative separator edge at start
        add_separator( g[start].face , start_null_face, neg_start_target, neg_sep_start, g[pos_face].site , g[neg_face].site );
        g[ g[start].face ].status = NONINCIDENT; // face is now done.
        assert( vd_checker->face_ok( g[start].face ) );
        
if (step==current_step) return; current_step++;

        SepTarget pos_end_target, neg_end_target;
        pos_end_target = find_separator_target( g[end].face ,  pos_sep_end);
        neg_end_target = find_separator_target( g[end].face , neg_sep_end);
        // add positive separator edge at end
        add_separator( g[end].face , end_null_face, pos_end_target, pos_sep_end, g[pos_face].site , g[neg_face].site );
        
if (step==current_step) return; current_step++;
        
        // add negative separator edge at end
        add_separator( g[end].face , end_null_face, neg_end_target, neg_sep_end, g[pos_face].site , g[neg_face].site );
        g[ g[end].face ].status = NONINCIDENT;
        assert( vd_checker->face_ok( g[end].face ) );

        if(debug) std::cout << "all separators  done.\n";
    } // end SEPARATORS

if (step==current_step) return; current_step++;
// add non-separator edges by calling add_edges on all INCIDENT faces
    {
        if(debug) std::cout << "adding edges.\n";
        BOOST_FOREACH( HEFace f, incident_faces ) {
            if ( g[f].status == INCIDENT )  {// end-point faces already dealt with in add_separator()
                if(debug) { 
                    std::cout << " add_edges f= " << f << "\n";
                    g.print_face(f);
                }
                add_edges( pos_face, f, neg_face, std::make_pair(seg_start,seg_end)); // each INCIDENT face is split into two parts: newface and f
            }
        }
        if(debug) std::cout << "adding edges. DONE.\n";
    }

if (step==current_step) return; current_step++;

// new vertices and edges inserted. remove the delete-set, repair faces.

    remove_vertex_set();

    if (debug) { std::cout << "will now repair pos/neg faces: " << pos_face << " " << neg_face << "\n"; }

    repair_face( pos_face, std::make_pair(seg_start,seg_end), 
                           std::make_pair(start_to_null,end_to_null),
                           std::make_pair(start_null_face,end_null_face) );
    assert( vd_checker->face_ok( pos_face ) );

    repair_face( neg_face, std::make_pair(seg_start,seg_end), 
                           std::make_pair(start_to_null,end_to_null),
                           std::make_pair(start_null_face,end_null_face) );
    assert( vd_checker->face_ok( neg_face ) );
    if (debug) { std::cout << "faces repaired.\n"; }

if (step==current_step) return; current_step++;

    // we are done and can remove split-vertices
    BOOST_FOREACH(HEFace f, incident_faces) {
        remove_split_vertex(f);
    }
    reset_status();


    assert( vd_checker->face_ok( start_null_face ) );
    assert( vd_checker->face_ok( end_null_face ) );
    assert( vd_checker->face_ok( pos_face ) );
    assert( vd_checker->face_ok( neg_face ) );    
    assert( vd_checker->is_valid() );

    //return true; 

}

/// \brief find vertex descriptors corresponding to \a idx1 and \a idx2
// given indices idx1 and idx2, return the corresponding vertex descriptors
// the vertex_map is populated in insert_point_site()
std::pair<HEVertex,HEVertex> VoronoiDiagram::find_endpoints(int idx1, int idx2) {
    std::map<int,HEVertex>::iterator it_start, it_end;
    it_start = vertex_map.find(idx1);
    it_end = vertex_map.find(idx2);
    assert( it_start != vertex_map.end() && it_end != vertex_map.end() ); // we must find idx1 and idx2 in the map
    return std::make_pair( it_start->second, it_end->second);
}


/// \brief prepare null-face 
// next_edge lies on an existing null face
// - Here we either insert a NEW NORMAL or SEPPOINT in the edge
// - OR we push and convert an existing vertex.
// the next_prev flag indicates if we are dealing with the next edge from the new segment-endpoint next_prev=true
// or if we are dealing with the previous edge (next_prev=false)
std::pair<HEVertex,HEFace> VoronoiDiagram::process_null_edge(Point dir, HEEdge next_edge , bool k3, bool next_prev) {
    assert( g[next_edge].type == NULLEDGE );
    HEVertex trg = g.target(next_edge);
    HEVertex src = g.source(next_edge);
    
    HEVertex adj = next_prev ? trg : src; // this is the vertex adjacent to the end-point, on the null face
    assert( g[ next_prev ? src : trg ].type == ENDPOINT );  // this is the end-point of the segment we are inserting
    
    HEVertex sep_point = HEVertex();
    double dir_mult = next_prev ? +1 : -1;
    Point   sep_dir = dir.xy_perp()*dir_mult;
    double sep_alfa = numeric::diangle(sep_dir.x,sep_dir.y); // alfa of potential SEPPOINT

    double new_k3; // k3-value of a new or pushed vertex
    if (next_prev) {
        new_k3 = k3 ? +1 : -1; 
    } else {
        new_k3 = k3 ? -1 : +1;
    }
    
    if (debug) { std::cout << "process_null_edge: ";
        g.print_edge(next_edge);
        std::cout << " next_prev=" << next_prev << "\n";  
    }
    
    if ( g[adj].type == ENDPOINT ) { // target is endpoint
        // insert a normal vertex, positioned at mid-alfa between src/trg.
        HEVertex new_v = g.add_vertex( VoronoiVertex(g[src].position,NEW,NORMAL,g[src].position) );
        double mid = numeric::diangle_mid( g[src].alfa, g[trg].alfa  );
        g[new_v].alfa = mid;
        modified_vertices.insert(new_v);
        g.add_vertex_in_edge( new_v, next_edge);
        g[new_v].k3=new_k3;

        if (debug) {
            std::cout << " e.trg=(ENDPOINT) \n";
            std::cout << " added NEW NORMAL vertex " << g[new_v].index << " in edge "; g.print_edge(next_edge);
        }
        return std::make_pair(HEVertex(), g.HFace() );
        
    } else {
        // Not an ENDPOINT vertex.
        assert( g[adj].type != ENDPOINT );
        double mid(0);
        bool seppoint_pred(false);
        bool parallel_pred(false);
        
        // set the two predicates 
        if (next_prev) {
            HEEdge next_next = g[next_edge].next; 
            HEEdge next_previous = g.previous_edge(next_edge);
            HEVertex next_trg = g.target(next_next); // source
            mid = numeric::diangle_mid( g[src].alfa, g[next_trg].alfa  ); // prev_src, trg
            seppoint_pred = ( g[next_trg].type != ENDPOINT );
            if (debug) {
                std::cout << "          next edge : "; g.print_edge(next_edge);
                std::cout << " next_previous edge : "; g.print_edge(next_previous);
            }
            HEVertex next_out_trg;
            bool next_out_found = null_vertex_target( g.target(next_edge) , next_out_trg ); 
            HEVertex prev_out_trg;
            bool prev_out_found = null_vertex_target( g.source(next_previous), prev_out_trg ); 
            if (next_out_found && prev_out_found) {
                
                if (debug) std::cout << " " << g[g.target(next_edge)].index << " has out-vertex " << g[next_out_trg].index << " status=" << g[next_out_trg].status << "\n";
                if (debug) std::cout << " " << g[g.source(next_previous)].index << " has out-vertex " << g[prev_out_trg].index << " status=" << g[prev_out_trg].status << "\n";
                
                parallel_pred = ( ( ( g[ next_out_trg ].status == OUT ) || ( g[ next_out_trg ].status == NEW ) || ( g[ next_out_trg ].status == UNDECIDED ) ) &&
                                  ( ( g[ prev_out_trg ].status == OUT ) || ( g[ prev_out_trg ].status == NEW ) || ( g[ prev_out_trg ].status == UNDECIDED ))
                                );
            }
        } else { // !next_prev
            HEEdge prev_prev = g.previous_edge(next_edge);
            HEEdge next_next2 = g[next_edge].next;
            HEVertex prev_src = g.source(prev_prev);
            mid = numeric::diangle_mid( g[prev_src].alfa, g[trg].alfa  );
            seppoint_pred = ( g[prev_src].type != ENDPOINT );

            HEVertex next_out_trg2;
            bool next_out_found2 = null_vertex_target( g.source(next_edge), next_out_trg2 ); 
            HEVertex prev_out_trg2;
            bool prev_out_found2 = null_vertex_target( g.target(next_next2), prev_out_trg2 );
            if (next_out_found2 && prev_out_found2) { 
                if (debug) std::cout << " " << g[g.source(next_edge)].index << " has out-vertex " << g[next_out_trg2].index << " status=" << g[next_out_trg2].status << "\n";
                if (debug) std::cout << " " << g[g.target(next_next2)].index << " has out-vertex " << g[prev_out_trg2].index << " status=" << g[prev_out_trg2].status << "\n";
      
                parallel_pred = ( ( ( g[ next_out_trg2 ].status == OUT ) || ( g[ next_out_trg2 ].status == NEW ) || ( g[ next_out_trg2 ].status == UNDECIDED ) ) &&
                                  ( ( g[ prev_out_trg2 ].status == OUT ) || ( g[ prev_out_trg2 ].status == NEW ) || ( g[ prev_out_trg2 ].status == UNDECIDED ) )
                                );
            }

        } // predicates now set.

        // parallel line-segments case
        if ( parallel_pred && g[adj].type == SEPPOINT ) {
            if (debug)  { std::cout << " identical SEPPOINT case!\n";
                //std::cout << " old alfa pred " << (sep_alfa == g[adj].alfa) << "\n";
            }
            // assign face of separator-edge
            // mark separator target NEW
            HEEdge sep_edge=HEEdge();
            BOOST_FOREACH(HEEdge e, g.out_edge_itr(adj) ) {
                assert( g.source(e) == adj );
                if ( g[e].type == SEPARATOR )
                    sep_edge = e;
            }
            assert(sep_edge!=HEEdge()); 
            if (debug) { std::cout << " existing SEPARATOR is "; g.print_edge(sep_edge); }
            HEEdge sep_twin = g[sep_edge].twin;
            HEFace sep_face =  g[sep_edge].face;
            Site* sep_site = g[sep_face].site;
            HEFace sep_twin_face = g[sep_twin].face;
            Site* sep_twin_site = g[sep_twin_face].site;
            HEEdge pointsite_edge;
            if (sep_site->isPoint() ) {
                if (debug) { std::cout << " PointSite SEPARATOR is "; g.print_edge(sep_edge); }
                pointsite_edge = sep_edge;
            }  else if (sep_twin_site->isPoint() ) {
                if (debug) { std::cout << " PointSite SEPARATOR is "; g.print_edge(sep_twin); }
                pointsite_edge = sep_twin;
            }
            
            // set the separator target to NEW
            if (debug) { std::cout << " setting SEPARATOR target "<< g[g.target(sep_edge)].index << " to NEW\n"; }
            HEVertex sep_target = g.target(sep_edge);
            g[sep_target].status = NEW;
            g[sep_target].k3 = new_k3;
            modified_vertices.insert(sep_target);
            
            return std::make_pair( HEVertex(), g[pointsite_edge].face ); // no new separator-point returned
        }

        HEVertex adj_out;
        bool adj_out_found = null_vertex_target(adj, adj_out);
        if (!adj_out_found) { exit(-1); }
        if (debug) std::cout << " not identical SEPPOINT case. either inserting new SEPPOINT, or pushing existing vertex which becomes SEPPOINT/NORMAL\n";
        if (debug) std::cout << " " << g[adj].index << " has out-vertex " << g[adj_out].index << " status=" << g[adj_out].status << "\n";
        
        if ( g[adj_out].status == OUT || g[adj_out].status == UNDECIDED) {
            if (debug) {
                std::cout << " inserting SEPPOINT in edge: "; g.print_edge(next_edge);
                std::cout << "   src alfa = " << g[src].alfa << "\n";
                std::cout << "   SEP==src alfa? = " << (g[src].alfa == sep_alfa) << "\n";
                std::cout << "   SEP alfa = " << sep_alfa << "\n";
                std::cout << "   trg alfa = " << g[trg].alfa << "\n"; 
            }
            sep_point = add_separator_vertex(src, next_edge, sep_dir);
            g[sep_point].k3 = new_k3;
            return std::make_pair( sep_point, g.HFace() );
        } else {
            // target is not endpoint so we push and convert the vertex

            if ( seppoint_pred  ) { 
                // the pushed vertex becomes a SEPPOINT
                if (debug) {
                    std::cout << " pushed vertex " << g[adj].index << " becomes SEPPOINT, ";
                    std::cout << "   sep_alfa = " << sep_alfa << "\n";
                }
                g[adj].alfa = sep_alfa;
                g[adj].type = SEPPOINT;
                g[adj].status = NEW;
                sep_point = adj;
            } else {
                // otherwise it becomes a normal NEW vertex
                if (debug) std::cout << " pushed vertex " << g[adj].index << " becomes NORMAL\n";
                g[adj].alfa = mid;
                g[adj].type = NORMAL;
                g[adj].status = NEW;
            }
            g[adj].k3 = new_k3;
            modified_vertices.insert(adj);
            return std::make_pair( sep_point, g.HFace() );
        }
    }
}

/// \brief add a ::SEPPOINT vertex into the given null-edge
///
/// \param endp the endpoint corresponding to the null-edge/null-face
/// \param edge the null-edge into which we insert the new vertex
/// \param sep_dir direction for setting alfa of the new vertex
HEVertex VoronoiDiagram::add_separator_vertex(HEVertex endp, HEEdge edge, Point sep_dir) {
    HEVertex sep = g.add_vertex( VoronoiVertex(g[endp].position,OUT,SEPPOINT) );
    g[sep].set_alfa(sep_dir);
    if (debug) {
        std::cout << " adding separator " << g[sep].index << " in null edge "; 
        g.print_edge(edge);
    }
    g.add_vertex_in_edge(sep,edge);
    modified_vertices.insert(sep);
    return sep;
}

/// \brief find an adjacent vertex to v, along an edge that is not a ::NULLEDGE
// return true if found, otherwise false.
bool VoronoiDiagram::null_vertex_target( HEVertex v , HEVertex& trg) {
    BOOST_FOREACH( HEEdge e, g.out_edge_itr( v ) ) { 
        if ( g[e].type != NULLEDGE ) {
            trg = g.target(e);
            return true;
        }
    }
    return false;
}


/// \brief either find an existing null-face, or create a new one.
/// \param start  the end of the segment for which we will find/create a null-face
/// \param other  the other end of the new segment
/// \param left   a point left of the new segment
/// \param dir    alfa-direction for positioning endpoint vertex on null-face
/// \param new_site    the new Site we are inserting 
///
/// \return HEVertex ::ENDPOINT-vertex for the new vertex
/// \return HEFace  null-face at endpoint (new or existing)
/// \return HEVertex positive ::SEPARATOR edge endpoint vertex (if a positive separator should be added)
/// \return HEVertex negative ::SEPARATOR edge endpoint vertex (if a negative separator should be added)
/// \return HEFace face-to-null. if a PointSite face should disappear, we return it here.
boost::tuple<HEVertex,HEFace,HEVertex,HEVertex,HEFace>
VoronoiDiagram::find_null_face(HEVertex start, HEVertex other, Point left, Point dir, Site* new_site) {
    HEVertex seg_start = HEVertex(); // new end-point vertices
    HEFace start_null_face; // either existing or new null-face at start-vertex
    
    HEVertex pos_sep_start = HEVertex(); // optional separator endpoints at start
    HEVertex neg_sep_start = HEVertex(); // invalid vertices are default
    
    HEFace face_to_null = g.HFace(); // invalid face is default
    
    // this works for LineSite
    bool k3_sign(true);
    if (new_site->isLine() ) {
        k3_sign = left.is_right( g[start].position , g[other].position); 
    } else if (new_site->isArc()) {
        k3_sign = Point(new_site->x(),new_site->y()).is_right( g[start].position , g[other].position );
    } else {
        assert(0);
    }
    // this is used below and in find_null_face()
    // k3_sign is already calculated in insert_line_segment() ??
    
    if (g[start].null_face != g.HFace() ) { // there is an existing null face
        if (debug) {
            std::cout << " endp= " << g[start].index << " has existing null_face :\n";
            g.print_face( g[start].null_face  );
        }
        start_null_face = g[start].null_face;

        // create a new segment ENDPOINT vertex with zero clearance-disk
        seg_start = g.add_vertex( VoronoiVertex(g[start].position,OUT,ENDPOINT,0) );
        // find the edge on the null-face where we insert seg_start
        HEEdge insert_edge = HEEdge();
        {
            HEEdge current2 = g[start_null_face].edge;
            HEEdge start_edge2 = current2;
            g[seg_start].set_alfa(dir);
            bool found = false;
            if (debug) std::cout << " Looking for endpoint edge:\n";
            do {
                bool face_incident = ( g[ g[ g[current2].twin ].face ].status == INCIDENT);
                if (debug) {
                    std::cout << "  incident= " << face_incident << "  "; g.print_edge(current2);
                }
                if ( face_incident ) { // pick any incident face!
                        insert_edge = current2;
                        found = true;
                }
                current2 = g[current2].next;
            } while (current2!=start_edge2 && !found ); // FIXME end early with !found 
            assert( insert_edge != HEEdge() );
            assert( found );
            /*
            if (!found) {
                std::cout << "find_null_face() FATAL ERROR. can't find edge to insert segment endpoint.\n";
                std::cout << " null face for endpoint " << g[start].index <<" is "; g.print_face(start_null_face);
                
                exit(-1);
            }*/
            if (debug) { std::cout << "  endpoint edge is "; g.print_edge(insert_edge); }
        }
        g.add_vertex_in_edge(seg_start,insert_edge); // insert endpoint in null-edge

        if (debug) { std::cout << "  new endpoint vertex " << g[seg_start].index << " inserted in edge "; g.print_edge(insert_edge); }

        // "process" the adjacent null-edges 
        HEEdge next_edge, prev_edge;
        boost::tie(next_edge,prev_edge) = g.find_next_prev(start_null_face, seg_start);
        //assert( g[prev_edge].next == next_edge );
        boost::tie( neg_sep_start, face_to_null ) = process_null_edge(dir, next_edge, k3_sign, true);
        boost::tie( pos_sep_start, face_to_null ) = process_null_edge(dir, prev_edge, k3_sign, false);
        if (debug) {
            //std::cout << "find_null_face() endp= " << g[start].index << " has existing null_face : " << g[start].null_face << "\n";
            g.print_face( g[start].null_face  );
        }
        return boost::make_tuple( seg_start, start_null_face, pos_sep_start, neg_sep_start, face_to_null);
    } else { // no existing null-face
        // create a new null face at start. the face has three vertices/edges:
        //
        //  neg_sep -> seg_endp -> pos_sep
        //
        start_null_face = g.add_face(); //  this face to the left of start->end edge  
        g[start_null_face].null = true;
          
        if (debug) std::cout << " find_null_face() endp= " << g[start].index <<  " creating new null_face " << start_null_face << "\n";
        seg_start = g.add_vertex( VoronoiVertex(g[start].position,OUT,ENDPOINT) );
        g[seg_start].zero_dist();
        g[seg_start].set_alfa(dir);
        g[seg_start].k3=0;
        pos_sep_start = g.add_vertex( VoronoiVertex(g[start].position,UNDECIDED,SEPPOINT) );
        neg_sep_start = g.add_vertex( VoronoiVertex(g[start].position,UNDECIDED,SEPPOINT) );
        
        g[pos_sep_start].zero_dist();
        g[neg_sep_start].zero_dist();
    
        if (k3_sign) {
            g[pos_sep_start].k3 = +1; 
            g[neg_sep_start].k3 = -1;
        } else {
            g[pos_sep_start].k3 = -1; 
            g[neg_sep_start].k3 = +1;
        }
        g[pos_sep_start].set_alfa( dir.xy_perp()*(+1) );        
        g[neg_sep_start].set_alfa( dir.xy_perp()*(-1) );
        if (debug) {
            std::cout << " k3_sign = " << k3_sign <<"\n"; 
            std::cout << " sep1 = " << g[pos_sep_start].index << " k3=" << g[pos_sep_start].k3 << "\n";
            std::cout << " sep2 = " << g[neg_sep_start].index << " k3=" << g[neg_sep_start].k3 << "\n";
        }
        // null-edges around the face
        HEEdge e1,e1_tw;
        boost::tie(e1,e1_tw) = g.add_twin_edges(seg_start,pos_sep_start);
        HEEdge e2,e2_tw;
        boost::tie(e2,e2_tw) = g.add_twin_edges(pos_sep_start,neg_sep_start);
        HEEdge e3,e3_tw;
        boost::tie(e3,e3_tw) = g.add_twin_edges(neg_sep_start,seg_start);
        
        // e1  ->  e2  ->  e3     on start_null_face, k=1
        // e1t <-  e2t <-  e3t   on g[start].face, k=1
        g.set_next_cycle( boost::assign::list_of(e1)(e2)(e3) , start_null_face, 1);
        HEFace start_face = g[start].face;
        HEEdge start_face_edge=g[start_face].edge; // crazy workaround, because set_next_cycles sets g[face].edge wrong here!
        g.set_next_cycle( boost::assign::list_of(e3_tw)(e2_tw)(e1_tw), g[start].face, 1);
        g[start_null_face].edge = e1;
        g[start_face].edge = start_face_edge;
        
        g[e1].type = NULLEDGE; g[e2].type = NULLEDGE; g[e3].type = NULLEDGE;
        g[e1_tw].type = NULLEDGE; g[e3_tw].type = NULLEDGE; g[e2_tw].type = NULLEDGE;

        g[start].null_face = start_null_face;
        g[start_null_face].site = g[start_face].site;
        
        return boost::make_tuple( seg_start, start_null_face, pos_sep_start, neg_sep_start, face_to_null);
    }
}

/// \brief add ::SEPARATOR edge on the face f, which contains the endpoint
/// \param f is the face of endp 
/// \param null_face null face of the endpoint
/// \param target target-data found by ??
/// \param sep_endp ??
/// \param s1 positive LineSite
/// \param s2 negative LineSite
void VoronoiDiagram::add_separator(HEFace f, HEFace null_face, 
                                   boost::tuple<HEEdge, HEVertex, HEEdge, bool> target,
                                   HEVertex sep_endp, Site* s1, Site* s2) {
    if ( sep_endp == HEVertex() ) // no separator
        return; // do nothing!
    
    if (debug) std::cout << "add_separator() f="<<f<<" endp=" << g[sep_endp].index << "\n";
    
    assert( (g[sep_endp].k3==1) || (g[sep_endp].k3==-1) );    
    g[sep_endp].zero_dist();
    
    HEEdge endp_next = HEEdge();
    HEEdge endp_prev = HEEdge();
    HEEdge endp_next_tw = HEEdge();
    HEEdge endp_prev_tw = HEEdge();

    boost::tie( endp_next_tw, endp_prev_tw ) = g.find_next_prev(null_face, sep_endp);
    if (debug) {
        std::cout << "  add_separator() endp_next_tw="; g.print_edge(endp_next_tw); 
        std::cout << "  add_separator() endp_prev_tw="; g.print_edge(endp_prev_tw); 
        std::cout << "  add_separator() g[endp_next_tw].twin="; g.print_edge(g[endp_next_tw].twin); 
        std::cout << "  add_separator() g[endp_prev_tw].twin="; g.print_edge(g[endp_prev_tw].twin); 
    }
    endp_prev = g[endp_next_tw].twin; // NOTE twin!
    endp_next = g[endp_prev_tw].twin; // NOTE twin!
    assert( endp_next != HEEdge() );
    assert( endp_prev != HEEdge() );

    // find NEW vertex on the old face f
    // this vertex has the correct alfa angle for this endp/separator
    HEEdge v_previous = boost::get<0>(target);
    HEVertex v_target = boost::get<1>(target);
    HEEdge    v_next  = boost::get<2>(target);
    bool    out_new_in = boost::get<3>(target);
    if (debug) {
        std::cout << "  add_separator() v_previous="; g.print_edge(v_previous) ;
        std::cout << "  add_separator() v_target="<< g[v_target].index <<"\n";
    }
    assert( (g[v_target].k3==1) || (g[v_target].k3==-1) );    
    assert( g[sep_endp].k3 == g[v_target].k3 );
    if (!s1->in_region( g[v_target].position ) ) {
        std::cout << " Error separator endpoint not in region of site " << s1->str() << "\n";
        std::cout << " site in_region_t = " << s1->in_region_t(  g[v_target].position ) << "\n";
        std::cout << " site in_region_raw = " << s1->in_region_t_raw(  g[v_target].position ) << "\n";
        
    }
    assert( s1->in_region( g[v_target].position ) ); // v1 and v2 should be in the region of the line-site
    assert( s2->in_region( g[v_target].position ) );
    
    // add new separator edge, and its twin
    HEEdge e2, e2_tw;
    boost::tie(e2,e2_tw) = g.add_twin_edges( sep_endp, v_target );
    g[e2].type    = SEPARATOR;
    g[e2_tw].type = SEPARATOR;
    
    // there are two cases. depending on how v_target (NEW) is found:
    // OUT-NEW-IN, when out_new_in = true
    // IN-NEW-OUT, when out_new_in = false
    // here we set the faces, sites, and next-pointers depending on the case
    if ( out_new_in ) {
        g[e2].k    = g[v_target].k3; // e2 is on the segment side
        g[e2_tw].k = +1;             // e2_tw is on the point-site side
        
        g[e2_tw].face = f; // point-site face
        g[e2_tw].null_face = f;
        g[e2_tw].has_null_face = true;
        
        g[f].edge = e2_tw;
        g[endp_prev].k = g[e2].k; // endp_prev is on the line-site side

        if (g[e2].k == -1) { // choose either s1 or s2 as the site
            g[e2].face = s2->face;
            g[s2->face].edge=e2;
            g[endp_prev].face = s2->face;
        } else {
            g[e2].face = s1->face;
            g[s1->face].edge=e2;
            g[endp_prev].face = s1->face;
        }

        g.set_next(v_previous,e2_tw);
        g.set_next(e2_tw, endp_next);

        g[endp_next].face = f;      // the null-edge
        g[endp_next].k = 1;

        //g.set_next(endp_prev,e2); // don't set!

        g.set_next(e2,v_next);
    } else {
        g[e2].k    = +1;             // e2 is on the point-site side
        g[e2_tw].k = g[v_target].k3; // e2_tw is on the segment side
        
        g[e2].face    = f; // point-site face
        g[e2].null_face    = f;
        g[e2].has_null_face = true;
        
        g[f].edge     = e2;
        g[endp_next].k = g[e2_tw].k; // endp_next is on the linesite-side
        if (g[e2_tw].k == -1) {
            g[e2_tw].face = s2->face;
            g[s2->face].edge=e2_tw;
            g[endp_next].face = s2->face;
        } else {
            g[e2_tw].face = s1->face;
            g[s1->face].edge=e2_tw;
            g[endp_next].face = s1->face;
        }
        g.set_next(v_previous,e2_tw);
        //g.set_next(e2_tw,endp_next);
        g[endp_prev].face = f;
        g[endp_prev].k = 1;

        g.set_next(endp_prev, e2);
        g.set_next(e2,v_next);
    }
    g[e2   ].set_sep_parameters( g[sep_endp].position, g[v_target].position );
    g[e2_tw].set_sep_parameters( g[sep_endp].position, g[v_target].position );
        
    if (debug) {
        std::cout << "add_separator(): ";
        std::cout << g[sep_endp].index << " - " << g[v_target].index << ". Done.\n";
    }
    assert( vd_checker->check_edge(e2) );
    assert( vd_checker->check_edge(e2_tw) );
}

/// find amount of clearance-disk violation on all vertices of face f 
/// \return vertex with the largest clearance-disk violation
HEVertex VoronoiDiagram::find_seed_vertex(HEFace f, Site* site)  {
    if (debug) { std::cout << "find_seed_vertex on f=" << f << "\n"; g.print_face(f); }
    double minPred( 0.0 ); 
    HEVertex minimalVertex = HEVertex();
    bool first( true );
    HEEdge current = g[f].edge;
    HEEdge start = current;
    do {        
        HEVertex q = g.target(current);
        if ( (g[q].status != OUT) && (g[q].type == NORMAL) ) {
            double h = g[q].in_circle( site->apex_point( g[q].position ) ); 
            if (debug) { 
                std::cout << g[q].index << " h= " << h << " dist=" << g[q].dist();
                std::cout << "apex="<< site->apex_point( g[q].position ) << "  ";
                std::cout << "\n";  
            }
            if ( first || ( (h<minPred) && (site->in_region(g[q].position) ) ) ) {
                minPred = h;
                minimalVertex = q;
                first = false;
            }
        }
        current = g[current].next;
    } while(current!=start);  
    assert( minPred < 0 );
    return minimalVertex;
}


/// \brief grow the delete-tree of ::IN vertices by "weighted breadth-first search"
///
/// we start at the seed and add vertices with detH<0 provided that:
/// - (C4) v should not be adjacent to two or more IN vertices (this would result in a loop/cycle!)
/// - (C5) for an incident face containing v: v is adjacent to an IN vertex on this face
///
/// C4 and C5 refer to the Sugihara&Iri 1992 "one million" paper.
///  We process ::UNDECIDED vertices adjacent to known ::IN-vertices in a "weighted breadth-first-search" manner
///  where vertices with a large fabs(detH) are processed first, since we assume the in-circle predicate
///  to be more reliable the larger fabs(in_circle()) is.
void VoronoiDiagram::augment_vertex_set(  Site* site ) {
    while( !vertexQueue.empty() ) {
        HEVertex v = HEVertex();
        double h(0);
        boost::tie( v, h ) = vertexQueue.top();
        assert( g.g[v].status == UNDECIDED );
        vertexQueue.pop(); 
        if ( h < 0.0 ) { // try to mark IN if h<0 and passes (C4) and (C5) tests and in_region(). otherwise mark OUT
            if ( predicate_c4(v) || !predicate_c5(v) || !site->in_region(g[v].position) ) {
                g[v].status = OUT; // C4 or C5 violated, so mark OUT
                if (debug) std::cout << g[v].index << " marked OUT (topo): c4="<< predicate_c4(v) << " c5=" << !predicate_c5(v) << " r=" << !site->in_region(g[v].position) << " h=" << h << "\n";
            } else {
                mark_vertex( v,  site); // h<0 and no violations, so mark IN. push adjacent UNDECIDED vertices onto Q.
                if (debug) { 
                    std::cout << g[v].index << " marked IN (in_circle) ( " << h << " )\n";
                    //std::cout << "  in_region?= " << site->in_region(g[v].position);
                    //std::cout << "  in_region_t= "<< site->in_region_t(g[v].position) << "\n";
                    //std::cout << "  in_region_t_raw= "<< (site->in_region_t_raw(g[v].position)-1.0) << "\n";
                }
            }
        } else {
            g[v].status = OUT; // detH was positive (or zero), so mark OUT
            if (debug) std::cout << g[v].index << " marked OUT (in_circle) ( " << h << " )\n";
        }
        modified_vertices.insert( v );
    }
    
    assert( vertexQueue.empty() );
    assert( vd_checker->all_in(v0) );
    if (debug) std::cout << "augment_vertex_set() DONE\n";

    // sanity-check?: for all incident faces the IN/OUT-vertices should be connected
}

/// mark vertex ::IN and mark adjacent faces ::INCIDENT
// push adjacent UNDECIDED vertices onto queue 
void VoronoiDiagram::mark_vertex(HEVertex& v,  Site* site) {
    g[v].status = IN;
    v0.push_back( v );
    modified_vertices.insert(v);
    
    if (site->isPoint())
        mark_adjacent_faces_p(v);
    else
        mark_adjacent_faces( v, site );

    // push the v-adjacent vertices onto the queue
    BOOST_FOREACH(HEEdge e, g.out_edge_itr( v )) {
        HEVertex w = g.target( e );
        if ( (g[w].status == UNDECIDED) && (!g[w].in_queue) ) {
                // when pushing onto queue we also evaluate in_circle predicate so that we process vertices in the correct order
                vertexQueue.push( VertexDetPair(w , g[w].in_circle(site->apex_point(g[w].position)) ) ); 
                g[w].in_queue=true;
                if (debug) std::cout << "  " << g[w].index << " queued (h=" << g[w].in_circle(site->apex_point(g[w].position)) << " )\n";
        }
    }
}

/// mark adjacent faces ::INCIDENT
// IN-Vertex v has three adjacent faces, mark nonincident faces incident
// and push them to the incident_faces queue
// NOTE: call this only when inserting point-sites
void VoronoiDiagram::mark_adjacent_faces_p( HEVertex v ) {
    assert( g[v].status == IN );
    BOOST_FOREACH(HEEdge e, g.out_edge_itr( v )) {
        HEFace adj_face = g[e].face;
        if ( g[adj_face].status  != INCIDENT ) {
            g[adj_face].status = INCIDENT; 
            incident_faces.push_back(adj_face);
        }
    }

}

/// mark adjacent faces ::INCIDENT
// call this when inserting line-sites
// since we call add_split_vertex we can't use iterators, because they get invalidated
// so use the slower adjacent_faces() instead.
void VoronoiDiagram::mark_adjacent_faces( HEVertex v, Site* site) {
    assert( g[v].status == IN );
    FaceVector new_adjacent_faces = g.adjacent_faces( v );
    
    assert(
        (g[v].type == APEX && new_adjacent_faces.size()==2 ) ||
        (g[v].type == SPLIT && new_adjacent_faces.size()==2 ) ||
        new_adjacent_faces.size()==3
    );

    BOOST_FOREACH( HEFace adj_face, new_adjacent_faces ) {
        if ( g[adj_face].status != INCIDENT ) {
            if ( site->isLine() )
                add_split_vertex(adj_face, site);

            g[adj_face].status = INCIDENT;
            incident_faces.push_back(adj_face);
        }
    }
}

/// \brief find and return edges on which we potentially need ::SPLIT vertices
///
/// walk around the face f
/// return edges whose endpoints are on separate sides of pt1-pt2 line
/// \todo ?not all edges found like this *need* SPLIT vertices? (but it does not hurt to insert SPLIT-vertices in this case)
EdgeVector VoronoiDiagram::find_split_edges(HEFace f, Point pt1, Point pt2) {
    assert( vd_checker->face_ok(f) );
    EdgeVector out;
    HEEdge current_edge = g[f].edge;
    HEEdge start_edge = current_edge;
    //int count=0;                             
    do { // FIND ALL! not just one.
        HEVertex trg = g.target( current_edge );
        HEVertex src = g.source( current_edge );
        bool src_is_right = g[src].position.is_right(pt1,pt2);
        bool trg_is_right = g[trg].position.is_right(pt1,pt2);
        if ( g[src].type == NORMAL || g[src].type == APEX || g[src].type == SPLIT) { //? check edge-type instead?
            if ( src_is_right != trg_is_right  ) { 
                out.push_back(current_edge);
                assert(vd_checker->check_edge(current_edge));
            }
        }
        current_edge = g[current_edge].next;   
        //count++;
        //assert(count<100000); // some reasonable max number of edges in face, to avoid infinite loop
    } while (current_edge!=start_edge);
    
    if (debug) {
        std::cout << " face " << f << " requires SPLIT vertices on edges: \n";
        BOOST_FOREACH( HEEdge e, out ) {
            std::cout << "  "; g.print_edge(e);
        }
    }
    return out;
}

/// \brief add one or many ::SPLIT vertices to the edges of the give face
///
/// these are projections/mirrors of the site of f with the new Site s acting as the mirror
///
/// ::SPLIT vertices are inserted to avoid deleting loops during augment_vertex_set()
void VoronoiDiagram::add_split_vertex(HEFace f, Site* s) {
    assert(!s->isPoint()); // no split-vertices when inserting point-sites
        
    Site* fs = g[f].site;
    
    // don't search for split-vertex on the start or end face
    if (fs->isPoint() && s->isLine()) {
        if ( fs->position() == s->start() || fs->position() == s->end() ) // FIXME: don't compare Points, instead compare vertex-index!
            return;
    }
        
    if ( fs->isPoint() && s->isLine() && s->in_region( fs->position() ) ) {
        // 1) find the correct edge
        Point pt1 = fs->position();
        Point pt2 = pt1-Point( s->a(), s->b() ); 
        
        assert( (pt1-pt2).norm() > 0 ); 
        
        EdgeVector split_edges = find_split_edges(f, pt1, pt2);
        // the sought edge should have src on one side of pt1-pt2
        // and trg on the other side of pt1-pt2
        
        BOOST_FOREACH(HEEdge split_edge, split_edges) {
            if ( (g[split_edge].type == SEPARATOR) || (g[split_edge].type == LINESITE) )
                return; // don't place split points on linesites or separators(?)

            // find a point = src + u*(trg-src)
            // with min_t < u < max_t
            // and minimum distance to the pt1-pt2 line
        #define TOMS748
        
        Point split_pt_pos;
        
        #ifdef TOMS748
            HEVertex split_src = g.source(split_edge);
            HEVertex split_trg = g.target(split_edge);
            if (debug) {
                std::cout << " split src=" << g[split_src].index << "("<< g[split_src].dist() << ")";
                std::cout << " trg=" << g[split_trg].index << "("<< g[split_trg].dist() << ") \n";
                std::cout << "is_right src=" << g[split_src].position.is_right(pt1,pt2) << "  trg="<< g[split_trg].position.is_right(pt1,pt2) << "\n";
            }
            SplitPointError errFunctr( g, split_edge, pt1, pt2); // error functor
            typedef std::pair<double, double> Result;
            boost::uintmax_t max_iter=500;
            boost::math::tools::eps_tolerance<double> tol(64); // bits of tolerance?
            double min_t = std::min( g[split_src].dist() , g[split_trg].dist() );
            double max_t = std::max( g[split_src].dist() , g[split_trg].dist() );
            // require that min_t and max_t bracket the root
            if ( errFunctr(min_t)*errFunctr(max_t) >= 0 )
                return;
                
            Result r1 = boost::math::tools::toms748_solve(errFunctr, min_t, max_t, tol, max_iter);
            split_pt_pos = g[split_edge].point( r1.first ); 
        #endif
        
            // alternative SPLIT-vertex positioning:
            // - create virtual line-site vs: same direction as s(lineSite), but goes through fs(pointSite)
            // - use solver to position SPLIT vertex. The sites are: (vs,fs, fs-adjacent)
        #ifndef TOMS748
            Site* vs = new LineSite(*s);
            vs->set_c( fs->position() ); // modify the line-equation so that the line goes trough fs->position()
            Solution sl = vpos->position( split_edge, vs );
            split_pt_pos = sl.p;
        #endif
        
            HEVertex v = g.add_vertex( VoronoiVertex(split_pt_pos, UNDECIDED, SPLIT, fs->position() ) );
            
        #ifndef TOMS748
            delete vs;
        #endif
        
            //std::cout << "toms748: " << split_pt << "\n";
            //std::cout << "solver:  " << sl.p << "\n";
            if (debug) {
                std::cout << " new split-vertex " << g[v].index << " t=" << r1.first;
                std::cout << " inserted into edge " << g[split_src].index << "-" << g[split_trg].index  << "\n";
            }
            
            assert( vd_checker->check_edge(split_edge) );
            // 3) insert new SPLIT vertex into the edge
            g.add_vertex_in_edge(v, split_edge);
        }
    }
}

/// find a ::SPLIT vertex on the Face f
// return true, and set v, if found.
bool VoronoiDiagram::find_split_vertex(HEFace f, HEVertex& v)  {
    VertexVector verts = g.face_vertices(f);
    BOOST_FOREACH(HEVertex q, verts) {
        if (g[q].type == SPLIT) {
            v = q;
            return true;
        }
    }
    return false;
}

/// \brief remove all ::SPLIT type vertices on the HEFace \a f
void VoronoiDiagram::remove_split_vertex(HEFace f) {

    if (debug) {
        std::cout << "remove_split_vertex( " << f << " )\n";
        g.print_face(f);
    }
    assert( vd_checker->face_ok( f ) );
    
    HEVertex v;
    while ( find_split_vertex(f,v) ) {
        assert(g[v].type == SPLIT); 
        if (debug) std::cout << " removing split-vertex " << g[v].index << "\n";
        
        g.remove_deg2_vertex( v );
        modified_vertices.erase(v);
        
        assert( vd_checker->face_ok( f ) );
    }
    
    assert( vd_checker->face_ok( f ) );
}

/// \brief add ::NEW vertices on ::IN-::OUT edges
/// 
/// generate new voronoi-vertices on all IN-OUT edges 
/// Note: used only by insert_point_site() !!
void VoronoiDiagram::add_vertices( Site* new_site ) {
    if (debug) std::cout << "add_vertices(): \n";
    assert( !v0.empty() );
    EdgeVector q_edges = find_in_out_edges();       // new vertices generated on these IN-OUT edges
    for( unsigned int m=0; m<q_edges.size(); ++m )  {   
        if (debug) {
            HEVertex src = g.source(q_edges[m]);
            HEVertex trg = g.target(q_edges[m]);
            std::cout << " Position NEW vertex on " << g[src].index << " - " << g[trg].index << "\n";
            vpos->solver_debug(true);
        }
        solvers::Solution sl = vpos->position( q_edges[m], new_site ); // vertex_positioner.cpp

        if ( vpos->dist_error( q_edges[m], sl, new_site) > 1e-3 ) {
            HEVertex src = g.source(q_edges[m]);
            HEVertex trg = g.target(q_edges[m]);
            std::cout << "ERROR while positioning new vertex  on edge\n";
            std::cout << g[ src ].index << "[" << g[ src ].type << "]" << "{" << g[ src ].status << "}" << "(t=" << g[ src ].dist() << ")";
            std::cout <<  " -[" << g[q_edges[m]].type << "]- "; 
            std::cout << g[ trg ].index << "[" << g[ trg ].type << "]" << "{" << g[ trg ].status << "}" << "(t=" << g[ trg ].dist() << ")";
            
            std::cout <<  "     derr =" << vpos->dist_error( q_edges[m], sl, new_site) << "\n";
            //exit(-1);
        }
        HEVertex q = g.add_vertex( VoronoiVertex( sl.p, NEW, NORMAL, new_site->apex_point( sl.p ), sl.k3 ) );
        modified_vertices.insert(q);
        g.add_vertex_in_edge( q, q_edges[m] );
        g[q].max_error = vpos->dist_error( q_edges[m], sl, new_site);
        if (debug) {
            HEVertex src = g.source(q_edges[m]);
            HEVertex trg = g.target(q_edges[m]);
            std::cout << " NEW vertex " << g[q].index << " k3= "<< g[q].k3 << " on edge " << g[src].index << " - " << g[trg].index << "\n";
            assert( (g[q].k3==1) || (g[q].k3==-1) );
        }
    }
    if (debug) std::cout << "add_vertices() done.\n";
}

/// \brief add a new face corresponding to the new Site
///
/// call add_new_edge() on all the incident_faces that should be split
HEFace VoronoiDiagram::add_face(Site* s) { 
    HEFace newface =  g.add_face(); 
    g[newface].site = s;
    s->face = newface;
    g[newface].status = NONINCIDENT;
    if (s->isPoint() )
        kd_tree->insert( kd_point( s->position(), newface ) );
        //fgrid->add_face( newface, s->position() ); 
    
    return newface;
}

/// two-argument version used by insert_point_site()
void VoronoiDiagram::add_edges(HEFace newface, HEFace f) {
    add_edges( newface, f, g.HFace(), std::make_pair(HEVertex(),HEVertex()) );
}

/// \brief add all ::NEW-::NEW edges
///
/// by adding a NEW-NEW edge, split the face f into one part which is newface, and the other part is the old f
/// for linesegment or arc sites we pass in both the k=+1 face newface and the k=-1 face newface2
/// the segment endpoints are passed to find_edge_data()
void VoronoiDiagram::add_edges(HEFace newface, HEFace f, HEFace newface2, std::pair<HEVertex, HEVertex> segment) {
    int new_count = num_new_vertices(f);
    if (debug) std::cout << " add_edges() on f=" << f << " with " << new_count << " NEW verts.\n";
    assert( new_count > 0 );
    assert( (new_count % 2) == 0 );
    //if ((new_count % 2) != 0) {
    //    std::cout << " add_edges() FATAL ERROR on f=" << f << " with " << new_count << " NEW verts.\n";
    //    exit(-1);
    //}
    int new_pairs = new_count / 2; // we add one NEW-NEW edge for each pair found
    VertexVector startverts; // this holds ed.v1 vertices for edges already added
    for (int m=0;m<new_pairs;m++) {
        EdgeData ed = find_edge_data(f, startverts, segment);
        add_edge( ed, newface, newface2);
        startverts.push_back( ed.v1 );
    }
    if (debug) std::cout << " add_edges() all edges on f=" << f << " added.\n";
}

/// \brief add a new edge to the diagram
// newface = the k=+1 positive offset face
// newface2 = the k=-1 negative offset face
void VoronoiDiagram::add_edge(EdgeData ed, HEFace newface, HEFace newface2) {
    HEEdge new_previous = ed.v1_prv;
    HEVertex new_source = ed.v1;         //-OUT-NEW(v1)-IN-...
    HEEdge twin_next = ed.v1_nxt;
    
    HEEdge twin_previous = ed.v2_prv;
    HEVertex new_target = ed.v2;         // -IN-NEW(v2)-OUT-
    HEEdge new_next = ed.v2_nxt;
        
    HEFace f = ed.f;
    Site* f_site = g[f].site;
    Site* new_site;
    HEFace new_face;
    if ( g[new_source].k3 == 1 ) { // find out if newface or newface2 should be used
        new_site = g[newface].site;
        new_face = newface; 
    } else {
        new_site = g[newface2].site; 
        new_face = newface2;
    }
        
    // both trg and src should be on same side of new site 
    assert( g[new_target].k3 == g[new_source].k3 );

    //                                           f
    // now connect:   new_previous -> new_source -> new_target -> new_next
    // and:              twin_next <- new_source <- new_target <- twin_previous 
    //                                           new_face   

    // check for potential apex-split
    // we need an apex-vertex if the source and target are on different branches of the new quadratic edge
    // we can set the src_sign and trg_sign with is_right where we compare to a line through the apex 
    bool src_sign=true, trg_sign=true;
    if (f_site->isPoint()  &&  new_site->isLine()  ) { // PL or PA
        Point pt2 = f_site->position();
        Point pt1 = new_site->apex_point(pt2); // projection of pt2 onto LineSite or ArcSite
        src_sign = g[new_source].position.is_right( pt1, pt2 );
        trg_sign = g[new_target].position.is_right( pt1, pt2 );
    } else if (f_site->isPoint()  &&   new_site->isArc()  ) {
        Point pt2 = f_site->position();
        // project p2 onto circle
        Point cen = Point( new_site->x(), new_site->y() );
        Point cen_pt2 = pt2 - cen;
        Point pt1 = cen + (new_site->r()/cen_pt2.norm())*cen_pt2;
        src_sign = g[new_source].position.is_right( pt1, pt2 );
        trg_sign = g[new_target].position.is_right( pt1, pt2 );        
    } else if (f_site->isPoint() && new_site->isPoint() ) { // PP
        src_sign = g[new_source].position.is_right( f_site->position(), new_site->position() );
        trg_sign = g[new_target].position.is_right( f_site->position(), new_site->position() );
    } else if (f_site->isLine() && new_site->isLine() )  { // LL
        //  a line-line bisector, sign should not matter because there is no sqrt()
        /*
        std::cout << "add_edge() LL-edge " << g[new_source].index << " - " << g[new_target].index ;
        std::cout << " f_site->k()= " << f_site->k() << " new_site->k()= "<< new_site->k() << "\n";
        std::cout << " f_site <-> src("<< g[new_source].index << ") = " << g[new_source].position.is_right( f_site->start(), f_site->end() ) << " " << g[new_source].position << "\n";
        std::cout << " f_site <-> trg("<< g[new_target].index << ") = " << g[new_target].position.is_right( f_site->start(), f_site->end() ) << " " << g[new_target].position <<  "\n";
        std::cout << " n_site <-> src("<< g[new_source].index << ") = " << g[new_source].position.is_right( new_site->start(), new_site->end() ) << "\n";
        std::cout << " n_site <-> trg("<< g[new_target].index << ") = " << g[new_target].position.is_right( new_site->start(), new_site->end() ) << "\n";
        */
        // this is essentially an in-region test 
        if ( (g[new_source].position != g[new_target].position) && // src and trg are different
              ( g[new_source].position != f_site->start() ) &&  // src/trg is not start or end
              ( g[new_source].position != f_site->end() ) &&
              ( g[new_target].position != f_site->start() ) &&
              ( g[new_target].position != f_site->end() ) &&
              ( (g[new_source].position -f_site->apex_point( g[new_source].position ) ).norm() > 1e-3 ) && // require some distance, 
              ( (g[new_target].position -f_site->apex_point( g[new_target].position ) ).norm() > 1e-3 )  // so that the is_right predicate is accurate
            ) {
                assert( !g[new_source].position.is_right( f_site->start(), f_site->end() ) );
                assert( !g[new_target].position.is_right( f_site->start(), f_site->end() ) );
                assert( !g[new_source].position.is_right( new_site->start(), new_site->end() ) );
                assert( !g[new_target].position.is_right( new_site->start(), new_site->end() ) );
        }
    } else if (f_site->isLine() && new_site->isArc() )  { // LA
        Point pt2 = Point( new_site->x(), new_site->y() );
        Point pt1 = f_site->apex_point(pt2);
        src_sign = g[new_source].position.is_right( pt1, pt2 );
        trg_sign = g[new_target].position.is_right( pt1, pt2 );
        // if one vertex is on a null-face, we cannot trust the sign
        if ( g[new_source].dist() == 0 || g[new_target].dist() == 0 ) {
            if ( g[new_source].dist() > g[new_target].dist() ) {
                src_sign = trg_sign;
            } else {
                trg_sign = src_sign;
            }
        }
    } else { // unhandled case!
        std::cout << " add_edge() WARNING: no code to deremine src_sign and trg_sign!\n";
        std::cout << " add_edge() f_site " << f_site->str() << "\n";
        std::cout << " add_edge() new_site " << new_site->str() << "\n"; // WARNING: no code to deremine src_sign and trg_sign!\n";
        assert(0);
    }
    
    // both src and trg are on the same side of the new site.
    // so no apex-split is required, just add a single edge.
    if ( src_sign == trg_sign ) {  // add a single src-trg edge
        if (debug) {
            std::cout << " add_edge " << g[new_source].index << " - " << g[new_target].index << "\n";
            std::cout << " f= " << f << " new_face= " << new_face << "\n";
            std::cout << " site= " << f_site->str() << " new_site=" << new_site->str() << "\n";
        }
        HEEdge e_new, e_twin;
        boost::tie(e_new,e_twin) = g.add_twin_edges( new_source, new_target );
        g[e_new].next = new_next;
        assert( g[new_next].k == g[new_previous].k );
        g[e_new].k = g[new_next].k; // the next edge is on the same face, so has the correct k-value
        g[e_new].face = f; // src-trg edge has f on its left
        g[new_previous].next = e_new;
        g[f].edge = e_new; 
        g[e_new].set_parameters( f_site, new_site, !src_sign ); 

        g[twin_previous].next = e_twin;
        g[e_twin].next = twin_next;
        g[e_twin].k = g[new_source].k3; 
        g[e_twin].set_parameters( f_site, new_site,  !src_sign ); // new_site, f_site, src_sign 
        g[e_twin].face = new_face; 
        g[new_face].edge = e_twin;

        assert( vd_checker->check_edge(e_new) && vd_checker->check_edge(e_twin) );
        /*
        std::cout << " added edge " << g[new_target].index << "(" << g[new_target].dist() <<")";
        std::cout << " - " << g[new_source].index << "(" << g[new_source].dist() << ")";
        std::cout << " f=" << g[e_new].face << " k=" << g[e_twin].k;
        std::cout << " twf=" << g[e_twin].face << " twk=" << g[e_new].k;
        std::cout << " type=" << g[e_twin].type <<   " \n";
        */
        //std::cout << " k3 target-source: "<<  g[new_target].k3 << " - " << g[new_source].k3 << "\n";
    } else {
        // need to do apex-split, and add two new edges
        //                         f               f  
        //   new_prv -> NEW -- e1 ----> APEX --e2 ---> NEW -> new_nxt
        //   twn_nxt <- NEW <- e1_tw -- APEX <-e2_tw-- NEW <- twn_prv    
        //                       new1/new2         new1/new2
        //   
        HEVertex apex = g.add_vertex( VoronoiVertex(Point(0,0), NEW,APEX) );
        if (debug) std::cout << " add_edge with APEX " << g[new_source].index << " - [" << g[apex].index << "] - " << g[new_target].index << "\n";
        
        HEEdge e1, e1_tw;
        HEEdge e2, e2_tw;
        boost::tie(e1, e1_tw) = g.add_twin_edges( new_source, apex );
        boost::tie(e2, e2_tw) = g.add_twin_edges( apex, new_target );
        g[e1].set_parameters( f_site, new_site, !src_sign);
        g[e2].set_parameters( f_site, new_site, !trg_sign);
        
        assert( g[new_previous].face == f );
        assert( g[new_next].face == g[new_previous].face );
        assert( g[new_next].k == g[new_previous].k );

        // new_previous -> e1 -> e2 -> new_next
        //g.set_next_chain( boost::assign::list_of(new_previous)(e1)(e2)(new_next), f, g[new_next].k );
        g[new_previous].next=e1; g[e1].next=e2; g[e2].next=new_next;
        g[e1].face=f; g[e2].face=f;
        g[e1].k=g[new_next].k; g[e2].k=g[new_next].k;
        g[f].edge=e1;
    // twin edges
        g[e1_tw].set_parameters(new_site, f_site, src_sign);
        g[e2_tw].set_parameters(new_site, f_site, trg_sign);

        assert( g[twin_previous].k == g[twin_next].k );  
        assert( g[twin_previous].face == g[twin_next].face );        
        // twin_prev -> e2_tw -> e1_tw -> twin_next   on new_face 

        //g.set_next_chain( boost::assign::list_of(twin_previous)(e2_tw)(e1_tw)(twin_next) );
        g[twin_previous].next=e2_tw; g[e2_tw].next=e1_tw; g[e1_tw].next=twin_next;
        //g[e1].face=f; g[e2].face=f;
        //g[e1].k=g[new_next].k; g[e2].k=g[new_next].k;

        //, new_face,  g[new_source].k3  );

                
        g[e1_tw].k = g[new_source].k3;
        g[e2_tw].k = g[new_source].k3;
        g[new_face].edge = e1_tw;
        g[e1_tw].face = new_face;
        g[e2_tw].face = new_face;
        
        assert( vd_checker->check_edge(e1) && vd_checker->check_edge(e1_tw) );
        assert( vd_checker->check_edge(e2) && vd_checker->check_edge(e2_tw) );
        
    // position the apex
        double min_t = g[e1].minimum_t(f_site,new_site);
        g[apex].position = g[e1].point(min_t);
        g[apex].init_dist(f_site->apex_point(g[apex].position));
        modified_vertices.insert( apex );
    }
}


/// \brief find the target of a new ::SEPARATOR edge
/// \param f the HEFace on which we search for the target vertex
/// \param endp the end-point of the null-face with the ::SEPARATOR source
///
/// we want to insert a ::SEPARATOR edge (endp, target) , on the give face f.
/// find and return the target vertex to which the new ::SEPARATOR edge should connect
/// also return the adjacent next/prev edges
///
/// flag==true when an ::OUT-::NEW-::IN vertex was found
///
/// flag==false when an ::IN-::NEW-::OUT vertex was found
///
boost::tuple<HEEdge,HEVertex,HEEdge,bool> VoronoiDiagram::find_separator_target(HEFace f, HEVertex endp) {
    if (endp==HEVertex()) // no separator
        return boost::make_tuple( HEEdge(), HEVertex(), HEEdge(), false) ;
    
    HEEdge current_edge = g[f].edge; // start on some edge of the face
    HEEdge start_edge = current_edge;
    bool found = false;
    HEVertex v_target = HEVertex();
    HEEdge v_previous, v_next;
    bool flag(true);
    if (debug) { 
        std::cout << " find_separator_target f=" << f << " endp= " << g[endp].index << "\n";
        g.print_face(f);
    }
    do {
        HEEdge next_edge = g[current_edge].next;
        HEVertex previous_vertex = g.source( current_edge );
        HEVertex current_vertex  = g.target( current_edge );
        HEVertex next_vertex     = g.target( next_edge );
        bool out_new_in = ( ((g[previous_vertex].status == OUT) || (g[previous_vertex].status == UNDECIDED)) && 
                             g[current_vertex].status == NEW && 
                             g[next_vertex].status == IN );
        bool in_new_out = ( g[previous_vertex].status == IN && 
                            g[current_vertex].status == NEW && 
                            (g[next_vertex].status == OUT || (g[next_vertex].status == UNDECIDED)) ); 
        if ( out_new_in || in_new_out ) {
            if (debug) {
                std::cout << "potential OUT/IN-NEW-IN/OUT: " << g[previous_vertex].index << "-" << g[current_vertex].index;
                std::cout << "-" << g[next_vertex].index << "\n";
            }
            if ( (g[endp].k3 == g[current_vertex].k3)  && endp!=current_vertex ) {
                    v_target = current_vertex;
                    v_previous = current_edge;
                    v_next = next_edge;
                    flag = out_new_in ? true : false;
                    found = true;                 
                    if (debug) std::cout << "FOUND!\n";
            }  else {
                if (debug) {
                    std::cout << " g[endp].k3  = " <<  g[endp].k3  <<"\n";
                    std::cout << " g[current_vertex].k3  = " << g[current_vertex].k3  <<"\n";
                    std::cout << " k3 match? = " <<  (g[endp].k3 == g[current_vertex].k3) <<"\n";
                    std::cout << " endp!=current_vertex ? = " <<  (endp!=current_vertex) <<"\n";
                }
            }
        }
        current_edge = g[current_edge].next;   
        //count++;
        //assert(count<10000); // some reasonable max number of edges in face, to avoid infinite loop
    } while (current_edge!=start_edge && !found);
    if (!found) {
        std::cout << "find_separator_target() FATAL ERROR\n";
        std::cout << " find_separator_target Unable to find target vertex on face f=" << f << " endp= " << g[endp].index << "\n";
        std::cout << " looking for OUT-NEW-IN: " << OUT << " - " << NEW << " - " << IN << "\n";
        std::cout << " looking for IN-NEW-OUT: " << IN << " - " << NEW << " - " << OUT << "\n";
        g.print_face(f);
        exit(-1);
    }
    assert(found);

    return boost::make_tuple(v_previous, v_target, v_next, flag);
}

/// \brief find EdgeData for a new edge
///
/// on a face which has ::IN and ::OUT-vertices, find the sequence
/// OUT-OUT-OUT-..-OUT-NEW(v1)-IN-...-IN-NEW(v2)-OUT-OUT
/// and return v1/v2 together with their previous and next edges
/// \param f face on which we search for vertices
/// \param startverts contains NEW-vertices already found, which are not valid for this call to find_edge_data
/// \param segment contains ENDPOINT vertices, when we are inserting a line-segment
/// (these vertices are needed to ensure finding correct points around sites/null-edges)

VoronoiDiagram::EdgeData VoronoiDiagram::find_edge_data(HEFace f, VertexVector startverts, std::pair<HEVertex,HEVertex> segment)  {
    EdgeData ed;
    ed.f = f;
    if (debug) {
        std::cout << "find_edge_data():\n";
        std::cout << " "; g.print_face(f);
    }
    HEEdge current_edge = g[f].edge; // start on some edge of the face
    HEEdge start_edge = current_edge;
    bool found = false;
    //int count=0;    
    if (debug) std::cout << "    finding OUT-NEW-IN vertex: \n";                         
    do { // find OUT-NEW-IN vertices in this loop
        HEEdge next_edge = g[current_edge].next;
        
        HEVertex previous_vertex = g.source( current_edge);
        HEVertex  current_vertex = g.target( current_edge );
        HEVertex     next_vertex = g.target( next_edge );
        bool previous_not_endpoint = (previous_vertex!=segment.first && previous_vertex!=segment.second);
        bool next_is_endpoint = (next_vertex==segment.first || next_vertex==segment.second);
        
        if ( (g[current_vertex].status==NEW) && (g[current_vertex].type != SEPPOINT) &&
             (  ( (g[previous_vertex].status==OUT || g[previous_vertex].status==UNDECIDED)  &&  
                     previous_not_endpoint ) 
                   ||
                ( next_is_endpoint )
             )
           ) {
            // slow? linear search through vector. but startverts.size() should not be too large..
            bool v_in_startverts =
                ( std::find(startverts.begin(), startverts.end(),  current_vertex) != startverts.end() );
            if (debug) {
                std::cout << "     " << g[current_vertex].index << "N=" << (g[current_vertex].status == NEW) ;
                std::cout << " !SEPP=" << (g[current_vertex].type != SEPPOINT) << "\n";
            }
            if ( !v_in_startverts ) {
                ed.v1 = current_vertex;
                ed.v1_prv = current_edge;
                ed.v1_nxt = g[current_edge].next;
                found = true;
            }
        }
        current_edge = g[current_edge].next;   
        //count++;
        //assert(count<10000); // some reasonable max number of edges in face, to avoid infinite loop
    } while (current_edge!=start_edge && !found);
    //if (!found) {
    //    std::cout << "ERROR: unable to find OUT-NEW-IN vertex on face:\n";
    //    g.print_face(f);
    //    std::cout << " The excluded vertices are: (size=" << startverts.size()<<")"; g.print_vertices(startverts);
    //}
    assert(found);
    if (debug) std::cout << " OUT-NEW-IN vertex is  " << g[ed.v1].index << "\n";

    // now search for v2
    //count=0; 
    start_edge = current_edge; // note that this search starts where we ended in the loop above!
    found=false;
    if (debug) std::cout << "    finding IN-NEW-OUT vertex: \n";   
    do { // find IN-NEW-OUT vertices in this loop
        HEVertex  current_vertex = g.target( current_edge );
        if ( g[current_vertex].status == NEW && g[current_vertex].type != SEPPOINT ) {
            if (debug) {
                std::cout << "     " << g[current_vertex].index << "N=" << (g[current_vertex].status == NEW) ;
                std::cout << " !SEPP=" << (g[current_vertex].type != SEPPOINT);
                std::cout << " !ed.v1=" << (current_vertex != ed.v1) <<"\n";
            }
            if (  current_vertex != ed.v1) { // -IN-NEW(v2)
                    ed.v2     = current_vertex;
                    ed.v2_prv = current_edge;
                    ed.v2_nxt = g[current_edge].next;
                    found = true;                 
            }
        }
        current_edge = g[current_edge].next;   
        //count++;
        //assert(count<10000); // some reasonable max number of edges in face, to avoid infinite loop
    } while (current_edge!=start_edge && !found);
    assert(found);
    if (debug) std::cout << " IN-NEW-OUT vertex is " << g[ed.v2].index << "\n";

    if (debug) std::cout << "find_edge_data() NEW-NEW vertex pair: " << g[ed.v1].index << " - " << g[ed.v2].index << "\n";
    return ed;
}

/// one-argument version of repair_face() used by insert_point_site()
void VoronoiDiagram::repair_face( HEFace f ) {
    repair_face(f,  std::make_pair(HEVertex(),HEVertex()), 
                    std::make_pair(g.HFace(), g.HFace() ),
                    std::make_pair(g.HFace(), g.HFace() ) );
}

/// \brief repair next-pointers of HEFace \a f
///
/// start on g[newface].edge, walk around the face and repair the next-pointers
/// this is called on the newly created face after all NEW-NEW edges have been added
void VoronoiDiagram::repair_face( HEFace f, std::pair<HEVertex,HEVertex> segment, 
                                            std::pair<HEFace,HEFace> nulled_faces,
                                            std::pair<HEFace,HEFace> null_face ) {
    if (debug) {
        std::cout << "repair_face ( " << f << " ) null1=" << null_face.first << " null2=" << null_face.second << "\n";
        if ( (segment.first!=HEVertex()) && (segment.second!=HEVertex()) )
            std::cout << " seg_start=" << g[segment.first].index << " seg_end=" << g[segment.second].index << "\n";
        std::cout << " nulled.first=" << nulled_faces.first << " nulled.second=" << nulled_faces.second << "\n";
    }
    HEEdge current_edge = g[f].edge;
    HEEdge start_edge = current_edge;
    //int count=0;
    do {
        assert( vd_checker->check_edge(current_edge) );
        HEVertex current_target = g.target( current_edge ); // an edge on the new face
        HEVertex current_source = g.source( current_edge );
        #ifndef NDEBUG
        bool found_next_edge= false;
        #endif
        if (debug) { 
            std::cout << " edge " << g[ g.source(current_edge) ].index << " - ";
            std::cout <<  g[ g.target(current_edge) ].index << "\n";
        }
        BOOST_FOREACH(HEEdge e, g.out_edge_itr(current_target)){
            HEVertex out_target = g.target( e );
            if(debug) {
                std::cout << "     candidate: " << g[ g.source(e) ].index << " - ";
                std::cout << g[ g.target(e) ].index << " f= "<< g[e].face << " \n";
            }
            if ( (out_target != current_source) && 
                 ( (g[out_target].status == NEW)    || 
                   (g[out_target].type == ENDPOINT) || 
                   (g[out_target].type == SEPPOINT) ) ) { // these are the possible vertices we want to go to
                
                
                // special cases where we do a brute-force face-assignment for a null-edge, or a separator
                if ( ((g[e].type == NULLEDGE) &&
                      (g[current_edge].type != NULLEDGE) && // only one null-edge in succession!
                         (
                           // from sep to end
                           ( (g[current_target].type==SEPPOINT) && (g[out_target].type == ENDPOINT) ) ||
                           // or from end -> end 
                           ( (g[current_source].type == ENDPOINT) && (g[current_target].type==ENDPOINT)  )
                           ||
                           (out_target == segment.first)
                           ||
                           (out_target == segment.second) 
                         ) &&
                        (g[e].face!=null_face.first) && // not along a null-face edge!
                        (g[e].face!=null_face.second)
                     ) 
                     ||
                     (g[e].face == nulled_faces.first) // edge previously belonged to point-site that has disappeared
                     ||
                     (g[e].face == nulled_faces.second)
                     )
                      {
                         
                    g[e].face = f; // override face-assignment!
                    g[e].k=g[current_edge].k; // override k-assignment!
                    if (debug) std::cout << " face and k-val override! f="<< f << " k=" << g[current_edge].k << "\n";
                }
                    
                // the next vertex should not where we came from
                // and it should be on the same face.
                if ( g[e].face == f ) {  
                    g[current_edge].next = e; // this is the edge we want to take
                    #ifndef NDEBUG
                    found_next_edge = true;
                    #endif
                    if(debug) {
                        std::cout << "         next: " << g[ g.source(e) ].index << " - ";
                        std::cout << g[ g.target(e) ].index << "\n";
                    }
                    assert( g[current_edge].k == g[e].k );
                    assert( vd_checker->current_face_equals_next_face( current_edge ) );
                }
            } 
        }
        //if (!found_next_edge) {
        //    std::cout << " repair_face( " << f << " ) error. could not find next-edge!\n";
        //    exit(-1);
        //}
        assert(found_next_edge); // must find a next-edge!
        //count++;
        //if (count>30)
        //    exit(-1);
        current_edge = g[current_edge].next; // jump to the next edge
    } while (current_edge != start_edge);
    
}

/// \brief remove the ::IN vertices of the delete-tree
///
/// removes the IN vertices stored in v0 (and associated IN-NEW edges)
void VoronoiDiagram::remove_vertex_set() {
    BOOST_FOREACH( HEVertex& v, v0 ) {      // it should now be safe to delete all IN vertices
        assert( g[v].status == IN );
        g.delete_vertex(v); // this also removes edges connecting to v
        modified_vertices.erase(v);
    }
}

/// \brief reset status of modified_vertices and incident_faces
///    
/// at the end after an incremental insertion of a new site,
/// reset status of modified_vertices to UNDECIDED and incident_faces to NONINCIDENT,
/// so that we are ready for the next insertion.
void VoronoiDiagram::reset_status() {
    BOOST_FOREACH( HEVertex v, modified_vertices ) {
        g[v].reset_status();
    }
    modified_vertices.clear();
    BOOST_FOREACH(HEFace& f, incident_faces ) { 
        g[f].status = NONINCIDENT; 
    }
    incident_faces.clear();
    v0.clear();
}

/// \brief find and return ::IN - ::OUT edges
/// 
/// given the set v0 of ::IN vertices, find and return the adjacent ::IN - ::OUT edges.
/// Later ::NEW vertices are inserted into each of the found ::IN - ::OUT edges
EdgeVector VoronoiDiagram::find_in_out_edges() { 
    assert( !v0.empty() );
    EdgeVector output; // new vertices generated on these edges
    BOOST_FOREACH( HEVertex& v, v0 ) {                                   
        assert( g[v].status == IN ); // all verts in v0 are IN
        BOOST_FOREACH(HEEdge e, g.out_edge_itr(v)){
            if ( g[ g.target( e ) ].status == OUT ) 
                output.push_back(e); // this is an IN-OUT edge
        }
    }
    if (debug) std::cout << "find_in_out_edges() " << output.size() << " IN-OUT edges \n";
    assert( !output.empty() );
    return output;
}

/// \brief adjacent in-count predicate for buildingdelete-tree
///
/// number of IN vertices adjacent to given vertex v
/// predicate C4 i.e. "adjacent in-count" from Sugihara&Iri 1992 "one million" paper
bool VoronoiDiagram::predicate_c4(HEVertex v) {
    int in_count=0;
    BOOST_FOREACH(HEEdge e, g.out_edge_itr(v)){
        HEVertex w = g.target( e );
        if ( g[w].status == IN ) {
            in_count++;
            if (in_count >= 2)
                return true;
        }
    }
    return false;
}

/// \brief connectedness-predicate for delete-tree building
///
/// do any of the three faces that are adjacent to the given IN-vertex v have an IN-vertex ?
/// predicate C5 i.e. "connectedness"  from Sugihara&Iri 1992 "one million" paper
bool VoronoiDiagram::predicate_c5(HEVertex v) {
    if (g[v].type == APEX || g[v].type == SPLIT ) { return true; } // ?
    FaceVector adjacent_incident_faces;

    BOOST_FOREACH(HEEdge e, g.out_edge_itr(v)){
        //HEFace f = ;
        if ( g[ g[e].face ].status == INCIDENT )
            adjacent_incident_faces.push_back( g[e].face );
    }

    assert( !adjacent_incident_faces.empty() );
    
    //bool all_found = true;
    BOOST_FOREACH( HEFace f, adjacent_incident_faces ) { // check each adjacent face f for an IN-vertex
        bool face_ok=false;
        HEEdge current = g[f].edge;
        HEEdge start = current;
        do {
            HEVertex w = g.target(current);
            if ( w != v && g[w].status == IN && g.has_edge(w,v) )  // v should be adjacent to an IN vertex on the face
                face_ok = true;
            else if ( w!=v && ( g[w].type == ENDPOINT || g[w].type == APEX  || g[w].type == SPLIT) ) // if we are next to an ENDPOINT, then ok(?)
                face_ok=true;
            current = g[current].next;
        } while(current!=start);  

        if (!face_ok)
            return false;
            //all_found=false;
    }
    return true; // if we get here we found all ok
    //return all_found; // if this returns false, we mark a vertex OUT, on topology grounds.
}

/// return number of ::SPLIT vertices
int VoronoiDiagram::num_split_vertices() const { 
    int count = 0;
    BOOST_FOREACH( const HEVertex v, g.vertices() ) {
        if (g[v].type == SPLIT)
            count++;
    }
    return count; 
}

/// count number of ::NEW vertices on the given face \a f
int VoronoiDiagram::num_new_vertices(HEFace f) {
    HEEdge current = g[f].edge;
    HEEdge start = current;
    int count=0;
    //int num_e=0;
    do {
        HEVertex v = g.target(current);
        if ( (g[v].status == NEW) && (g[v].type != SEPPOINT) )
            count++;
        //num_e++;
        //assert( num_e <3000);
        current = g[current].next;
    } while(current!=start);  
    return count;
}

/// filter the graph using given Filter \a flt
void VoronoiDiagram::filter( Filter* flt) {
    flt->set_graph(&g);
    BOOST_FOREACH(HEEdge e, g.edges() ) {
        if ( ! (*flt)(e) )
            g[e].valid = false;
    }
}

/// \brief reset filtering by setting all edges valid
void VoronoiDiagram::filter_reset() { // this sets valid=true for all edges 
    BOOST_FOREACH(HEEdge e, g.edges() ) {
        g[e].valid = true;
    }
}
    
/// run topology/geometry check on diagram
bool VoronoiDiagram::check() {
    if( vd_checker->is_valid() ) {
        if (debug) std::cout << "diagram check OK.\n";
        return true;
    } else {
        if (debug) std::cout << "diagram check ERROR.\n";
        return false;
    }
}

/// string repr
std::string VoronoiDiagram::print() const {
    std::ostringstream o;
    o << "VoronoiDiagram \n";
    o << " num_vertices    = "<< g.num_vertices() << "\n";
    o << " num_edges       = "<< g.num_edges() <<"\n";
    o << " num_point_sites = "<< num_point_sites() <<"\n";
    o << " num_line_sites  = "<< num_line_sites() <<"\n";
    o << " num_split_vertices  = "<< num_split_vertices() <<"\n";
    return o.str();
}

} // end namespace
// end file voronoidiagram.cpp
