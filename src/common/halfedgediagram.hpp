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

#include <vector>
#include <list>

#include <boost/graph/adjacency_list.hpp>
#include <boost/foreach.hpp> 
#include <boost/iterator/iterator_facade.hpp>
#include <boost/assign/list_of.hpp>

// bundled BGL properties, see: http://www.boost.org/doc/libs/1_44_0/libs/graph/doc/bundles.html

// dcel notes from http://www.holmes3d.net/graphics/dcel/

// vertex (boost::out_edges)
//  -leaving pointer to HalfEdge that has this vertex as origin
//   if many HalfEdges have this vertex as origin, choose one arbitrarily

// HalfEdge
//  - origin pointer to vertex (boost::source)
//  - face to the left of halfedge
//  - twin pointer to HalfEdge (on the right of this edge)
//  - next pointer to HalfEdge
//     this edge starts from h->twin->origin and ends at next vertex in h->face
//     traveling ccw around boundary
//     (allows face traverse, follow h->next until we arrive back at h)

// Face
//  - edge pointer to HalfEdge
//    this edge has this Face object as face
//    half-edge can be any one on the boundary of face
// special "infinite face", face on "outside" of boundary
// may or may not store edge pointer





namespace hedi  { 
/*! 
 * \namespace hedi
 * \brief Half-edge diagram
 */

/// \brief half-edge diagram, based on the boost graph-library
///
/// half_edge_diagram is a half-edge diagram class.
/// Templated on Vertex/Edge/Face property classes which allow
/// attaching information to vertices/edges/faces that is 
/// required for a particular algorithm.
/// 
/// Inherits from boost::adjacency_list
/// minor additions allow storing face-properties.
///
/// the hedi namespace contains functions for manipulating HEDIGraphs
///
/// For a general description of the half-edge data structure see e.g.:
///  - http://www.holmes3d.net/graphics/dcel/
///  - http://openmesh.org/index.php?id=228
template <class TOutEdgeList, 
          class TVertexList,
          class TDirected, 
          class TVertexProperties,
          class TEdgeProperties,
          class TFaceProperties,
          class TGraphProperties,
          class TEdgeList 
          >
class half_edge_diagram {
public:
    /// type of face descriptor
    typedef unsigned int Face; 
    /// underlying boost graph type
    typedef typename boost::adjacency_list< TOutEdgeList,            
                                            TVertexList,            
                                            TDirected,   
                                            TVertexProperties,             
                                            TEdgeProperties,                
                                            TGraphProperties,
                                            TEdgeList
                                            > BGLGraph;
    /// edge descriptor
    typedef typename boost::graph_traits< BGLGraph >::edge_descriptor   Edge;
    /// vertex descriptor
    typedef typename boost::graph_traits< BGLGraph >::vertex_descriptor Vertex;
    /// vertex iterator type
    typedef typename boost::graph_traits< BGLGraph >::vertex_iterator   VertexItr;
    /// out edge iterator type
    typedef typename boost::graph_traits< BGLGraph >::out_edge_iterator OutEdgeItr;
    /// edge iterator type
    typedef typename boost::graph_traits< BGLGraph >::edge_iterator     EdgeItr; 

    /// vertex descriptor
    typedef typename boost::graph_traits< BGLGraph >::vertex_descriptor      vertex_descriptor;
    /// edge descriptor
    typedef typename boost::graph_traits< BGLGraph >::edge_descriptor        edge_descriptor;
    /// edge iterator
    typedef typename boost::graph_traits< BGLGraph >::edge_iterator          edge_iterator;
    /// out edge iterator
    typedef typename boost::graph_traits< BGLGraph >::out_edge_iterator      out_edge_iterator;
    /// in edge iterator
    typedef typename boost::graph_traits< BGLGraph >::in_edge_iterator       in_edge_iterator;
    /// vertex iterator
    typedef typename boost::graph_traits< BGLGraph >::vertex_iterator        vertex_iterator;
    /// directed or underected graph
    typedef typename boost::graph_traits< BGLGraph >::directed_category      directed_category;
    /// allow or disallow parallel edges
    typedef typename boost::graph_traits< BGLGraph >::edge_parallel_category edge_parallel_category;
    /// ?
    typedef typename boost::graph_traits< BGLGraph >::traversal_category     traversal_category;
    /// vertex size type
    typedef typename boost::graph_traits< BGLGraph >::vertices_size_type     vertices_size_type;
    /// edge size type
    typedef typename boost::graph_traits< BGLGraph >::edges_size_type        edges_size_type;
    /// degree size type
    typedef typename boost::graph_traits< BGLGraph >::degree_size_type       degree_size_type;
    /// adjacency iterator
    typedef typename boost::graph_traits< BGLGraph >::adjacency_iterator     adjacency_iterator;

    /// vector of vertices
    typedef std::vector<Vertex> VertexVector;
    /// vector of faces
    typedef std::vector<Face>   FaceVector;
    /// vector of edges
    typedef std::vector<Edge>   EdgeVector;  
    
    /// access to Face properties
    inline TFaceProperties& operator[](Face f) { return faces[f]; }
    /// const access to Face properties
    inline const TFaceProperties& operator[](Face f) const { return faces[f]; } 
    /// access to Edge properties
    inline TEdgeProperties& operator[](Edge e) { return g[e]; }
    /// const access to Edge properties
    inline const TEdgeProperties& operator[](Edge e) const { return g[e]; }
    /// access to Vertex properties
    inline TVertexProperties& operator[](Vertex v)  { return g[v]; }
    /// const access to Vertex properties
    inline const TVertexProperties& operator[](Vertex v) const  { return g[v]; }

//DATA
    /// container for face properties
    std::vector< TFaceProperties > faces; // this could maybe be held as a GraphProperty of the BGL-graph?
    /// underlying BGL graph
    BGLGraph g;
    
// NOTE: there is no HEDIGraph constructor, we use the default one..

/// dtor
virtual ~half_edge_diagram(){
    // sites are associated with faces. go through all faces and delete the site
    //std::cout << "~half_edge_diagram()...";
    /*
    BOOST_FOREACH( TFaceProperties fprop, faces ) {
        if (fprop.site)
            delete fprop.site;
    }
    g.clear();
    */
    //std::cout << "DONE.";
}

// One-liner wrappers around boost-graph-library functions:

/// return an invalid face_descriptor
Face HFace() { return std::numeric_limits<Face>::quiet_NaN(); }
/// add a blank vertex and return its descriptor
Vertex add_vertex() { return boost::add_vertex( g ); }
/// add a vertex with given properties, return vertex descriptor
Vertex add_vertex(const TVertexProperties& prop) { return boost::add_vertex( prop, g ); }
/// return the target vertex of the given edge
Vertex target(const Edge e ) const { return boost::target( e, g ); }
/// return the source vertex of the given edge
Vertex source(const Edge e ) const { return boost::source( e, g ); }
/// return degree of given vertex
unsigned int degree( Vertex v)  { return boost::degree( v, g); }
/// return number of faces in graph
unsigned int num_faces() const { return faces.size(); }
/// return number of vertices in graph
unsigned int num_vertices() const { return boost::num_vertices( g ); }
/// return number of edges in graph
unsigned int num_edges() const { return boost::num_edges( g ); }
/// return number of edges on Face f
unsigned int num_edges(Face f) { return face_edges(f).size(); }
/// add an edge between vertices v1-v2
Edge add_edge(Vertex v1, Vertex v2) { return boost::add_edge( v1, v2, g).first; }
/// add an edge with given properties between vertices v1-v2
Edge add_edge( Vertex v1, Vertex  v2, const TEdgeProperties& prop ) { return boost::add_edge( v1, v2, prop, g).first; }
/// return begin/edge iterators for out-edges of Vertex \a v
std::pair<OutEdgeItr, OutEdgeItr> out_edge_itr( Vertex v ) { return boost::out_edges( v, g ); } // FIXME: change name to out_edges!!
/// return true if v1-v2 edge exists
inline bool has_edge( Vertex v1, Vertex v2) { return boost::edge( v1, v2, g ).second; }
/// return v1-v2 Edge
Edge edge( Vertex v1, Vertex v2) { assert(has_edge(v1,v2)); return boost::edge( v1, v2, g ).first; }
/// clear given vertex. this removes all edges connecting to the vertex.
void clear_vertex( Vertex v ) { boost::clear_vertex( v, g ); }
/// remove given vertex. call clear_vertex() before this!
void remove_vertex( Vertex v ) { boost::remove_vertex( v , g ); }
/// remove given edge
void remove_edge( Edge e ) { boost::remove_edge( e , g ); }
/// delete a vertex. clear and remove.
void delete_vertex(Vertex v) { clear_vertex(v); remove_vertex(v); }

/// insert Vertex \a v into the middle of Edge \a e
void add_vertex_in_edge( Vertex v, Edge e) {
    // the vertex v is inserted into the middle of edge e
    // edge e and its twin are replaced by four new edges: e1,e2 and their twins te2,te1
    // before:             face
    //                      e
    // previous-> source  ------> target -> next
    //  tw_next<- tw_trg  <-----  tw_src <- tw_previous
    //                      twin 
    //                    twin_face
    //
    // after:               face
    //                    e1   e2
    // previous-> source  -> v -> target -> next
    //  tw_next<- tw_trg  <- v <- tw_src <- tw_previous
    //                    te2  te1
    //                    twin_face
    //

    Edge e_twin = g[e].twin;
    assert( e_twin != Edge() );
    Vertex esource = boost::source(e,g); 
    Vertex etarget = boost::target(e,g); 
    Face face = g[e].face;
    Face twin_face = g[e_twin].face;
    Edge previous = previous_edge(e);
    Edge twin_previous = previous_edge(e_twin);
    
    assert( g[previous].face == g[e].face );
    assert( g[twin_previous].face == g[e_twin].face );
    
    Edge e1 = boost::add_edge( esource, v, g).first;
    Edge te2 = boost::add_edge( v, esource,  g).first;
    g[e1].twin = te2; g[te2].twin = e1;
    //boost::tie(e1,te2) = add_twin_edges( esource, v ); 
    //Edge e2, te1;
    //boost::tie(e2,te1) = add_twin_edges( v, etarget );    
    Edge e2 = boost::add_edge( v, etarget, g).first;
    Edge te1 = boost::add_edge( etarget, v,  g).first;
    g[e2].twin = te1; g[te1].twin = e2;


    // next-pointers
    g[previous].next = e1; g[e1].next=e2; g[e2].next = g[e].next;
    //set_next_chain( boost::assign::list_of(previous)(e1)(e2)(g[e].next) );
    //set_next_chain( boost::assign::list_of(twin_previous)(te1)(te2)(g[e_twin].next) );
    g[twin_previous].next = te1; g[te1].next=te2; g[te2].next = g[e_twin].next;    
    // this copies params, face, k, type
    g[e1] = g[e];       g[e2] = g[e];       // NOTE: we use EdgeProperties::operator= here to copy !
    g[te1] = g[e_twin]; g[te2] = g[e_twin];
    // update the faces 
    faces[face].edge = e1;
    faces[twin_face].edge = te1;
    // finally, remove the old edge
    //remove_twin_edges(esource, etarget);
    boost::remove_edge( e , g );
    boost::remove_edge( e_twin , g );
}
/// ad two edges, one from \a v1 to \a v2 and one from \a v2 to \a v1
std::pair<Edge,Edge> add_twin_edges(Vertex v1, Vertex v2) {
    //Edge e1,e2;
    //bool b;
    //boost::tie( e1 , b ) = boost::add_edge( v1, v2, g);
    //boost::tie( e2 , b ) = boost::add_edge( v2, v1, g);
    Edge e1 = boost::add_edge( v1, v2, g).first;
    Edge e2 = boost::add_edge( v2, v1, g).first;
    //twin_edges(e1,e2);
    g[e1].twin = e2;
    g[e2].twin = e1;
    return std::make_pair(e1,e2);
}

/// make e1 the twin of e2 (and vice versa)
void twin_edges( Edge e1, Edge e2 ) {
    if (target(e1) != source(e2)) {
        std::cout << " error target(e1)= " << g[target(e1)].index << " != " << g[source(e2)].index << " = source(e2) \n";
        std::cout << "target(e1) = " << target(e1) << "\n";
        std::cout << "source(e2) = " << source(e2) << "\n";
    }
    assert( target(e1) == source(e2) );
    assert( source(e1) == target(e2) );
    
    g[e1].twin = e2;
    g[e2].twin = e1;
}

/// add a face 
Face add_face() {
    TFaceProperties f_prop;
    faces.push_back( f_prop); 
    Face index = faces.size()-1;
    faces[index].idx = index;
    return index;    
}

/// add a face, with given properties
Face add_face(const TFaceProperties& prop) {
    faces.push_back( prop ); 
    Face index = faces.size()-1;
    faces[index].idx = index;
    return index;    
}

/// return all vertices in a vector of vertex descriptors
VertexVector vertices()  const {
    VertexVector vv;
    VertexItr it_begin, it_end, itr;
    boost::tie( it_begin, it_end ) = boost::vertices( g );
    for ( itr=it_begin ; itr != it_end ; ++itr ) {
        vv.push_back( *itr );
    }
    return vv;
}

/// return all vertices adjecent to given vertex
VertexVector adjacent_vertices(  Vertex v) const {
    VertexVector vv;
    BOOST_FOREACH( Edge e, out_edges( v ) ) {
        vv.push_back( target( e ) );
    }
    return vv;
}

/// return all vertices of given face
VertexVector face_vertices(Face face_idx) const {
    VertexVector verts;
    Edge startedge = faces[face_idx].edge; // the edge where we start
    Vertex start_target = boost::target( startedge, g ); 
    verts.push_back(start_target);
    Edge current = g[startedge].next;
    int count=0;
    EdgeVector f_edges; // for debug.
    f_edges.push_back(current);
    do {
        Vertex current_target = boost::target( current, g ); 
        //assert( current_target != start_target );
        verts.push_back(current_target);
        f_edges.push_back(current);
        assert( g[current].face == g[ g[current].next ].face );
        current = g[current].next;
        
        if (count >= 3000000 ) {
            std::cout << " ERROR too many vertices on face! count=" << count << "\n";
            std::cout << " verts.size() = " << verts.size();
            std::cout << " edges.size()=" << f_edges.size() <<"\n";
            for (unsigned int n=0;n<verts.size()-10;n++) {
                std::cout << n << "   : " << g[ verts[n] ].index << "\n"; 
            }
        }
        assert( count < 3000000 ); // stop at some max limit
        count++;
    } while ( current != startedge );
    return verts;
}

/// return edges of face f as a vector
/// NOTE: it is faster to write a do-while loop in client code than to call this function!
EdgeVector face_edges( Face f) const {
    Edge start_edge = faces[f].edge;
    Edge current_edge = start_edge;
    EdgeVector out;
    std::cout << " edges on face " << f << " :\n ";
    do {
        Vertex src = source(current_edge);
        Vertex trg = target(current_edge);
           std::cout << out.size() << " " << g[src].index << "[" << g[src].type <<"]";
           std::cout << " - " << g[trg].index << "[" << g[trg].type <<"]" <<"\n ";
        out.push_back(current_edge);
        current_edge = g[current_edge].next;
    } while( current_edge != start_edge );
    return out;
}

/// return out_edges of given vertex
EdgeVector out_edges( Vertex v) const { 
    EdgeVector ev;
    OutEdgeItr it, it_end;
    boost::tie( it, it_end ) = boost::out_edges( v, g );
    for ( ; it != it_end ; ++it ) {
        ev.push_back(*it);
    }
    return ev;
}

/// return all edges as a vector
// FIXME: provide std::pair<edge_iterator,edge_iterator> version of this function also?
EdgeVector edges() const {
    EdgeVector ev;
    EdgeItr it, it_end;
    boost::tie( it, it_end ) = boost::edges( g );
    for ( ; it != it_end ; ++it ) {
        ev.push_back(*it);
    }
    return ev;
}

/// return the previous edge. traverses all edges in face until previous found.
Edge previous_edge( Edge e ) const {
    Edge previous = g[e].next;
    while ( g[previous].next != e ) {
        previous = g[previous].next;
    }
    return previous;
}

/// return adjacent faces to the given vertex
FaceVector adjacent_faces( Vertex q ) const {
    std::set<Face> face_set;
    OutEdgeItr itr, itr_end;
    boost::tie( itr, itr_end) = boost::out_edges(q, g);
    for ( ; itr!=itr_end ; ++itr ) {
        face_set.insert( g[*itr].face );
    }
    //assert( face_set.size() == 3); // true for normal vertices, but SPLIT/APEX are degree 2..
    FaceVector fv(face_set.begin(), face_set.end());
    return fv;
}

/// inserts given vertex into edge e, and into the twin edge e_twin
/// maintain next-pointers, face-assignment, and k-values
void insert_vertex_in_edge(Vertex v, Edge e ) {
    // the vertex v is in the middle of edge e
    //                    face
    //                    e1   e2
    // previous-> source  -> v -> target -> next
    //            tw_trg  <- v <- tw_src <- tw_previous
    //                    te2  te1
    //                    twin_face
    
    Edge twin = g[e].twin;
    Vertex src = boost::source( e , g);
    Vertex trg = boost::target( e , g);
    Vertex twin_source = boost::source( twin , g);
    Vertex twin_target = boost::target( twin , g);
    assert( src == twin_target );    
    assert( trg == twin_source );
    
    Face face = g[e].face;
    Face twin_face = g[twin].face;
    Edge previous = previous_edge(e);
    assert( g[previous].face == g[e].face );
    Edge twin_previous = previous_edge(twin);
    assert( g[twin_previous].face == g[twin].face );
    
    Edge e1 = add_edge( src, v ); // these replace e
    Edge e2 = add_edge( v, trg );
    
    // preserve the left/right face link
    g[e1].face = face;
    g[e2].face = face;
    // next-pointers
    g[previous].next = e1;
    g[e1].next = e2;
    g[e2].next = g[e].next;
    
    Edge te1 = add_edge( twin_source, v  ); // these replace twin
    Edge te2 = add_edge( v, twin_target  );
    
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
    
    // update the faces 
    faces[face].edge = e1;
    faces[twin_face].edge = te1;
    
    // finally, remove the old edge
    boost::remove_edge( e   , g);
    boost::remove_edge( twin, g);
}

/// remove given v1-v2 edge
void remove_edge( Vertex v1, Vertex v2) { 
    assert( has_edge(v1,v2) );
    typedef typename std::pair<Edge, bool> EdgeBool;
    EdgeBool result = boost::edge(v1, v2, g );    
    boost::remove_edge( result.first , g );
}

/// remove given v1-v2 edge and its twin
void remove_twin_edges( Vertex v1, Vertex v2) { 
    assert( has_edge(v1,v2) );
    assert( has_edge(v2,v1) );
    typedef typename std::pair<Edge, bool> EdgeBool;
    EdgeBool result1 = boost::edge(v1, v2, g ); 
    EdgeBool result2 = boost::edge(v2, v1, g );    
    boost::remove_edge( result1.first , g );
    boost::remove_edge( result2.first , g );
}

/// remove a degree-two Vertex from the middle of an Edge
// preserve edge-properties (next, face, k)
void remove_deg2_vertex( Vertex v ) {
    //                    face1 e[1]
    //    v1_prev -> v1 -> SPLIT -> v2 -> v2_next
    //    v1_next <- v1 <- SPLIT <- v2 <- v2_prev
    //                  e[0]  face2
    //
    // is replaced with a single edge:
    //                    face1
    //    v1_prev -> v1 ----------> v2 -> v2_next
    //    v1_next <- v1 <---------- v2 <- v2_prev
    //                     face2
    
    EdgeVector v_edges = out_edges(v);
    assert( v_edges.size() == 2);
    assert( source(v_edges[0]) == v && source(v_edges[1]) == v );
     
    Vertex v1 = target( v_edges[0] );
    Vertex v2 = target( v_edges[1] );
    Edge v1_next = g[ v_edges[0] ].next;
    Edge v1_prev = previous_edge( g[ v_edges[0] ].twin );
    Edge v2_next = g[ v_edges[1] ].next;
    Edge v2_prev = previous_edge( g[ v_edges[1] ].twin );
    Face face1 = g[ v_edges[1] ].face;
    Face face2 = g[ v_edges[0] ].face;
    
    Edge new1, new2;
    boost::tie(new1,new2) = add_twin_edges(v1,v2);
    set_next(new1,v2_next);
    set_next(new2,v1_next);
    set_next(v2_prev,new2);
    set_next(v1_prev,new1);
    faces[face1].edge = new1;
    faces[face2].edge = new2;
    g[new1] = g[ v_edges[1] ]; // NOTE: uses EdgeProperties::operator= to copy edge properties
    g[new2] = g[ v_edges[0] ]; //  this sets: params, type, k, face
    remove_twin_edges(v,v1);
    remove_twin_edges(v,v2);
    remove_vertex(v);
}
/// set next-pointer of e1 to e2
void set_next(Edge e1, Edge e2) {
    if (target(e1) != source(e2) ){
        std::cout << " ERROR target(e1) = " << g[target(e1)].index << " source(e2)= " << g[source(e2)].index << "\n"; 
    }
    assert( target(e1) == source(e2) );
    g[e1].next = e2;
}

/// form a face from the edge-list:
/// e1->e2->...->e1
/// for all edges, set edge.face=f, and edge.k=k
void set_next_cycle( std::list<Edge> list, Face f, double k) {
    typename std::list<Edge>::iterator begin,it,nxt,end;
    it= list.begin();
    begin = it;
    faces[f].edge = *it;
    end= list.end();
    for( ; it!=end ; it++ ) {
        nxt = it;
        nxt++;
        if ( nxt != end )
            set_next(*it,*nxt);
        else
            set_next(*it,*begin);
            
        g[*it].face = f;
        g[*it].k = k;
    }
}

/// set next-pointers for the given list (but don't close to form a cycle)
// also set face and k properties for edge
void set_next_chain( std::list<Edge> list, Face f, double k) {
    typename std::list<Edge>::iterator it,nxt,end;
    it= list.begin();
    faces[f].edge = *it;
    //set_next_chain(list);
    //begin = it;    
    end= list.end();
    for( ; it!=end ; it++ ) {
        nxt = it;
        nxt++;
        if ( nxt != end )
            set_next(*it,*nxt);
            
        g[*it].face = f;
        g[*it].k = k;
    }
}

/// set next-pointers for the list
void set_next_chain( std::list<Edge> list ) {
    typename std::list<Edge>::iterator it,nxt,end;
    it= list.begin();
    end= list.end();
    for( ; it!=end ; it++ ) {
        nxt = it;
        nxt++;
        if ( nxt != end )
            set_next(*it,*nxt);
    }
}

/// on a face, search and return the left/right edge from endp
std::pair<Edge,Edge> find_next_prev(Face f, Vertex endp) {
    Edge current = faces[f].edge;
    Edge start_edge = current;
    Edge next_edge = current; // this causes unintialized warning: Edge next_edge = Edge();
    Edge prev_edge = current; // uninitialized warning on ubuntu 11.04 ?
    do {
        Vertex src = source(current);
        Vertex trg = target(current);
        if (src==endp)
            next_edge = current;
        if (trg==endp)
            prev_edge = current;
        current = g[current].next;
    } while (current!=start_edge);
    assert( next_edge != Edge() );
    assert( prev_edge != Edge() );
    //if (debug) {
    //    std::cout << " find_next_prev() next_edge = "; g.print_edge(next_edge); 
    //    std::cout << " find_next_prev() prev_edge = "; g.print_edge(prev_edge);
    //}
    return std::make_pair(next_edge, prev_edge);
}

/// print all faces of graph
void print_faces() {
    for( Face f=0;f<g.num_faces();f++) {
        print_face(f);
    }
}

/// print out vertices on given Face
void print_face(Face f) {
    std::cout << " Face " << f << ": ";
    Edge current = faces[f].edge;
    Edge start=current;
    int num_e=0;
    do {
        Vertex v = source(current);
        std::cout << g[v].index  << "(" << g[v].status  << ")-f"<< g[current].face << "-";
        num_e++;
        assert(num_e<300);
        current = g[current].next;
    } while ( current!=start );
    std::cout << "\n";
}

/// print given edges
void print_edges(EdgeVector& q) {
    BOOST_FOREACH( Edge e, q ) {
        Vertex src = source(e);
        Vertex trg = target(e);
        std::cout << g[src].index << "-" << g[trg].index << "\n";
    }
}

/// print edge
void print_edge(Edge e) {
    Vertex src = source(e);
    Vertex trg = target(e);
    std::cout << g[src].index << "-f" << g[e].face << "-" << g[trg].index << "\n";
}

/// print given vertices
void print_vertices(VertexVector& q) {
    BOOST_FOREACH( Vertex v, q) {
        std::cout << g[v].index << "["<< g[v].type << "]" << " ";
    }
    std::cout << std::endl;
}

}; // end HEDIGraph class definition


} // end hedi namespace
// end halfedgediagram.hpp
