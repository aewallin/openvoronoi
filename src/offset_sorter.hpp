/*
 *  Copyright (c) 2010-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com).
 *  
 *  This file is part of Openvoronoi 
 *  (see https://github.com/aewallin/openvoronoi).
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 2.1 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <string>
#include <iostream>
#include <fstream> // std::filebuf

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "graph.hpp"
#include "site.hpp"
#include "offset.hpp"

namespace ovd
{

/*struct OffsetVertex {
    Point p;  ///< position (start)
    double r; ///< arc radius (line-vertex is indicated by radius of -1)
    Point c;  ///< arc center
    bool cw;  ///< clockwise (or not)
    HEFace f; ///< corresponding face in the vd-graph
    /// ctor
    OffsetVertex(Point pi, double ri, Point ci, bool cwi, HEFace fi): p(pi), r(ri), c(ci), cw(cwi), f(fi) {}
    /// ctor
    OffsetVertex(Point pi): p(pi), r(-1.), cw(false), f(0) {}
};

/// a single offset loop
struct OffsetLoop {
    std::list<OffsetVertex> vertices;
    double offset_distance;
    void push_back(OffsetVertex v) {
        vertices.push_back(v);
    }
};
*/

//typedef std::list<OffsetLoop> OffsetLoops;

/// a graph for holding pocketing-loops.
/// each vertex in this graph corresponds to an offset loop
/// each edge corresponds a possible link between two loops
typedef boost::adjacency_list< boost::vecS,             // out-edge storage
                                boost::vecS,            // vertex set storage
                                boost::directedS,       // directed tag
                                OffsetLoop,             // vertex properties
                                boost::no_property,     // edge properties
                                boost::no_property,     // graph-properties 
                                boost::listS            // edge storage
                           > MachiningGraph;
/// vertices in MachiningGraph
typedef boost::graph_traits<MachiningGraph>::vertex_descriptor  MGVertex;
typedef boost::graph_traits<MachiningGraph>::edge_descriptor    MGEdge;
/// vertex iterator for MachiningGraph
typedef boost::graph_traits<MachiningGraph>::vertex_iterator   MGVertexItr;

/// sort-predicate for sorting offset loops by offset distance
class OffsetLoopCompare {
public:
    /// sort predicate
    bool operator() (OffsetLoop l1, OffsetLoop l2) {
        return (l1.offset_distance > l2.offset_distance);
    }
};

/// pair offsetloop with the corresponding MachininGraph Vertex
typedef std::pair<MGVertex,OffsetLoop> VertexOffsetLoop;

/// predicate for sorting loops by offset-distance
class VertexOffsetLoopCompare {
public:
    /// sort predicate
    bool operator() (VertexOffsetLoop l1, VertexOffsetLoop l2) {
        return (l1.second.offset_distance < l2.second.offset_distance);
    }
};

/// for writing labels in graphviz-files
template <class Name>
class label_writer {
public:
    label_writer(Name _name) : name(_name) {}     ///< ctor
    /// called by graphviz file-output
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const {
        // vertex-descriptor (an int), and offset-distance)
      out << "[label=\"" << v  << " ("<< name[v].offset_distance <<")\"]";
    }
private:
    Name name; ///< graph that is written to file
};

/// this class sorts offset-loops into a MachiningGraph which should
/// contain the loops in a sensible order for pocket machining
class OffsetSorter {
public:
    OffsetSorter(HEGraph& gi): vdg(gi) {} ///< ctor
    void add_loop(OffsetLoop l) { all_loops.push_back(l); } ///< add an OffsetLoop
    /// sort offset loops
    void sort_loops() {
        // push loops to a set, so they come out in decreasing offset-distance order
        // , i.e. max offset-distance (innermost) loop first.
        BOOST_FOREACH( const OffsetLoop l, all_loops ) {
            distance_sorted_loops.insert(l);
        }

        // go through the loops, in distance order, and add them as vertices to the MachiningGraph 
        BOOST_FOREACH( OffsetLoop l, distance_sorted_loops ) {
            std::cout << "MachiningGraph adding loop at " << l.offset_distance << "\n";
            // each offset loop corresponds to a vertex in the machining-graph
            MGVertex new_vert = boost::add_vertex(g);
            g[new_vert] = l;
            vertex_order.push_back( new_vert ); // output std::vector
        }
        
        // now add edges between vertices
        BOOST_FOREACH( MGVertex v, vertex_order) {
            std::cout << "connecting loop " << v << " at " << g[v].offset_distance << "\n";
            connect_vertex(v); // attempt to connect the new vertex to existing vertices in the graph
        }
        write_dotfile();
    }
    
    /// try to connect the new vertex to existing vertices    
    void connect_vertex(MGVertex v) {

        std::vector< MGVertex > ext_loops; // loops outside the current one
        double current_offset=0;
        bool first = true;
        
        BOOST_FOREACH( MGVertex trg, vertex_order ) {
            if ( (trg != v ) && (g[trg].offset_distance < g[v].offset_distance) ) { 
                // don't connect to self, or to inside loops
                if (first ) {
                    current_offset = g[trg].offset_distance;
                    first = false;
                    std::cout << " first loop outside " << g[v].offset_distance << " is "  << current_offset << "\n";
                }
                // consider only loops just outside of the current one
                if ( g[trg].offset_distance == current_offset ) {
                    if ( inside(trg,v) ) {
                        std::cout << "   connecting " << v << " -> " << trg <<  "\n";
                        ext_loops.push_back( trg );
                    }
                }
            }
        }
        
        // go through the found loops and connect 
        BOOST_FOREACH( MGVertex trg, ext_loops ) {
            boost::add_edge(v,trg,g);
        }
    }
    
    /// return true if the in Vertex is interior to the out Vertex
    bool inside(MGVertex out, MGVertex in) {
        OffsetLoop in_loop = g[in];
        OffsetLoop out_loop = g[out];
        
        std::vector<HEFace> in_loop_faces;
        std::vector<HEFace> out_loop_faces;
        
        //std::cout << " " << in_loop.offset_distance << " in_loop faces: ";
        bool first = true;
        BOOST_FOREACH( OffsetVertex in_ofs_vert, in_loop.vertices ) {
            if (first) {
                first = false;
            } else {
                //std::cout << in_ofs_vert.f << " ";
                in_loop_faces.push_back( in_ofs_vert.f );
            }
        }
        //std::cout << "\n";
        //std::cout << " " << out_loop.offset_distance << " out_loop faces: ";
        first = true;
        BOOST_FOREACH( OffsetVertex out_ofs_vert, out_loop.vertices ) {
            if (first) {
                first = false;
            } else {
                //std::cout << out_ofs_vert.f << " ";
                out_loop_faces.push_back( out_ofs_vert.f );
            }
        }
        //std::cout << "\n";
        
        std::set<HEVertex> in_enclosed = loop_enclosed_vertices(in_loop_faces);
        std::cout << "  IN " << in << " enclosed vertices: ";
        BOOST_FOREACH(HEVertex tv, in_enclosed) {
            std::cout << vdg[tv].index << " ";
        }
        std::cout << "\n";
        
        std::set<HEVertex> out_enclosed = loop_enclosed_vertices(out_loop_faces);
        std::cout << "  OUT " << out << " enclosed vertices: ";
        
        BOOST_FOREACH(HEVertex tv, out_enclosed) {
            std::cout << vdg[tv].index << " ";
        }
        std::cout << "\n";
        
        BOOST_FOREACH( HEVertex inv, in_enclosed) {
            if ( out_enclosed.find( inv ) != out_enclosed.end() )
                return true;
        }
        /*
        BOOST_FOREACH( OffsetVertex in_ofs_vert, in_loop.vertices ) {
            BOOST_FOREACH( OffsetVertex out_ofs_vert, out_loop.vertices ) {
                if (in_ofs_vert.f!=0 && out_ofs_vert.f!=0 && in_ofs_vert.f == out_ofs_vert.f) {
                    std::cout << "  Match!\n";
                    return true;
                }
            }
        } 
        std::cout << "  NO Match!\n";
        */
        return false;
    }

/// find the vd-vertices enclosed by the current offset loop
std::set<HEVertex> loop_enclosed_vertices( std::vector<HEFace> in_loop_faces) {
        std::vector< VertexVector > in_loop_vertices;
        // now go through each in_loop face, and collect all interior vd-vertices
        BOOST_FOREACH(HEFace f, in_loop_faces) {
            VertexVector in_vertices = vdg.face_vertices(f);
            in_loop_vertices.push_back(in_vertices);
            //std::cout << " face " << f << " has " << in_vertices.size() << " vertices\n";
        }
        
        std::set< HEVertex > in_loop_enclosed_vertices;
        // now figure out which of the vertices are enclosed
        for(unsigned int i=0;i<in_loop_vertices.size();++i) {
            // i chooses the face we work on
            // now choose a vertex
            VertexVector test_face_verts = in_loop_vertices[i];
            for (unsigned int iv=0;iv<in_loop_vertices[i].size();++iv) { 
                HEVertex test_vert = test_face_verts[iv];
                for(unsigned int j=0;j<in_loop_vertices.size();++j) { // choose first comp-face
                    if ( j!=i ) {
                        VertexVector comp1_face_verts = in_loop_vertices[j];
                        for(unsigned int k=0;k<in_loop_vertices.size();++k) { // second comp-face
                            if ( j!=k && i!=k ) {
                                VertexVector comp2_face_verts = in_loop_vertices[k];
                                // if we have three faces that each contains a Vertex then it is interior.
                                bool comp1_ok=false;
                                bool comp2_ok=false;
                                BOOST_FOREACH(HEVertex comp1_vert, comp1_face_verts) {
                                    if (test_vert == comp1_vert)
                                        comp1_ok = true;
                                }
                                BOOST_FOREACH(HEVertex comp2_vert, comp2_face_verts) {
                                    if (test_vert == comp2_vert)
                                        comp2_ok = true;
                                }
                                if (comp1_ok && comp2_ok)
                                    in_loop_enclosed_vertices.insert(test_vert);
                            }
                        }   
                    }
                }
            }
        }
    //assert( !in_loop_enclosed_vertices.empty() );
    return in_loop_enclosed_vertices;
}    
    /// write the MachiningGraph to a .dot file for visualization
    void write_dotfile() {        
        std::filebuf fb;
        fb.open ("test.dot",std::ios::out);
        std::ostream out(&fb);
        label_writer<MachiningGraph> lbl_wrt(g);
        boost::write_graphviz( out, g, lbl_wrt);
    }
       
protected:
    /// set of Loops, sorted by decreasing offset-distance
    std::multiset<OffsetLoop, OffsetLoopCompare> distance_sorted_loops; 
    std::vector<MGVertex> vertex_order; ///< the output of this algorithm, vertices in sorted order
    OffsetLoops all_loops; ///< all loops we deal with
    MachiningGraph g; ///< machining-graph constructed when this algorithm runs
    HEGraph& vdg; ///< vd-graph
};


} // end ovd namespace
// end file offset_sorter.hpp
