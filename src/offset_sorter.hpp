/*  
 *  Copyright 2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

typedef boost::adjacency_list< boost::vecS,            // out-edge storage
                                boost::vecS,            // vertex set storage
                                boost::directedS,       // directed tag
                                OffsetLoop,     // vertex properties
                                boost::no_property,     // edge properties
                                boost::no_property,     // graph-properties 
                                boost::listS            // edge storage
                           > MachiningGraph;
typedef boost::graph_traits<MachiningGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<MachiningGraph>::vertex_iterator   VertexItr;

class OffsetLoopCompare {
public:
    bool operator() (OffsetLoop l1, OffsetLoop l2) {
        return (l1.offset_distance > l2.offset_distance);
    }
};

typedef std::pair<Vertex,OffsetLoop> VertexOffsetLoop;

class VertexOffsetLoopCompare {
public:
    bool operator() (VertexOffsetLoop l1, VertexOffsetLoop l2) {
        return (l1.second.offset_distance < l2.second.offset_distance);
    }
};

template <class Name>
class label_writer {
public:
    label_writer(Name _name) : name(_name) {}
    template <class VertexOrEdge>
    void operator()(std::ostream& out, const VertexOrEdge& v) const {
      out << "[label=\"" << v << " (t="<< name[v].offset_distance <<")\"]";
    }
private:
    Name name;
};
  
class OffsetSorter {
public:
    OffsetSorter() {}
    void add_loop(OffsetLoop l) {
        all_loops.push_back(l);
    }
    void sort_loops() {
        // push loops to a set, so they come out in decreasing offset-distance order, i.e. max offset-distance (innermost) loop first.
        BOOST_FOREACH( const OffsetLoop l, all_loops ) {
            dist_sorted_loops.insert(l);
        }
        // go through the loops and build the machining graph
        
        BOOST_FOREACH( OffsetLoop l, dist_sorted_loops ) {
            std::cout << "loop at " << l.offset_distance << "\n";
            Vertex new_vert = boost::add_vertex(g);
            g[new_vert] = l;
            connect_vertex(new_vert);
        }
        
        write_dotfile();
    }
    void connect_vertex(Vertex v) {
        // try to connect the new vertex to existing vertices
        VertexItr it_begin, it_end, itr;
        boost::tie( it_begin, it_end ) = boost::vertices( g );
        
        std::set< VertexOffsetLoop, VertexOffsetLoopCompare> outside_loops;
        for ( itr=it_begin ; itr != it_end ; ++itr ) {
            if ( v != *itr ) { // don't connect to self
                if ( g[v].offset_distance < g[*itr].offset_distance ){ // connect only if v is contained in *itr, based on offset-distance
                    outside_loops.insert( std::make_pair(*itr, g[*itr]) );
                }
            }
        }
        if (!outside_loops.empty()) {
            VertexOffsetLoop trg_pair = *(outside_loops.begin());
            Vertex trg = trg_pair.first;
            boost::add_edge(v,trg,g);
        }
    }
    void write_dotfile() {
        // write graph to file to see it
        std::filebuf fb;
        fb.open ("test.dot",std::ios::out);
        std::ostream out(&fb);
        label_writer<MachiningGraph> lbl_wrt(g);
        boost::write_graphviz( out, g, lbl_wrt);
    }
protected:
    std::multiset<OffsetLoop, OffsetLoopCompare> dist_sorted_loops; // set of Loops, sorted by decreasing offset-distance
    OffsetLoops all_loops;
    MachiningGraph g;
};


} // end ovd namespace
// end file offset_sorter.hpp
