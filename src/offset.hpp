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

#include "graph.hpp"

namespace ovd
{

/// \brief From a voronoi-diagram, generate offset curve(s).
class Offset {
public:
    Offset(HEGraph& gi): g(gi) {
        face_done.clear();
        face_done.assign( g.num_faces(), 1 );
    }
    void print() {
        std::cout << "Offset: verts: " << g.num_vertices() << "\n";
        std::cout << "Offset: edges: " << g.num_edges() << "\n";
        std::cout << "Offset: faces: " << g.num_faces() << "\n";
    }
    void offset(double t) {
        std::cout << " generating offset at t = " << t << "\n";
        set_flags(t);
        // find a face at which to start offset-loop
        HEFace start= g.HFace();
        for(HEFace f=0; f<g.num_faces() ; f++) {
            if (face_done[f]==0 ) {
                start=f;
                break;
            }
        }
        std::cout << "start offset at " << start << "\n";
        HEEdge start_edge =  find_next_offset_edge( g[start].edge , t); // the first edge on the start-face
        
        HEEdge current_edge = start_edge;
        do {
            HEEdge next_edge = find_next_offset_edge( g[current_edge].next, t); // the following edge
            std::cout << "offset-output: "; print_edge(current_edge); std::cout << " to "; print_edge(next_edge); std::cout << "\n";
            //face_done[start]=1; // this is WRONG, need to check all offsets done first, not only one (for non-convex cells)
            current_edge = g[next_edge].twin;
        } while (current_edge != start_edge);
    }
    
    // starting at e, find the edge that brackets t
    HEEdge find_next_offset_edge(HEEdge e, double t) {
        // find the first edge that has an offset
        HEEdge start=e;
        HEEdge current=start;
        HEEdge ofs_edge=e; // = HEEdge();
        do {
            HEVertex src = g.source(current);
            HEVertex trg = g.target(current);
            double src_r = g[src].dist();
            double trg_r = g[trg].dist();
            if (t_bracket(src_r,trg_r,t)) {
                ofs_edge = current;
                break;
            }
            current =g[current].next;
        } while( current!=start );
        //std::cout << "offset_edge = "; print_edge(ofs_edge);
        return ofs_edge;
    }
    
    
    void set_flags(double t) {
        for(HEFace f=0; f<g.num_faces() ; f++) {
            //print_face(f);
            HEEdge start = g[f].edge;
            HEEdge current = start;
            do {
                HEVertex src = g.source(current);
                HEVertex trg = g.target(current);
                double src_r = g[src].dist();
                double trg_r = g[trg].dist();
                if (t_bracket(src_r,trg_r,t)) {
                    //print_edge(current); // store a potential start-edge for the face here!
                    if ( face_done[f] )
                        face_done[f] = 0; // this is a face that requires an offset!
                }
                current = g[current].next;
            } while ( current!=start );
        }
        for(HEFace f=0; f<g.num_faces() ; f++) {
            std::cout << (int)face_done[f];
        }
        std::cout << "\n";
    }
    
    bool t_bracket(double a, double b, double t) {
        double min_t = std::min(a,b);
        double max_t = std::max(a,b);
        return ( (min_t<t) && (t<max_t) );
    }
    void print_face(HEFace f) { // move to common printing-class?
        std::cout << " Face " << f << ": ";
        HEEdge current = g[f].edge;
        HEEdge start=current;
        int num_e=0;
        do {
            HEVertex v = g.source(current);
            std::cout << g[v].index  << "(" << g[v].status  << ")-f"<< g[current].face << "-";
            num_e++;
            assert(num_e<30);
            current = g[current].next;
        } while ( current!=start );
        std::cout << "\n";
    }
    void print_edge(HEEdge e) {
        HEVertex src = g.source(e);
        HEVertex trg = g.target(e);
        std::cout << g[src].index << "-f" << g[e].face << "-" << g[trg].index << "\n";
    }
private:
    Offset(); // don't use.
    HEGraph& g;
    std::vector<unsigned char> face_done;
};


} // end namespace

// end file offset.hpp
