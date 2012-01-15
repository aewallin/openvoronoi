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

#include <boost/python.hpp>

#include "graph.hpp"
#include "site.hpp"

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
    boost::python::list offset(double t) {
        offset_list = boost::python::list(); // clear the list
        //std::cout << "Offset::offset(t= " << t << ")\n";
        set_flags(t);
        HEFace start;        
        while (find_start_face(start)) {
            offset_walk(start,t);
            //print_status();
        }
        return get_offsets();
    }
    bool find_start_face(HEFace& start) {
        for(HEFace f=0; f<g.num_faces() ; f++) {
            if (face_done[f]==0 ) {
                start=f;
                return true;
            }
        }
        return false;
    }

    void offset_walk(HEFace start,double t) {
        std::cout << " offset_walk() starting on face " << start << "\n";
        HEEdge start_edge =  find_next_offset_edge( g[start].edge , t); // the first edge on the start-face
        boost::python::list loop;
        HEEdge current_edge = start_edge;
        //loop.append( g[current_edge].point(t) );
        boost::python::list pt;
        pt.append( g[current_edge].point(t) );
        pt.append( -1 );
        loop.append(pt);
        do {
            HEEdge next_edge = find_next_offset_edge( g[current_edge].next, t); // the following edge
            //std::cout << "offset-output: "; print_edge(current_edge); std::cout << " to "; print_edge(next_edge); std::cout << "\n";
            HEFace current_face = g[current_edge].face;
            Site* s = g[current_face].site;
            // ask the Site for offset-geometry here.
            
            Ofs* o = s->offset( g[current_edge].point(t), g[next_edge].point(t) );
            //std::cout << o->str();
            // figure out cw or ccw arcs?
            bool cw(true);
            if (!s->isLine() ) // point and arc-sites
                cw = find_cw( o->start(), o->center(), o->end() );
            
            boost::python::list lpt;
                lpt.append( g[next_edge].point(t) );
                lpt.append( o->radius() );
                lpt.append( o->center() );
                lpt.append( cw );
            loop.append(lpt);
            face_done[current_face]=1; // this is WRONG, need to check all offsets done first, not only one (for non-convex cells)
            current_edge = g[next_edge].twin;
        } while (current_edge != start_edge);
        offset_list.append(loop);
    }
    
    // figure out cw or ccw for an arc
    bool find_cw(Point s, Point c, Point e) {
        // arc from current to next edge
        // center at 
        return c.is_right(s,e);
    }
    
    boost::python::list get_offsets() {
        return offset_list;
    }
    // starting at e, find the edge that brackets t
    HEEdge find_next_offset_edge(HEEdge e, double t) {
        // find the first edge that has an offset
        HEEdge start=e;
        HEEdge current=start;
        HEEdge ofs_edge=e;
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
            //g.print_face(f);
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
        //print_status();
    }
    
    bool t_bracket(double a, double b, double t) {
        double min_t = std::min(a,b);
        double max_t = std::max(a,b);
        return ( (min_t<t) && (t<max_t) );
    }
    void print_status() {
        for(HEFace f=0; f<g.num_faces() ; f++) {
            std::cout << (int)face_done[f];
        }
        std::cout << "\n";
    }
private:
    Offset(); // don't use.
    HEGraph& g;
    std::vector<unsigned char> face_done;
    boost::python::list offset_list;
};


} // end namespace

// end file offset.hpp
