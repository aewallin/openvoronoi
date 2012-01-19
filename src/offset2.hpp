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

typedef std::vector<HEFace> FaceLoop; // store the t-value too?
typedef std::vector<FaceLoop> OffsetLoops;

// experimental alternative offset approach.
class FaceOffset {
public:
    FaceOffset(HEGraph& gi): g(gi) {
        face_done.clear();
        face_done.assign( g.num_faces(), 1 );
    }
    void offset(double t) {
        offset_list = OffsetLoops(); // clear the list
        set_flags(t); // mark faces as todo or done, based on the t-value, and validity of edges (after filtering).
        HEFace start;        
        while (find_start_face(start)) { // while there are faces that still require offsets
            offset_walk(start,t); // start on the face, and do an offset loop
        }
        //return offset_list;
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
        //std::cout << " offset_walk() starting on face " << start << "\n";
        bool out_in_mode= false; 
        HEEdge start_edge =  find_next_offset_edge( g[start].edge , t, out_in_mode); // the first edge on the start-face
        HEEdge current_edge = start_edge;
        
        FaceLoop loop; // store the output in this loop
        
        // add the first point to the loop.
        //boost::python::list pt;
        //pt.append( g[current_edge].point(t) );
        //pt.append( -1 ); // radius, center, cw
        //loop.append(pt);
        //loop.append
        do {
            out_in_mode = edge_mode(current_edge, t);
            // find the next edge
            HEEdge next_edge = find_next_offset_edge( g[current_edge].next, t, out_in_mode); 
            //std::cout << "offset-output: "; print_edge(current_edge); std::cout << " to "; print_edge(next_edge); std::cout << "\n";
            HEFace current_face = g[current_edge].face;
            { // append the offset-element of current_face to the output
                loop.push_back(current_face);
                /*
                Site* s = g[current_face].site;
                Ofs* o = s->offset( g[current_edge].point(t), g[next_edge].point(t) ); // ask the Site for offset-geometry here.
                bool cw(true);
                if (!s->isLine() ) // point and arc-sites produce arc-offsets, for which cw must be set.
                    cw = find_cw( o->start(), o->center(), o->end() ); // figure out cw or ccw arcs?
                // add offset to output
                boost::python::list lpt;
                    lpt.append( g[next_edge].point(t) );
                    lpt.append( o->radius() );
                    lpt.append( o->center() );
                    lpt.append( cw );
                loop.append(lpt);
                */
            }
            face_done[current_face]=1; // although we may revisit current_face (if it is non-convex), it seems safe to mark it "done" here.
            current_edge = g[next_edge].twin; // NOTE!! twin!! we change from one face to another here.
        } while (current_edge != start_edge);
        offset_list.push_back(loop); // append the created loop to the output
    }
    bool edge_mode(HEEdge e, double t) {
        HEVertex src = g.source(e);
        HEVertex trg = g.target(e);
        double src_r = g[src].dist();
        double trg_r = g[trg].dist();
        if ( (src_r<t) && (t<trg_r) ) {
            return true;
        } else if ((trg_r<t) && (t<src_r) ) {
            return false;
        } else {
            assert(0);
            return false;
        }
    }
    // figure out cw or ccw for an arc
    bool find_cw(Point start, Point center, Point end) {
        // arc from current to next edge
        // center at 
        return center.is_right(start,end); // this only works for arcs smaller than a half-circle !
    }
    
    // starting at e, find the next edge on the face that brackets t
    // we can be in one of two modes.
    // if mode=false then we are looking for an edge where src_t<t<trg_t
    // if mode=true we are looning for an edge where trg_t<t<src_t
    HEEdge find_next_offset_edge(HEEdge e, double t, bool mode) {
        // find the first edge that has an offset
        HEEdge start=e;
        HEEdge current=start;
        HEEdge ofs_edge=e;
        do {
            HEVertex src = g.source(current);
            HEVertex trg = g.target(current);
            double src_r = g[src].dist();
            double trg_r = g[trg].dist();
            if ( !mode && (src_r<t) && (t<trg_r) ) {
                ofs_edge = current;
                break;
            } else if (mode && (trg_r<t) && (t<src_r) ) {
                ofs_edge = current;
                break;
            }
            current =g[current].next;
        } while( current!=start );
        //std::cout << "offset_edge = "; g.print_edge(ofs_edge);
        return ofs_edge;
    }
    
    
    void set_flags(double t) {
        
        // go through all faces and set flag=0 if the face requires an offset.
        for(HEFace f=0; f<g.num_faces() ; f++) {
            HEEdge start = g[f].edge;
            HEEdge current = start;
            do {
                HEVertex src = g.source(current);
                HEVertex trg = g.target(current);
                double src_r = g[src].dist();
                double trg_r = g[trg].dist();
                if (t_bracket(src_r,trg_r,t)) {
                    if ( face_done[f] ) // if 1
                        face_done[f] = 0; // , set to 0. this is a face that requires an offset!
                }
                current = g[current].next;
            } while ( current!=start );
        }
        
        // again go through faces again, and set flag=1 in any edge on the face is invalid
        // this is required because an upstream filter will set valid=false on some edges, but not all, on a face where we do not want offsets.
        for(HEFace f=0; f<g.num_faces() ; f++) {
            HEEdge start = g[f].edge;
            HEEdge current = start;
            do {
                if ( !g[current].valid ) {
                    face_done[f] = 1; // don't offset faces with invalid edges
                }
                current = g[current].next;
            } while ( current!=start );
        }
        

        
        //print_status();
    }
    // is t in (a,b) ?
    bool t_bracket(double a, double b, double t) {
        double min_t = std::min(a,b);
        double max_t = std::max(a,b);
        return ( (min_t<t) && (t<max_t) );
    }
    
    void print() {
        std::cout << "Offset has " << offset_list.size() << " loops.\n";
        int n(0);
        BOOST_FOREACH(FaceLoop l, offset_list ) {
            std::cout << " Loop " << n++ << " has " << l.size() << " faces:\n";
            BOOST_FOREACH(HEFace f, l) {
                std::cout << "   " << f << "\n";
            }
        }
    }
private:
    FaceOffset(); // don't use.
    HEGraph& g;
    // hold a 0/1 flag for each face, indicating if an offset for this face has been produced or not.
    std::vector<unsigned char> face_done;
    OffsetLoops offset_list;
};


} // end namespace

// end file offset2.hpp
