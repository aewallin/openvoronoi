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
#ifndef VODI_PY_H
#define VODI_PY_H


#include "voronoidiagram.hpp"
#include "facegrid.hpp"

namespace ovd
{

/// \brief python wrapper for VoronoiDiagram
class VoronoiDiagram_py : public VoronoiDiagram {
    public:
        VoronoiDiagram_py() : VoronoiDiagram() {}
        /// create diagram with given far-radius and number of bins
        VoronoiDiagram_py(double far, unsigned int n_bins) 
            : VoronoiDiagram( far, n_bins) {}
        
        int insert_point_site1(const Point& p) {
            return insert_point_site(p);
        }
        int insert_point_site2(const Point& p, int step) {
            return insert_point_site(p,step);
        }
        
        bool insert_line_site2(int idx1, int idx2) {
            return insert_line_site( idx1, idx2);
        }
        bool insert_line_site3(int idx1, int idx2, int step) {
            return insert_line_site( idx1, idx2, step);
        }
        
        /// return list of generators to python
        boost::python::list getGenerators()  {
            boost::python::list plist;
            BOOST_FOREACH( HEVertex v, g.vertices() ) {
                if ( g[v].type == POINTSITE || g[v].type == ENDPOINT ) {
                    boost::python::list pd;
                    pd.append( g[v].position );
                    pd.append( g[v].dist() );
                    pd.append( g[v].status );
                    pd.append( g[v].index );
                    plist.append(pd);
                }
            }
            return plist;
        }
        
        /// return list of vd vertices to python
        boost::python::list getVoronoiVertices()  {
            boost::python::list plist;
            BOOST_FOREACH( HEVertex v, g.vertices() ) {
                if ( g[v].type == NORMAL  || g[v].type == APEX || g[v].type == OUTER) {
                    boost::python::list pd;
                    pd.append( g[v].position );
                    pd.append( g[v].dist() );
                    pd.append( g[v].status );
                    pd.append( g[v].index );
                    plist.append(pd);
                }
            }
            return plist;
        }
        /// return list of the three special far-vertices to python
        boost::python::list getFarVoronoiVertices()  {
            boost::python::list plist;
            BOOST_FOREACH( HEVertex v, g.vertices() ) {
                if ( g.degree( v ) == 4 ) {
                    boost::python::list pd;
                    pd.append( g[v].position );
                    pd.append( g[v].dist() );
                    plist.append(pd);
                }
            }
            return plist;
        }
        /// return list of vd-edges to python
        boost::python::list getVoronoiEdges()  {
            boost::python::list edge_list;
            BOOST_FOREACH( HEEdge edge, g.edges() ) { // loop through each edge
                    boost::python::list edge_data;
                    boost::python::list point_list; // the endpoints of each edge
                    HEVertex v1 = g.source( edge );
                    HEVertex v2 = g.target( edge );
                    if ( (g[edge].type == SEPARATOR) || (g[edge].type == LINESITE) || (g[edge].type == OUTEDGE) || (g[edge].type == LINELINE) ) {
                        point_list.append( g[v1].position );
                        point_list.append( g[v2].position );
                    } else if ( g[edge].type == PARABOLA || (g[edge].type == LINE)  ) {
                        double t_src = g[v1].dist();
                        double t_trg = g[v2].dist();
                        double t_min = std::min(t_src,t_trg);
                        double t_max = std::max(t_src,t_trg);
                        //std::cout << g[v1].index << " t=" << g[v1].dist() << " - " << g[v2].index << " t=" << g[v2].dist() << "\n";
                        int Nmax = 20;
                        double dt = (t_max-t_min)/(Nmax-1);
                        for (int n=0;n<Nmax;n++) {
                            double t = t_min + n*dt;
                            Point pt = g[edge].point(t);
                            //std::cout << t << " " << pt << "\n";
                            point_list.append(pt);
                        }
                    } else {
                        assert(0);
                    } 
                    edge_data.append( point_list );
                    edge_data.append( g[edge].type );
                    edge_data.append( g[v1].status ); // source status
                    edge_data.append( g[v2].status ); // target status
                    edge_list.append( edge_data );
            }
            return edge_list;
        }
        /// return edges and generators to python
        boost::python::list getEdgesGenerators()  {
            boost::python::list edge_list;
            BOOST_FOREACH( HEEdge edge, g.edges() ) {
                    boost::python::list point_list; // the endpoints of each edge
                    HEVertex v1 = g.source( edge );
                    HEVertex v2 = g.target( edge );
                    Point src = g[v1].position;
                    Point tar = g[v2].position;
                    int src_idx = g[v1].index;
                    int trg_idx = g[v2].index;
                    point_list.append( src );
                    point_list.append( tar );
                    point_list.append( src_idx );
                    point_list.append( trg_idx );
                    edge_list.append(point_list);
            }
            return edge_list;
        }
    
        // for animation/visualization only, not needed in main algorithm
        EdgeVector find_in_in_edges() { 
            assert( !v0.empty() );
            EdgeVector output; // new vertices generated on these edges
            BOOST_FOREACH( HEVertex v, v0 ) {                                   
                assert( g[v].status == IN ); // all verts in v0 are IN
                BOOST_FOREACH( HEEdge edge, g.out_edges( v ) ) {
                    HEVertex adj_vertex = g.target( edge );
                    if ( g[adj_vertex].status == IN ) 
                        output.push_back(edge); // this is an IN-IN edge
                }
            }
            return output;
        }

        
};


} // end namespace
#endif
// end voronoidiagram_py.h
