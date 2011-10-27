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
        
        /*
        // for visualizing the closest face (returns it's generator)
        Point getClosestFaceGenerator( const Point p ) {
            HEFace closest_face = fgrid->grid_find_closest_face( p );
            return g[closest_face].generator;
        }
        
        // for visualizing seed-vertex
        Point getSeedVertex( const Point p ) {
            HEFace closest_face = fgrid->grid_find_closest_face( p );
            HEVertex v = find_seed_vertex(closest_face, p);
            return g[ v ].position;
        }
        
        // for visualizing the delete-set
        boost::python::list getDeleteSet( Point p ) { // no const here(?)
            boost::python::list out;
            HEFace closest_face = fgrid->grid_find_closest_face( p );
            HEVertex v_seed = find_seed_vertex(closest_face, p);
            augment_vertex_set(v_seed, p);
            v0.erase( v0.begin() ); // remove the seed vertex, as it is already visualized
            BOOST_FOREACH( HEVertex v, v0) {
                boost::python::list vert;
                vert.append( g[ v ].position );
                vert.append( g[ v ].status );
                out.append( vert );
            }
            reset_status();            
            return out;
        }
        /// for visualizing the delete-edges
        boost::python::list getDeleteEdges( Point p ) {
            boost::python::list out;
            HEFace closest_face = fgrid->grid_find_closest_face( p );
            HEVertex v_seed = find_seed_vertex(closest_face, p);
            augment_vertex_set(v_seed, p);
            EdgeVector del = find_in_in_edges();
            BOOST_FOREACH( HEEdge e, del) {
                boost::python::list edge;
                HEVertex src = g.source(e);
                HEVertex trg = g.target(e);
                edge.append( g[ src ].position );
                edge.append( g[ trg ].position );
                out.append( edge );
            }
            reset_status();            
            return out;
        }
        /// for visualizing the edges to be modified
        boost::python::list getModEdges( Point p ) {
            boost::python::list out;
            HEFace closest_face = fgrid->grid_find_closest_face( p );
            HEVertex v_seed = find_seed_vertex(closest_face, p);
            augment_vertex_set(v_seed, p);
            EdgeVector del = find_in_out_edges();
            BOOST_FOREACH( HEEdge e, del) {
                boost::python::list edge;
                HEVertex src = g.source(e);
                HEVertex trg = g.target(e);
                edge.append( g[ src ].position );
                edge.append( g[ trg ].position );
                out.append( edge );
            }
            reset_status();            
            return out;
        }
        */
        
        /// return list of generators to python
        boost::python::list getGenerators()  {
            boost::python::list plist;
            for ( HEFace f=0;f< g.num_faces();++f ) {
                if ( g[f].site->isPoint() )
                    plist.append( g[f].site->position()  );
            }
            return plist;
        }
        
        /// return list of vd vertices to python
        boost::python::list getVoronoiVertices()  {
            boost::python::list plist;
            BOOST_FOREACH( HEVertex v, g.vertices() ) {
                if ( g[v].type == NORMAL ) {
                    boost::python::list pd;
                    pd.append( g[v].position );
                    pd.append( g[v].dist() );
                    pd.append( g[v].status );
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
                    if ( (g[edge].type == LINE) || (g[edge].type == SEPARATOR) || (g[edge].type == LINESITE) ) {
                        point_list.append( g[v1].position );
                        point_list.append( g[v2].position );
                    } else if ( g[edge].type == PARABOLA ) {
                        double t_src = g[v1].dist();
                        double t_trg = g[v2].dist();
                        double t_min = std::min(t_src,t_trg);
                        double t_max = std::max(t_src,t_trg);
                        int Nmax = 20;
                        double dt = (t_max-t_min)/Nmax;
                        for (int n=0;n<Nmax;n++) {
                            double t = t_min + n*dt;
                            Point pt = g[edge].point(t);
                            point_list.append(pt);
                        }
                    } else {
                        assert(0);
                    } 
                    edge_data.append( point_list );
                    edge_data.append( g[edge].type );
                    edge_list.append(edge_data);
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
