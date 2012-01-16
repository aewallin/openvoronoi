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

#ifndef VODI_G_HPP
#define VODI_G_HPP

#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "common/point.hpp"
#include "common/halfedgediagram.hpp"
#include "vertex.hpp"
#include "site.hpp"
#include "edge.hpp"

// this file contains typedefs used by voronoidiagram.hpp
namespace ovd {

// vecS is slightly faster than listS
// vecS   5.72us * n log(n)
// listS  6.18 * n log(n)
#define OUT_EDGE_CONTAINER boost::listS 
//#define OUT_EDGE_CONTAINER boost::vecS 

// note: cannot use vecS since remove_vertex invalidates iterators/edge_descriptors (?)
#define VERTEX_CONTAINER boost::listS
#define EDGE_LIST_CONTAINER boost::listS

// type of edge-descriptors in the graph
typedef boost::adjacency_list_traits<OUT_EDGE_CONTAINER, 
                                     VERTEX_CONTAINER, 
                                     boost::bidirectionalS, 
                                     EDGE_LIST_CONTAINER >::edge_descriptor HEEdge;
                                     
// type of face-descriptors in the graph 
// (if there were a traits-class for HEDIGraph we could use it here, instead of "hard coding" the type)
typedef unsigned int HEFace;    

/// Status of faces in the voronoi diagram
/// INCIDENT faces contain one or more IN-vertex
/// NONINCIDENT faces contain only OUT-vertices
enum VoronoiFaceStatus {INCIDENT, NONINCIDENT};

/// properties of a face in the voronoi diagram
/// each face stores one edge on the boundary of the face
struct FaceProps {
    FaceProps() {
        site = 0;
    }
    virtual ~FaceProps() {}
    /// create face with given edge, generator, and type
    FaceProps( HEEdge e , Site* s, VoronoiFaceStatus st) : edge(e), site(s), status(st) {}
    /// operator for sorting faces
    bool operator<(const FaceProps& f) const {return (this->idx<f.idx);}
    /// face index
    HEFace idx;
    /// one edge that bounds this face
    HEEdge edge;
    /// the site/generator for this face (either PointSite, LineSite, or ArcSite)
    Site* site;
    /// face status (either incident or nonincident)
    VoronoiFaceStatus status;
};

// the type of graph with which we construct the voronoi-diagram
typedef hedi::half_edge_diagram< OUT_EDGE_CONTAINER,     // out-edges storage
                       VERTEX_CONTAINER,         // vertex set stored here
                       boost::bidirectionalS,    // bidirectional graph.
                       VoronoiVertex,            // vertex properties
                       EdgeProps,                // edge properties
                       FaceProps,                // face properties
                       boost::no_property,       // graph properties
                       EDGE_LIST_CONTAINER       // edge storage
                       > HEGraph;
// NOTE: if these listS etc. arguments ever change, they must also be updated
// above where we do: adjacency_list_traits

typedef boost::graph_traits< HEGraph::BGLGraph >::vertex_descriptor  HEVertex;
typedef boost::graph_traits< HEGraph::BGLGraph >::vertex_iterator    HEVertexItr;
typedef boost::graph_traits< HEGraph::BGLGraph >::edge_iterator      HEEdgeItr;
typedef boost::graph_traits< HEGraph::BGLGraph >::out_edge_iterator  HEOutEdgeItr;
typedef boost::graph_traits< HEGraph::BGLGraph >::adjacency_iterator HEAdjacencyItr;
typedef boost::graph_traits< HEGraph::BGLGraph >::vertices_size_type HEVertexSize;

// these containers are used, for simplicity, instead of iterators (like in BGL) when accessing
// adjacent vertices, edges, faces.
// FIXME: it may be faster to rewrite the code so it uses iterators, as does the BGL.
typedef std::vector<HEVertex> VertexVector;
typedef std::vector<HEFace>   FaceVector;
typedef std::vector<HEEdge>   EdgeVector;  

} // end namespace
#endif


