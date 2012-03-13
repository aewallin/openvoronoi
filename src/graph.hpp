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

#pragma once

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
#define VERTEX_CONTAINER boost::listS
#define EDGE_LIST_CONTAINER boost::listS

/// edge-descriptors in the graph
typedef boost::adjacency_list_traits<OUT_EDGE_CONTAINER, 
                                     VERTEX_CONTAINER, 
                                     boost::bidirectionalS, 
                                     EDGE_LIST_CONTAINER >::edge_descriptor HEEdge;
                                     
///face-descriptors in the graph 
// (if there were a traits-class for HEDIGraph we could use it here, instead of "hard coding" the type)
typedef unsigned int HEFace;    

/// Status of faces in the voronoi diagram
enum VoronoiFaceStatus {
    INCIDENT,    /*!< INCIDENT faces contain one or more IN-vertex */
    NONINCIDENT  /*!< NONINCIDENT faces contain only OUT/UNDECIDED-vertices */
    };

/// \brief properties of a face in the voronoi diagram
/// 
struct FaceProps {
    FaceProps() {
        site = 0;
        null = false;
    }
    virtual ~FaceProps() {}
    /// create face with given edge, generator, and type
    FaceProps( HEEdge e , Site* s, VoronoiFaceStatus st) : edge(e), site(s), status(st), null(false) {}
    /// operator for sorting faces
    bool operator<(const FaceProps& f) const {return (this->idx<f.idx);}
    HEFace idx;     ///< face index
    HEEdge edge;     ///< one edge that bounds this face
    Site* site;     ///< the Site for this face 
    VoronoiFaceStatus status;     ///< face status (either ::INCIDENT or ::NONINCIDENT)
    bool null;     ///< flag to indicate null-face

};

/// the type of graph with which we construct the voronoi-diagram
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

typedef boost::graph_traits< HEGraph::BGLGraph >::vertex_descriptor  HEVertex;       ///< vertex descriptor
typedef boost::graph_traits< HEGraph::BGLGraph >::vertex_iterator    HEVertexItr;    ///< vertex iterator
typedef boost::graph_traits< HEGraph::BGLGraph >::edge_iterator      HEEdgeItr;      ///< edge iterator
typedef boost::graph_traits< HEGraph::BGLGraph >::out_edge_iterator  HEOutEdgeItr;   ///< out edge iterator
typedef boost::graph_traits< HEGraph::BGLGraph >::adjacency_iterator HEAdjacencyItr; ///< adj iterator
typedef boost::graph_traits< HEGraph::BGLGraph >::vertices_size_type HEVertexSize;   ///< vertex size

// these containers are used, for simplicity, instead of iterators (like in BGL) when accessing
// adjacent vertices, edges, faces.
// FIXME: it may be faster to rewrite the code so it uses iterators, as does the BGL.
typedef std::vector<HEVertex> VertexVector; ///< vector of vertices 
typedef std::vector<HEFace>   FaceVector;  ///< vector of faces
typedef std::vector<HEEdge>   EdgeVector;  ///< vector of edges

} // end ovd namespace
