/*  
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  
 *  This file is part of OpenCAMlib.
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

#ifndef FACEGRID_H
#define FACEGRID_H

#include <set>
#include <vector>

#include <boost/multi_array.hpp> 

#include "graph.hpp"
#include "common/point.hpp"

namespace ovd
{

typedef std::vector<FaceProps>                   FacePropVector;
typedef boost::multi_array< FacePropVector* , 2> Grid;
typedef Grid::index                              GridIndex;

typedef unsigned int HEFace;


/// this class is used to quickly find the nearest neighbor among the 
/// existing generators when adding a new generator to a VoronoiDiagram
/// the nearest-neighbor search algorithm is a grid-search
/// a) find the grid-cell that is closest to the the new point
/// b) starting at the closest cell, search for a neighbor
/// c) when one neighbor is found, add points from all cells that are closer or at same distance as found point
/// d) among the found cells, do a brute-force search for the nearest neighbor 
class FaceGrid {
    public:
        FaceGrid();
        /// constructor with parameters: far-limit and number of bins
        FaceGrid(double far, unsigned int n_bins);
        virtual ~FaceGrid();
        /// add face with FaceProps
        void add_face(FaceProps props);
        /// find and return the closest face to a given Point p. Naive algorithm.
        HEFace find_closest_face(const Point& p);
        /// find and return the closest face to a given Point p. Faster grid-search algorithm.
        HEFace grid_find_closest_face(const Point& p);

    private:     
        /// convert an x-coordinate into the corresponding grid index
        GridIndex get_grid_index( double x );
        /// go through all faces in the set and find the closest one to Point p
        HEFace find_closest_in_set( const Point& p );
        /// add neighbors at distance dist from (row,col) to the set
        void insert_faces_from_neighbors(  GridIndex row, GridIndex col , GridIndex dist );
        /// add faces from (row,col) to the set
        void insert_faces_from_bucket(  GridIndex row, GridIndex col );
    // DATA
        /// all generators should fit within this radius from the origin
        double far_radius;
        /// the width of a bin (grid cell)
        double binwidth;
        /// the number of bins in the X and Y direction
        GridIndex nbins;
        /// a 2D grid where each cell contains a 
        /// pointer to a vector of FaceProps, corresponding to the faces in the cell
        Grid* grid;
        typedef std::vector<FaceProps> FaceSetContainer;
        FaceSetContainer face_set;
};

}
#endif
// end facegrid.h
