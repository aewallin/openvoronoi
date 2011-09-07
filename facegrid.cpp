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

//#define NDEBUG 
#include <cassert>

#include "facegrid.hpp"

namespace ovd {

FaceGrid::FaceGrid() {
    assert(0); // DO NOT use. There are problems because operator= is not defined...
}

FaceGrid::~FaceGrid() {
    for ( GridIndex m=0 ; m<nbins ; ++m ) {
        for ( GridIndex n=0 ; n<nbins ; ++n ) {
            delete (*grid)[m][n];
        }
    }
    delete grid;
}

FaceGrid::FaceGrid(double far, unsigned int n_bins) {
    far_radius = 3.1*far;
    nbins = n_bins;
    binwidth = 2*far_radius/nbins;
    grid = new Grid( boost::extents[n_bins][n_bins] );
    for ( GridIndex m=0 ; m<nbins ; ++m ) {
        for ( GridIndex n=0 ; n<nbins ; ++n ) {
            (*grid)[m][n] = new FacePropVector();
        }
    }
}

// insert f_prop into the correct bucket
void FaceGrid::add_face(FaceProps f_prop) {
    GridIndex row = get_grid_index( f_prop.generator.x );
    GridIndex col = get_grid_index( f_prop.generator.y );
    FacePropVector* bucket = (*grid)[row][col];
    bucket->push_back( f_prop );
} 

GridIndex FaceGrid::get_grid_index( double x ) {
    GridIndex idx;
    idx = (int)( floor( (x+far_radius)/binwidth ) );                
    assert( (idx >= 0) && (idx <= nbins) );     
    return idx;
}

/// simple implementation to find the closest face to the new generator p
HEFace FaceGrid::find_closest_face(const Point& p) {
    HEFace closest_face;
    face_set.clear(); // the set we are searching in
    for ( GridIndex m=0 ; m<nbins ; ++m ) {
        for ( GridIndex n=0 ; n<nbins ; ++n ) {
            insert_faces_from_bucket(  m, n ); // add ALL
        }
    }
    closest_face = find_closest_in_set(  p );
    return closest_face;
}
    
// grid-based search for the closest face to generator p
// grid search algorithm:
// - find the closest grid-cell
// - move in a spiral outward to find the first point
// - when the first point is found we need only search grid cells within a radius = distance-to-found-point
// - find the cells within a radius = distance-to-found-point
// -- in these cells, search all points and find the closest one
HEFace FaceGrid::grid_find_closest_face(const Point& p) {
    face_set.clear();
    GridIndex row = get_grid_index( p.x );
    GridIndex col = get_grid_index( p.y );
    insert_faces_from_bucket( row, col ); // the closest bucket
    GridIndex dist = 0;
    do {
        dist++;
        insert_faces_from_neighbors(  row, col , dist );
    } while (face_set.empty());
    GridIndex max_dist = (int)( ceil( sqrt(2.0)*dist ) ); // expand up to this radius, to be sure to find the closest point
    for (GridIndex d = dist; d<=max_dist;d++)
        insert_faces_from_neighbors(  row, col , d );

    return find_closest_in_set( p ); // among the faces found, find the closest one
}


// go through the HEFace set and return the one closest to p
HEFace FaceGrid::find_closest_in_set(  const Point& p ) {
    HEFace closest_face(0);
    double closest_distance = 30*far_radius; // a big number...
    double d;
    BOOST_FOREACH( FaceProps f, face_set ) {
        d = (f.generator - p).norm_sq();
        if (d<closest_distance ) {
            closest_distance=d;
            closest_face=f.idx;
        }
    }
    return closest_face;
}
    
void FaceGrid::insert_faces_from_neighbors(  GridIndex row, GridIndex col , GridIndex dist ) {
    // insert faces from neighbors of (row,col) at distance dist
    GridIndex min_row = (row >= dist ? row-dist : 0); // N
    GridIndex max_row = (row <= nbins-dist-1 ? row+dist : nbins-1); // S
    GridIndex max_col = (col <= nbins-dist-1 ? col+dist : nbins-1); // E
    GridIndex min_col = (col >= dist  ? col-dist : 0); // W
    for (GridIndex c = min_col; c<=max_col; c++) {
        insert_faces_from_bucket( min_row , c );
        insert_faces_from_bucket( max_row , c );
    }
    for (GridIndex r = min_row; r<=max_row; r++) {
        insert_faces_from_bucket( r, min_col );
        insert_faces_from_bucket( r, max_col );
    }
}
    
void FaceGrid::insert_faces_from_bucket( GridIndex row, GridIndex col ) {
    FacePropVector* bucket = (*grid)[row][col];
    BOOST_FOREACH( FaceProps f, *bucket ) {
        //face_set.insert(f);
        face_set.push_back(f);
    }
}

} // end namespace
// end of facegrid.cpp
