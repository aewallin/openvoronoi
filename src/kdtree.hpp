/*  
 *  Copyright 2012 Anders E. Wallin (anders.e.e.wallin "at" gmail.com)
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

namespace kdtree {

/// a node in the KD-tree
template<class point_type>
struct kd_node {
    /// create node
    /// \param p position of node
    /// \param d direction of cut/split at this node
    /// \param l left child kd_node
    /// \param r right child kd_node
    kd_node( point_type p, int d, kd_node* l, kd_node* r) :
        pos(p), dir(d), left(l), right(r) {}
    virtual ~kd_node() {
        if (left)
            delete left;
        if (right)
            delete right;
    }
    point_type pos; ///< position of the node
    int dir; ///< direction of cut at this node
    kd_node *left;     ///< pointer to left child node 
    kd_node *right;    ///< pointer to right child node 
};

/// hyper rectangle
template<class point_type>
struct  kd_hyperrect {
    /// ctor
    kd_hyperrect(int dimi, const point_type mini, const point_type maxi) : 
        dim(dimi) {
        min = mini;
        max = maxi;
    }
    /// extend hyperrectangle to include the given point
    void extend(const point_type& pos ) {
        for (int i=0; i < dim; i++) {
            if (pos[i] < min[i]) {
                min[i] = pos[i];
            }
            if (pos[i] > max[i]) {
                max[i] = pos[i];
            }
        }
    }
    
    int dim; ///< number of dimensions
    point_type min; ///< minimum point (e.g. lower left corner)
    point_type max; ///< maximum point (e.g. upper right corner)
};

/// \brief KD-tree for nearest neigbor search
template<class point_type>
class KDTree {
public:
    /// ctor
    KDTree(int dim = 3) : dim_(dim), root_(0), rect_(0) {
    }
    virtual ~KDTree() {
        if (rect_)
            delete rect_;
        if (root_)
            delete root_;
    }
    /// insert given point into tree
    int insert( const point_type pos) {
        return insert(root_, pos);
    }
    /// return a point in the tree that is nearest to the given point
    /// false is returned if the tree is empty. 
    std::pair<point_type,bool> nearest( const point_type& pos) {
        if (!root_) return std::make_pair(pos,false);
        num_nearest_i_calls=0;
        kd_hyperrect<point_type> rect = *rect_;
        kd_node<point_type>* result;
        /* Our first guesstimate is the root node */
        result = root_;
        /* Search for the nearest neighbour recursively */
        double result_dist = root_->pos.dist(pos);
        kd_nearest_i( root_, pos, result, result_dist , rect);
        return std::make_pair( result->pos, true);
    }
    /// for debug, return the number of function calls made during a search
    int get_num_calls() {return num_nearest_i_calls;}
    /// print output of tree
    void print_tree() {
        if (root_)
            print_node(root_,0);
    }
private:
    /// insert given point into tree
    int insert(kd_node<point_type>* node, const point_type pos ) {
        if ( insert_rec( node, pos,  0) ) {
            return -1; // error (?)
        }

        if ( rect_ == 0 ) 
            rect_ = new kd_hyperrect<point_type>(dim_,pos,pos); 
        else
            rect_->extend( pos );

        return 0;
    }
    
    /// recursive point-insertion function
    int insert_rec( kd_node<point_type>*& node, const point_type pos, int dir) {
        if (node == 0) {
            node = new kd_node<point_type>(pos,dir,0,0);
            if (root_==0) {
                root_=node;
            }
            return 0;
        } else {
        
            // else an existing node
            //node = *nptr;
            int new_dir = (node->dir + 1) % dim_;
            if( pos[node->dir] < node->pos[ node->dir] )  {
                //std::cout << " left "; print_node(node,0);
                return insert_rec( node->left, pos,  new_dir  );
            } else {
                //std::cout << " right "; print_node(node,0);
                return insert_rec( node->right, pos, new_dir );
            }
        }
    }
    /// square
    double sq(double x) {return x*x;}
    

    /// find nearest point
    void kd_nearest_i(kd_node<point_type> *node, const point_type pos, 
                       kd_node<point_type>*& result, double& result_dist_sq , kd_hyperrect<point_type>& rect) {
        int dir = node->dir;
        num_nearest_i_calls++;
        //int i;
        //double dummy;
        kd_node<point_type> *nearer_subtree, *farther_subtree;
        double nearer_hyperrect_coord, farther_hyperrect_coord;

        /* Decide whether to go left or right in the tree */
        //dummy = pos[dir] - node->pos[dir];
        //if (dummy <= 0) {
            bool b;
        if (pos[dir] <= node->pos[dir]) {
            nearer_subtree = node->left;
            farther_subtree = node->right;
            nearer_hyperrect_coord  = rect.max[dir];
            farther_hyperrect_coord = rect.min[dir];
            b= true;
        } else {
            nearer_subtree = node->right;
            farther_subtree = node->left;
            nearer_hyperrect_coord  = rect.min[dir];
            farther_hyperrect_coord = rect.max[dir];
            b=false;
        }

        if (nearer_subtree) {
            //std::cout << " recurse into nearer subtree \n";
            /* Slice the hyperrect to get the hyperrect of the nearer subtree */
            double dummy = nearer_hyperrect_coord;
            //nearer_hyperrect_coord = node->pos[dir];
            if (b)
                rect.max[dir] = node->pos[dir];
            else
                rect.min[dir] = node->pos[dir];
            /* Recurse down into nearer subtree */
            kd_nearest_i(nearer_subtree, pos, result, result_dist_sq, rect);
            /* Undo the slice */
            if (b)
                rect.max[dir] = dummy; //node->pos[dir];
            else
                rect.min[dir] = dummy; //node->pos[dir];
            //rect->min[dir] = //nearer_hyperrect_coord = dummy;
        }

        /* Check the distance of the point at the current node, compare it
         * with our best so far */
        double dist_sq = node->pos.dist( pos );
        //for(int i=0; i < dim_; i++) {
        //    dist_sq += sq(node->pos[i] - pos[i]);
        //}
        if (dist_sq < result_dist_sq) {
            result = node;
            result_dist_sq = dist_sq;
        }

        if (farther_subtree) {
            //std::cout << " recurse into farther subtree \n";
            /* Get the hyperrect of the farther subtree */
            double dummy = farther_hyperrect_coord;
            if (b)
                rect.min[dir] = node->pos[dir];
            else
                rect.max[dir] = node->pos[dir];
            //farther_hyperrect_coord = node->pos[dir]; // this changes rect_ !!
            /* Check if we have to recurse down by calculating the closest
             * point of the hyperrect and see if it's closer than our
             * minimum distance in result_dist_sq. */
            if (hyperrect_dist_sq(rect, pos) < result_dist_sq) 
            {
                /* Recurse down into farther subtree */
                kd_nearest_i(farther_subtree, pos, result, result_dist_sq, rect);
                      //double dist_sq = node->pos.dist( pos );
            //for(int i=0; i < dim_; i++) {
            //    dist_sq += sq(node->pos[i] - pos[i]);
            //}
            /*
            if (dist_sq < result_dist_sq) {
                result = node;
                result_dist_sq = dist_sq;
            }*/
            }
            /* Undo the slice on the hyperrect */
            if (b)
                rect.min[dir] = dummy;
            else
                rect.max[dir] = dummy;
                
            //*farther_hyperrect_coord = dummy;
        }
    }
    /// return distance squared from hyperrect to given point
    double hyperrect_dist_sq(kd_hyperrect<point_type>& rect, const point_type& pos){
        double result = 0;
        for (int i=0; i < dim_; i++) {
            if (pos[i] < rect.min[i]) {
                result += sq(rect.min[i] - pos[i]);
            } else if (pos[i] > rect.max[i]) {
                result += sq(rect.max[i] - pos[i]);
            }
        }
        return result;
    }

    /// text output
    void print_node(kd_node<point_type>* n, int d) {
        for(int i=0;i<d;++i)
            std::cout << " ";
        std::cout << "d=" << n->dir << " node at ";
        for (int i=0;i<dim_;++i)
            std::cout << n->pos[i] << " "; // << n->pos[1] << "  " << " " << n->pos[2] << "\n";
        std::cout << "\n";
        if (n->left)
            print_node(n->left,d+1);
        if (n->right)
            print_node(n->right,d+1);
    }

    int num_nearest_i_calls; ///< for debug, keep track of function calls
    int dim_; ///< number of dimensions 
    kd_node<point_type>* root_; ///< root of the tree
    kd_hyperrect<point_type>* rect_; ///< hyperrectangle
};

} // kdtree namespace
