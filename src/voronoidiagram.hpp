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

#include <queue>
#include <set>
#include <boost/tuple/tuple.hpp>

#include "common/point.hpp"
#include "graph.hpp"
#include "vertex_positioner.hpp"
#include "filter.hpp"
#include "kdtree.hpp"

/*! \mainpage OpenVoronoi
 *
 * \author Anders E. Wallin (anders.e.e.wallin "at" gmail.com)
 * \section intro_sec Introduction
 *
 * OpenVoronoi is a c++ library with python bindings (using boost::python) for calculating 2D voronoi-diagrams of point, 
 * line-segment, and circular-arc(not implement yet!) sites.
 * An incremental topology-oriented algorithm is used.
 * 
 * See github for build/install instructions https://github.com/aewallin/openvoronoi
 * 
 * DEB packages for Ubuntu are available at https://launchpad.net/~anders-e-e-wallin/+archive/cam
 * 
 * Output gallery: https://picasaweb.google.com/106188605401091280402/OpenVoronoiExamples
 * 
 * \section Utilities
 * - Offset
 * - Medial-Axis
 * - SVG output
 * 
 * \section Tests
 * Tests are written for CTest.
 * - "ctest -R cpp" runs only the c++ tests (these are fast)
 * - "ctest" runs all tests (some may be slow)
 * Some tests use truetype-tracer (font geometry source), and some use CGAL_RPG (random polygon generator). 
 * 
 * \section coverage Code Coverage Testing
 * - compile using CMAKE_BUILD_TYPE=Coverage  (uses "-fprofile-arcs -ftest-coverage" )
 * - install the library "sudo make install"
 * - Run the custom target coverage-report with "make coverage-report". It will do the following:
 *  - reset lcov counters "lcov --directory ./ --zerocounters"
 *  - run CTest c++ tests with "ctest -R cpptest"
 *  - generate an info-file with "lcov --directory ./ --capture --output-file app.info"
 *  - generate html output with "genhtml --output-directory coverage --title OpenVoronoi Test Coverage app.info"
 * - point your browser to build/doc/index.html to see the output
 * 
 * \section debian Debian source package
 * - See the files in src/deb for more information. 
 * - A debian source package in build/Debian is built with the spackage target, run with "make spackage"
 *  - remember to set the "Release" build-type in CMakeLists.txt
 *  - disable building of tests (these require truetypetracer(ttt) and randompolygon(rpg) which pbuilder/PPA does not find
 * - The source-package can be tested with pbuilder
 *  - To test-build the package (assuming you are on a precise distribution). This will take a long time.
 *   - "sudo pbuilder build openvoronoi_12.02.257-ubuntu1~precise1.dsc" 
 *  - To test-build for other distributions:
 *   - "sudo pbuilder build --distribution lucid openvoronoi_12.02.257-ubuntu1~lucid1.dsc"
 * - The source-package(s) can be uploaded to the Launchpad PPA with dput (this requires that you have write-access to the PPA)
 *  - "dput ppa:anders-e-e-wallin/cam *.changes"
 */
 

namespace ovd
{
/*! 
 * \namespace ovd 
 * \brief OpenVoronoi classes and functions
 */
 
 
class VoronoiDiagramChecker;

/// \brief KD-tree for 2D point location
///
/// a kd-tree is used for nearest-neighbor search when inserting point sites
struct kd_point {
    /// default ctor
    kd_point() : p(0,0), face(0) { }
    /// kd-point with given position and HEFace
    kd_point(Point pt, HEFace f) : p(pt), face(f) { }
    /// kd-point at given position
    kd_point(Point pt) : p(pt), face(0) { }
    /// distance (suared) to given point
    double dist(const kd_point& pt) const {
        return (p.x-pt.p.x)*(p.x-pt.p.x) + (p.y-pt.p.y)*(p.y-pt.p.y); 
    }
    /// return x or y coordinate of Point
    double operator[](unsigned int i) const {
        return i == 0 ? p.x : p.y; 
    }
    /// return x or y coordinate of Point
    double& operator[](unsigned int i) { 
        return i == 0 ? p.x : p.y; 
    }
    Point p; ///< position of 2D PointSite
    HEFace face; ///< the HEFace correspoinding to the PointSite
};

/// type of the KD-tree used for nearest-neighbor search
typedef kdtree::KDTree<kd_point> kd_type; 

/// \brief Voronoi diagram.
///
/// see http://en.wikipedia.org/wiki/Voronoi_diagram
/// 
/// the dual of a voronoi diagram is the delaunay diagram(triangulation).
///  voronoi-faces are dual to delaunay-vertices.
///  vornoi-vertices are dual to delaunay-faces 
///  voronoi-edges are dual to delaunay-edges
class VoronoiDiagram {
public:
    VoronoiDiagram(double far, unsigned int n_bins);
    virtual ~VoronoiDiagram();
    int insert_point_site(const Point& p, int step=0);
    bool insert_line_site(int idx1, int idx2, int step=99); // default step should make algorithm run until the end!
    void insert_arc_site(int idx1, int idx2, const Point& c, bool cw, int step=99);
    
    /// return the far radius
    double get_far_radius() const {return far_radius;}
    /// return number of point sites in diagram
    int num_point_sites() const {return num_psites-3;} // the three initial vertices don't count
    /// return number of line-segments sites in diagram
    int num_line_sites() const {return num_lsites;}
    /// return number of arc-sites in diagram
    int num_arc_sites() const {return num_asites;}
    /// return number of voronoi-vertices
    int num_vertices() const { return g.num_vertices()-num_point_sites(); }
    /// return number of faces in graph
    int num_faces() const { return g.num_faces(); }
    int num_split_vertices() const;
    /// return reference to graph \todo not elegant. only used by vd2svg ?
    HEGraph& get_graph_reference() {return g;}
    
    std::string print() const;
    /// reset vertex index count \todo not very elegant...
    static void reset_vertex_count() { VoronoiVertex::reset_count(); }
    /// turn on debug output
    void debug_on() {debug=true;} 
    /// set silent mode on/off
    void set_silent(bool b) {
        silent=b;
        vpos->set_silent(silent);
    } 
    bool check(); 
    void filter( Filter* flt);
    void filter_reset();
protected:
    /// type for item in VertexQueue. pair of vertex-desxriptor and in_circle predicate
    typedef std::pair<HEVertex, double> VertexDetPair;
    /// \brief comparison-predicate for VertexQueue
    ///
    /// in augment_vertex_set() we grow the delete-tree by processing vertices
    /// one-by-one from a priority_queue. This is the priority_queue sort predicate.
    /// We handle vertices with a large fabs( in_circle() ) first, since we 
    /// believe their predicate to be more reliable.
    class abs_comparison {
    public:
      /// return true if absolute-value of lhs.second is smaller than rhs.second
      bool operator() (const VertexDetPair& lhs, const VertexDetPair&rhs) const {
        return ( fabs(lhs.second) < fabs(rhs.second) );
      }
    };

    /// priority_queue for vertex for processing 
    // sorted by decreasing fabs() of in_circle-predicate, so that the vertices whose IN/OUT status we are 'most certain' about are processed first
    typedef std::priority_queue< VertexDetPair , std::vector<VertexDetPair>, abs_comparison > VertexQueue;
    
    /// \brief data required for adding a new edge
    ///
    /// used in add_edge() for storing information related to
    /// the new edge.
    struct EdgeData {
        HEEdge v1_prv; ///< edge prior to v1
        HEVertex v1;   ///< NEW edge source 
        HEEdge v1_nxt; ///< edge following v1 
        HEEdge v2_prv; ///< edge prior to v2 
        HEVertex v2;   ///< NEW edge target 
        HEEdge v2_nxt; ///< edge following v2 
        HEFace f;      ///< face of v1 and v2 
    };

    void initialize();
    HEVertex   find_seed_vertex(HEFace f, Site* site);
    EdgeVector find_in_out_edges(); 
    EdgeData   find_edge_data(HEFace f, VertexVector startverts, std::pair<HEVertex,HEVertex> segment);
    EdgeVector find_split_edges(HEFace f, Point pt1, Point pt2);
    bool       find_split_vertex(HEFace f, HEVertex& v);
    std::pair<HEVertex,HEVertex> find_endpoints(int idx1, int idx2);
    bool null_vertex_target( HEVertex v , HEVertex& trg);
    void augment_vertex_set( Site* site);        
    bool predicate_c4(HEVertex v);
    bool predicate_c5(HEVertex v);
    void mark_adjacent_faces(HEVertex v, Site* site);
    void mark_adjacent_faces_p( HEVertex v );
    void mark_vertex(HEVertex& v,  Site* site); 
    void   add_vertices( Site* site );
    HEFace add_face(Site* site);
    void   add_edges(HEFace new_f1, HEFace f);        
    void   add_edges(HEFace new_f1, HEFace f, HEFace new_f2, std::pair<HEVertex,HEVertex> seg);
    void   add_edge(EdgeData ed, HEFace new1, HEFace new2=0);
    void   add_separator(HEFace f, HEFace nf, boost::tuple<HEEdge, HEVertex, HEEdge,bool> target, HEVertex endp, Site* s1, Site* s2);
    void   add_split_vertex(HEFace f, Site* s);
    boost::tuple<HEVertex,HEFace,HEVertex,HEVertex,HEFace> find_null_face(HEVertex start, HEVertex other, Point l, Point dir, Site* new_site);
    boost::tuple<HEEdge,HEVertex,HEEdge,bool> find_separator_target(HEFace f, HEVertex endp);
    std::pair<HEVertex,HEFace> process_null_edge(Point dir, HEEdge next_edge , bool k3, bool next_prev);
    HEVertex add_separator_vertex(HEVertex endp, HEEdge edge, Point sep_dir);
    void repair_face( HEFace f );
    void repair_face( HEFace f , std::pair<HEVertex,HEVertex> segs,
                                 std::pair<HEFace,HEFace> nulled_faces,
                                 std::pair<HEFace,HEFace> null_faces );
    void remove_vertex_set();
    void remove_split_vertex(HEFace f);
    void reset_status();
    int num_new_vertices(HEFace f);
// HELPER-CLASSES
    VoronoiDiagramChecker* vd_checker; ///< sanity-checks on the diagram are done by this helper class
    kd_type* kd_tree; ///< kd-tree for nearest neighbor search during point Site insertion
    VertexPositioner* vpos; ///< an algorithm for positioning vertices
// DATA
    typedef std::map<int,HEVertex> VertexMap; ///< type for vertex-index to vertex-descriptor map
    typedef std::pair<int,HEVertex> VertexMapPair; ///< associate vertex index with vertex descriptor
    
    VertexMap vertex_map; ///< map from int handles to vertex-descriptors, used in insert_line_site()
    VertexQueue vertexQueue; ///< queue of vertices to be processed
    HEGraph g; ///< the half-edge diagram of the vd
    double far_radius; ///< sites must fall within a circle with radius far_radius
    int num_psites; ///< the number of point sites
    int num_lsites; ///< the number of line-segment sites
    int num_asites; ///< the number of arc-sites
    FaceVector incident_faces; ///< temporary variable for ::INCIDENT faces, will be reset to ::NONINCIDENT after a site has been inserted
    std::set<HEVertex> modified_vertices; ///< temporary variable for in-vertices, out-vertices that need to be reset after a site has been inserted
    VertexVector v0; ///< IN-vertices, i.e. to-be-deleted
    bool debug; ///< turn debug output on/off
    bool silent; ///< no warnings emitted when silent==true
private:
    VoronoiDiagram(); // don't use.
};

/// \brief error-functor to locate ::SPLIT vertices
///
/// for passing to numerical boost::toms748 root-finding algorithm
class SplitPointError {
public:
    /// \param gi graph
    /// \param split_edge the edge on which we want to position a SPLIT vertex
    /// \param pt1 first point of split-line
    /// \param pt2 second point of split-line
    SplitPointError(HEGraph& gi, HEEdge split_edge, Point pt1, Point pt2) :
    g(gi),  edge(split_edge), p1(pt1), p2(pt2)
    {}
    
    /// \return signed distance to the pt1-pt2 line from edge-point at given offset \a t
    double operator()(const double t) {
        Point p = g[edge].point(t);
        // line: pt1 + u*(pt2-pt1) = p
        //   (p-pt1) dot (pt2-pt1) = u* (pt2-pt1) dot (pt2-pt1)
        
        double u = (p-p1).dot(p2-p1) / ( (p2-p1).dot(p2-p1) );
        Point proj = p1 + u*(p2-p1);
        double dist = (proj-p).norm();
        double sign;
        if ( p.is_right(p1,p2) )
            sign = +1;
        else
            sign = -1;
            
        return sign*dist;
    }
private:
    HEGraph& g;      ///< reference to vd-graph
    HEEdge edge;     ///< the HEEdge on which we position the new SPLIT vertex
    Point p1;     ///< first point of the split-line
    Point p2;    ///< second point of the split-line
};

} // end ovd namespace

// end voronoidiagram.hpp
