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

#ifndef VODI_VERTEX_HPP
#define VODI_VERTEX_HPP

#include <map>
#include <cmath>

#include "point.hpp"


namespace ovd {


/// As we incrementally construct the diagram voronoi-vertices can have one of these four different states. 
/// The status is updated as follows:
/// OUT-vertices will not be deleted
/// IN-vertices will be deleted
/// UNDECIDED-vertices have not been examied yet
/// NEW-vertices are constructed on OUT-IN edges
enum VoronoiVertexStatus {OUT, IN, UNDECIDED, NEW };

/// This is the permanent type of a vertex in the diagram. 
/// OUTER vertices are special vertices added in init(), should have degree==4
/// VERTEXGEN are vertex generators, should have degree==0
/// NORMAL are normal voronoi-vertices, should have degree==6  (degree 3 graph with double-edges)
enum VoronoiVertexType {OUTER, NORMAL, VERTEXGEN, ENDPOINT};

/// a map of this type is used by topology-checker to check that all vertices
/// have the expected (correct) degree (i.e. number of edges)
typedef std::map<VoronoiVertexType, unsigned int> VertexDegreeMap;

/// equation-parameters
/// the offset in direction k by a distance t of a general site (point,line,circle) can be expressed as
/// q ( x*x + y*y - t*t ) + a x + b y + c + k t = 0
/// the parameters (q,a,b,k,c) are set as:
/// line:   (0,   a,   b,    k, c          )    line ax+by+c=0  where a*a+b*b=1
/// circle: (1, -2x, -2y, -2kr, x*x+y*y-r*r)    circle center at (x,y) and radius r
/// point:  (1, -2x, -2y,    0, x*x+y*y    )    point at (x,y)
struct Eqp {
    double a;
    double b;
    double c;
    double k;
};

/// Base-class for a voronoi-diagram site, or generator.
class Site {
public:
    Site() {}
    virtual ~Site() {}
    /// return closest point on site to given point p
    virtual Point apex_point(const Point& p) = 0;
    virtual const Point position() const = 0;
    Eqp eqp() {return eq;} 
    bool is_linear() {return isLine(); }
    
    virtual double x() const {assert(0); return 0;}
    virtual double y() const {assert(0); return 0;}
    virtual double r() const {assert(0); return 0;}
    virtual double k() const {assert(0); return 0;}

    virtual double a() const {assert(0); return 0;}
    virtual double b() const {assert(0); return 0;}
    virtual double c() const {assert(0); return 0;}
    
    virtual std::string str() const {assert(0); return "Site";}
    virtual bool isPoint() const { return false;}
    virtual bool isLine() const {  return false;}
    virtual bool in_region(const Point& p) {return false;}
protected:
    Eqp eq;
};

/// point, or vertex site.
class PointSite : public Site {
public:
    PointSite( const Point& p): _p(p) {
        //eq.q = 1;
        eq.a = -2*p.x;
        eq.b = -2*p.y;
        eq.k = 0;
        eq.c = p.x*p.x + p.y*p.y;
    }
    ~PointSite() {}
    virtual Point apex_point(const Point& p) { return _p; }
    virtual const Point position() const {
        return _p;
    }
    virtual double x() const {return _p.x;}
    virtual double y() const {return _p.y;}
    virtual double r() const {return 0;}
    virtual double k() const {return 0;}
    virtual bool isPoint() const {return true;}
    virtual std::string str() const {return "PointSite";}
    virtual bool in_region(const Point& p) {return true;}
private:
    PointSite() {} // don't use!
    Point _p;
};

/// line-segment site
class LineSite : public Site {
public:
    /// create line-site between start and end Point.
    LineSite( const Point& s, const Point& e): _start(s), _end(e) {
        //eq.q = 0;
        eq.a = _end.y - _start.y;
        eq.b = _start.x - _end.x;
        eq.k = 1; // ??
        eq.c = _end.x*_start.y - _start.x*_end.y;
        // now normalize
        double d = sqrt( eq.a*eq.a + eq.b*eq.b );
        eq.a /= d;
        eq.b /= d;
        eq.c /= d;
    }
    ~LineSite() {}
    /// closest point on start-end segment to given point.
    /// project onto line and return either the projected point
    /// or one endpoint of the linesegment
    virtual Point apex_point(const Point& p) {
        Point s_p = p-_start;
        Point s_e = _end - _start;
        double t = s_p.dot(s_e) / s_e.dot(s_e);
        if (t<0)
            return _start;
        if (t>1)
            return _end;
        else {
            return _start + t*(_end-_start);
        }
    }
    virtual const Point position() const {
        assert(0);
        return _start; // FIXME!!
    }
    virtual std::string str() const {return "LineSite";}
    virtual bool in_region(const Point& p) {
        Point s_p = p-_start;
        Point s_e = _end - _start;
        double t = s_p.dot(s_e) / s_e.dot(s_e);
        
        // rounding... UGLY
        double eps = 1e-11;
        if (fabs(t) < eps) 
            t= 0;
        else if ( fabs(t-1.0) < eps )
            t= 1;
        
        //std::cout << "in_region t= " << t << "\n";
        return ( (t>=0) && (t<=1) );
    }
    virtual bool isLine() const {return true;}
    virtual double a() const { return eq.a; }
    virtual double b() const { return eq.b; }
    virtual double c() const { return eq.c; }

private:
    LineSite() {} // don't use!
    Point _start;
    Point _end;
};

/// arc or circle site
class ArcSite : public Site {
public:
    ArcSite( const Point& s, const Point& e, const Point& c, bool dir): _start(s), _end(e), _center(c), _dir(dir) {
        _radius = (_center - _start).norm();
        //eq.q = 1;
        eq.a = -2*_center.x;
        eq.b = -2*_center.y;
        eq.k = -2*_radius; // ?? k
        eq.c = _center.x*_center.x + _center.y*_center.y - _radius*_radius;
    }
    ~ArcSite() {}
    Point apex_point(const Point& p) {
        return Point(0,0); // FIXME
    }
    virtual double x() const {return _center.x;}
    virtual double y() const {return _center.y;}
    virtual double r() const {return _radius;}
    virtual double k() const {return 1;} // ?
    virtual std::string str() const {return "ArcSite";}
private:
    ArcSite() {} // don't use!
    Point _start;
    Point _end;
    Point _center;
    bool _dir; // CW or CCW
    double _radius; // redundant?
};


/// properties of a vertex in the voronoi diagram
/// an object of this type is held in the BGL-graph for each vertex.
class VoronoiVertex {
public:
    VoronoiVertex();
    /// construct vertex at position p with type t
    VoronoiVertex( Point p, VoronoiVertexStatus st);
    /// vertex with given position, status, and type
    VoronoiVertex( Point p, VoronoiVertexStatus st, VoronoiVertexType t);
    void init();
    void reset();
    /// index of vertex
    int index;
    /// vertex status. when the incremental algorithm runs
    /// vertices are marked: undecided, in, out, or new
    VoronoiVertexStatus status;
    VoronoiVertexType type;
    bool in_queue;
    /// the position of the vertex
    Point position;
    
    typedef unsigned int HEFace; 
    HEFace face; // the face associated with this vertex, if type==VERTEXGEN
    friend class VoronoiDiagramChecker;
    
    /// initialize clerance-disk
    void init_dist(const Point& p) {
        r = dist(p);
    }
    /// update clearance-disk
    //void update_dist(const Point& p) {
    //    double d = dist(p);
    //    if (d<r)
    //        r=d;
    //}
    /// return distance to a point from this vertex
    double dist(const Point& p) const { return (position-p).norm(); }
    /// return clearance-disk radius
    double dist() const { return r; }
    /// in-circle predicate 
    double in_circle(const Point& p) const { return dist(p) - r; }
    /// if this vertex is a PointSite, then we store a pointer to the site here.
    Site* site;
    double k3; // the offset-direction to the newly inserted site..
protected:
    /// global vertex count
    static int count;
    /// map for checking topology correctness
    static VertexDegreeMap expected_degree;
    /// clearance-disk radius, i.e. the closest site is at this distance
    double r;
    
};


} // end ocl namespace
#endif
