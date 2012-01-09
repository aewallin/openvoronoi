/* 
 *  Copyright 2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#include <qd/qd_real.h> 

#include "common/point.hpp"

namespace ovd {

/// equation-parameters
/// the offset in direction k by a distance t of a general site (point,line,circle) can be expressed as
/// q ( x*x + y*y - t*t ) + a x + b y + c + k t = 0
/// the parameters (q,a,b,k,c) are set as:
/// line:   (0,   a,   b,    k, c          )    line ax+by+c=0  where a*a+b*b=1
/// circle: (1, -2x, -2y, -2kr, x*x+y*y-r*r)    circle center at (x,y) and radius r
/// point:  (1, -2x, -2y,    0, x*x+y*y    )    point at (x,y)
template<class Scalar>
struct Eq {
    bool q; // true for quadratic, false for linear
    Scalar a;
    Scalar b;
    Scalar c;
    Scalar k;
    
    Eq<Scalar>() {
        a = Scalar(0);
        b = Scalar(0);
        c = Scalar(0);
        k = Scalar(0);
        q = false;
    }
    
    template<class Scalar2>
    Eq<Scalar>& operator=(const Eq<Scalar2>& other) {
        q = other.q;
        a = other.a;
        b = other.b;
        c = other.c;
        k = other.k;
        return *this;
    }
    
    Eq<Scalar>& operator-=(const Eq<Scalar>& other) {
        a-=other.a;
        b-=other.b;
        c-=other.c;
        k-=other.k;
        return *this;
    }
        
    const Eq<Scalar> operator-(const Eq<Scalar>& other) const {
        return Eq<Scalar>(*this) -= other;
    }
    bool operator==(const Eq<Scalar>& other) {
        return ( a==other.a && b==other.b && c==other.c );
    }
    Scalar operator[](int idx) const {
        switch (idx) {
            case 0:
                return a;
            case 1:
                return b;
            case 2:
                return k;
            default:
                assert(0);
                return Scalar(0);
        }
    }
    
};

/// Base-class for a voronoi-diagram site, or generator.
class Site {
public:
    Site() {}
    virtual ~Site() {}
    /// return closest point on site to given point p
    virtual Point apex_point(const Point& p) = 0;
    inline virtual const Point position() const {assert(0); return Point(0,0);}
    virtual const Point start() const {assert(0); return Point(0,0);}
    virtual const Point end() const {assert(0); return Point(0,0);}
    Eq<double> eqp() {return eq;} 
    Eq<double> eqp(double kk) {
        Eq<double> eq2(eq);
        eq2.k *= kk;
        return eq2;
    } 
    Eq<qd_real> eqp_qd(double kk) const {
        Eq<qd_real> eq2;
        eq2=eq;
        eq2.k *= kk;
        return eq2;
    }
    
    bool is_linear() {return isLine(); }
    bool is_quadratic() {return isPoint();}
    
    virtual double x() const {
        std::cout << " WARNING: never call Site !\n";
        assert(0); 
        return 0;
    }
    virtual double y() const {
        std::cout << " WARNING: never call Site !\n";
        assert(0); 
        return 0;
    }
    virtual double r() const {
        std::cout << " WARNING: never call Site !\n";
        assert(0); 
        return 0;
    }
    virtual double k() const {
        std::cout << " WARNING: never call Site !\n";
        assert(0); 
        return 0;
    }

    virtual double a() const {
        std::cout << " WARNING: never call Site !\n";
        assert(0); 
        return 0;
    }
    virtual double b() const {
        std::cout << " WARNING: never call Site !\n";
        assert(0); 
        return 0;
    }
    virtual double c() const {
        std::cout << " WARNING: never call Site !\n";
        assert(0); 
        return 0;
    }
    virtual void set_c(const Point& ) {}
    
    virtual std::string str() const {assert(0); return "Site";}
    virtual std::string str2() const {assert(0); return "Site";}
    inline virtual bool isPoint() const { return false;}
    inline virtual bool isLine() const  { return false;}
    virtual bool in_region(const Point& ) const {
        std::cout << " WARNING: never call Site !\n";
        return false;
    }
    virtual double in_region_t(const Point& ) const {
        std::cout << " WARNING: never call Site !\n";
        return 0;
    } 
    
    typedef unsigned int HEFace;    
    HEFace face;
protected:
    Eq<double> eq;
};

/// point, or vertex site.
class PointSite : public Site {
public:
    PointSite( const Point& p, HEFace f=0): _p(p) {
        face = f;
        eq.q = true;
        eq.a = -2*p.x;
        eq.b = -2*p.y;
        eq.k = 0;
        eq.c = p.x*p.x + p.y*p.y;
    }
    ~PointSite() {}
    virtual Point apex_point(const Point& ) { return _p; }
    inline virtual const Point position() const { return _p; }
    virtual double x() const {return _p.x;}
    virtual double y() const {return _p.y;}
    virtual double r() const {return 0;}
    virtual double k() const {return 0;}
    inline virtual bool isPoint() const {return true;}
    virtual std::string str() const {return "PointSite";}
    virtual std::string str2() const {
        std::string out = "PointSite: ";
        out.append( _p.str() );
        return out;
    }
    virtual bool in_region(const Point& ) const {return true;}
    virtual double in_region_t(const Point& p) const {return -1;} 
private:
    PointSite() {} // don't use!
    Point _p;
};

/// line-segment site
class LineSite : public Site {
public:
    /// create line-site between start and end Point.
    LineSite( const Point& s, const Point& e, double koff, HEFace f = 0): _start(s), _end(e) {
        face = f;
        eq.q = false;
        eq.a = _end.y - _start.y;
        eq.b = _start.x - _end.x;
        eq.k = koff; // ??
        eq.c = _end.x*_start.y - _start.x*_end.y;
        // now normalize
        double d = sqrt( eq.a*eq.a + eq.b*eq.b );
        eq.a /= d;
        eq.b /= d;
        eq.c /= d;
        assert( fabs( eq.a*eq.a + eq.b*eq.b -1.0 ) < 1e-5);
    }
    LineSite( Site& s ) {
        eq = s.eqp();
        face = s.face;
        _start = s.start();
        _end = s.end();
    }
    ~LineSite() {}
    /*
    virtual void set_c(const Point& p) {
        eq.c = -( eq.a * p.x + eq.b * p.y );
    }*/
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
    virtual std::string str() const {return "LineSite";}
    virtual std::string str2() const {
        std::string out = "LineSite: ";
        out.append( _start.str() );
        out.append( " - " );
        out.append( _end.str() );
        return out;
    }
    virtual bool in_region(const Point& p) const{
        double t = in_region_t(p);
        return ( (t>=0) && (t<=1) );
    }
    virtual double in_region_t(const Point& p) const {
        Point s_p = p-_start;
        Point s_e = _end - _start;
        double t = s_p.dot(s_e) / s_e.dot(s_e);
        double eps = 1e-5;
        if (fabs(t) < eps)  // rounding... UGLY
            t = 0.0;
        else if ( fabs(t-1.0) < eps )
            t = 1.0;
        return t;
    }
    inline virtual bool isLine() const {return true;}
    virtual double a() const { return eq.a; }
    virtual double b() const { return eq.b; }
    virtual double c() const { return eq.c; }
    virtual double k() const {
        assert( eq.k==1 || eq.k==-1 );
        return eq.k;
    }
    virtual const Point start() const {return _start;}
    virtual const Point end() const {return _end;}
private:
    LineSite() {} // don't use!
    Point _start;
    Point _end;
};

/// arc or circle site
class ArcSite : public Site {
public:
    ArcSite( const Point& s, const Point& e, const Point& center, bool dir): _start(s), _end(e), _center(center), _dir(dir) {
        _radius = (_center - _start).norm();
        eq.q = true;
        eq.a = -2*_center.x;
        eq.b = -2*_center.y;
        eq.k = -2*_radius; 
        eq.c = _center.x*_center.x + _center.y*_center.y - _radius*_radius;
    }
    ~ArcSite() {}
    Point apex_point(const Point& p) {
        return p+Point(0,0); // FIXME
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


} // end namespace
