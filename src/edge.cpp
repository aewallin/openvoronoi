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

#include "edge.hpp"
#include "common/numeric.hpp"

using namespace ovd::numeric;

namespace ovd {

/// create edge with zero edge-parameters
EdgeProps::EdgeProps() {
    x[0]=0;x[1]=0;x[2]=0;x[3]=0;x[4]=0;x[5]=0;x[6]=0;x[7]=0;
    y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;y[5]=0;y[6]=0;y[7]=0;
    has_null_face = false;
    valid=true;
}

#define stringify( name ) # name
/// table for type_str()
const char* edgeTypeNames[] = {
  stringify( LINE ),
  stringify( LINELINE ),
  stringify( PARA_LINELINE ),
  stringify( OUTEDGE ),
  stringify( PARABOLA ),
  stringify( ELLIPSE ),
  stringify( HYPERBOLA ),
  stringify( SEPARATOR ),
  stringify( NULLEDGE ),
  stringify( LINESITE )
};

/// return type of edge as string
std::string EdgeProps::type_str() const {
    std::ostringstream o;
    o << edgeTypeNames[type];
    return o.str();
}

// the edge is not parametrized by t-value as normal edges
// so we need a projection of sl onto the edge instead
/*
Point EdgeProps::projection_point(Solution& sl) const {
    assert( type == PARA_LINELINE );
    // edge given by
    // x = x[0] + x[1]*t
    // y = y[0] + y[1]*t
    //
    // p = p0 + v*t
    // (p-p0) = v*t
    // t = (p-p0).dot(v) / v.dot(v)
    Point p0(x[0],y[0]);
    Point v(x[1],y[1]);
    std::cout << " edge is from " << p0 << "\n";
    std::cout << " edge direction: " << v << "\n";
    
    double t = (sl.p - p0).dot(v) / v.dot(v);
    // clamp to [0,1]
    if ( t>1)
        t=1;
    else if (t<0)
        t=0;
    std::cout << " projection of solution " << sl.p << " is " << (p0+v*t) << "\n";
    return (p0+v*t);
}*/

/// \brief return point on edge at given offset-distance t
///
/// the eight-parameter formula for a point on the edge is:
/// x = x1 - x2 - x3*t +/- x4 * sqrt( square(x5+x6*t) - square(x7+x8*t) )
Point EdgeProps::point(double t) const {
    double discr1 =  chop( sq(x[4]+x[5]*t) - sq(x[6]+x[7]*t), 1e-14 );
    double discr2 =  chop( sq(y[4]+y[5]*t) - sq(y[6]+y[7]*t), 1e-14 );
    if ( (discr1 >= 0) && (discr2 >= 0) ) {
        double psig = sign ? +1 : -1;
        double nsig = sign ? -1 : +1;
        double xc = x[0] - x[1] - x[2]*t + psig * x[3] * sqrt( discr1 );
        double yc = y[0] - y[1] - y[2]*t + nsig * y[3] * sqrt( discr2 );
        if (xc!=xc) { // test for NaN!
            std::cout << "Edge::point() ERROR: " << xc << " , " << yc << " t=" << t << "\n";
            print_params();
            assert(0);
            return Point(0,0);
        }
        return Point(xc,yc);
    } else {
        std::cout << " warning bisector sqrt(-1) discr1=" << discr1 << " discr2=" << discr2 << "!\n";
        std::cout << " t= " << t << "\n";
        // assert(0);
        return Point(x[0] - x[1] - x[2]*t ,y[0] - y[1] - y[2]*t); // coordinates without sqrt()
    }
}

/// dispatch to setter functions based on type of \a s1 and \a s2
void EdgeProps::set_parameters(Site* s1, Site* s2, bool sig) {
    sign = sig; // sqrt() sign for edge-parametrization
    if (s1->isPoint() && s2->isPoint())        // PP
        set_pp_parameters(s1,s2);
    else if (s1->isPoint() && s2->isLine())    // PL
        set_pl_parameters(s1,s2);
    else if (s2->isPoint() && s1->isLine())  {  // LP
        set_pl_parameters(s2,s1);
        sign = !sign;
    } else if (s1->isLine() && s2->isLine())     // LL
        set_ll_parameters(s2,s1);
    else if (s1->isPoint() && s2->isArc() ) // PA
        set_pa_parameters(s1,s2);
    else if (s2->isPoint() && s1->isArc() ) { // AP
        sign = !sign;
        set_pa_parameters(s2,s1);
        
    } else if (s1->isLine() && s2->isArc() ) // LA
        set_la_parameters(s1,s2);
    else if (s2->isLine() && s1->isArc() ) // AL
        set_la_parameters(s2,s1);
    else
        assert(0);
        // AA
}

/// assignment of edge-parameters
EdgeProps& EdgeProps::operator=(const EdgeProps &other) {
    if (this == &other)
        return *this;
    sign = other.sign;
    x[0] = other.x[0];
    x[1] = other.x[1];
    x[2] = other.x[2];
    x[3] = other.x[3];
    x[4] = other.x[4];
    x[5] = other.x[5];
    x[6] = other.x[6];        
    x[7] = other.x[7];
    y[0] = other.y[0];
    y[1] = other.y[1];
    y[2] = other.y[2];
    y[3] = other.y[3];
    y[4] = other.y[4];
    y[5] = other.y[5];
    y[6] = other.y[6];        
    y[7] = other.y[7];
    face = other.face; 
    null_face = other.null_face;
    has_null_face = other.has_null_face;
    k=other.k; 
    type = other.type;
    valid = other.valid;
    // NOTE we do *not* set: twin, next    
    return *this;
}

/// set edge parameters for PointSite-PointSite edge
void EdgeProps::set_pp_parameters(Site* s1, Site* s2) {
    assert( s1->isPoint() && s2->isPoint() );
    double d = (s1->position() - s2->position()).norm();
    double alfa1 = (s2->x() - s1->x()) / d;
    double alfa2 = (s2->y() - s1->y()) / d;
    double alfa3 = -d/2;
    
    type = LINE;
    x[0]=s1->x();       
    x[1]=alfa1*alfa3; // 
    x[2]=0;  
    x[3]=-alfa2;       
    x[4]=0;             
    x[5]=+1;          
    x[6]=alfa3;       
    x[7]=0;
    y[0]=s1->y();     
    y[1]=alfa2*alfa3; 
    y[2]=0; 
    y[3]=-alfa1;       
    y[4]=0;           
    y[5]=+1;          
    y[6]=alfa3;       
    y[7]=0;
}

/// set ::PARABOLA edge parameters (between PointSite and LineSite).
void EdgeProps::set_pl_parameters(Site* s1, Site* s2) {
    assert( s1->isPoint() && s2->isLine() );
    
    type = PARABOLA;
    double alfa3 = s2->a()*s1->x() + s2->b()*s1->y() + s2->c(); // signed distance to line
    
    // figure out kk, i.e. offset-direction for LineSite
    //double kk = 1.0;
    /*
    if (alfa3>0.0) {
        std::cout << " alfa3>0 ! \n";
        kk = -1.0;
        assert(0); // this branch never taken?
    } else {
    */    
        //std::cout << " alfa3<0 ! \n";
    //    sign = !sign;
    //}
    
    x[0]=s1->x();       // xc1
    x[1]=s2->a()*alfa3; // alfa1*alfa3
    x[2]=s2->a(); //*kk;    // -alfa1 = - a2 * k2?
    x[3]=s2->b();       // alfa2 = b2
    x[4]=0;             // alfa4 = r1 (PointSite has zero radius)
    x[5]=+1;            // lambda1 (allways positive offset from PointSite)
    x[6]=alfa3;         // alfa3= a2*xc1+b2*yc1+d2?
    x[7]=+1; //kk;            // -1 = k2 side of line??

    y[0]=s1->y();       // yc1
    y[1]=s2->b()*alfa3; // alfa2*alfa3
    y[2]=s2->b(); //*kk;    // -alfa2 = -b2
    y[3]=s2->a();       // alfa1 = a2
    y[4]=0;             // alfa4 = r1 (PointSite has zero radius)
    y[5]=+1;            // lambda1 (allways positive offset from PointSite)
    y[6]=alfa3;         // alfa3
    y[7]=+1; //kk;            // -1 = k2 side of line??
}

/// set ::SEPARATOR edge parameters
void EdgeProps::set_sep_parameters(Point& endp, Point& p) {
    type = SEPARATOR;
    double dx = p.x - endp.x;
    double dy = p.y - endp.y;
    double d = (p-endp).norm();
    assert( d > 0 );
    x[0]=endp.x;
    x[2]=-dx/d; // negative of normalized direction from endp to p
    y[0]=endp.y;
    y[2]=-dy/d;
    
    x[1]=0;x[3]=0;x[4]=0;x[5]=0;x[6]=0;x[7]=0;
    y[1]=0;y[3]=0;y[4]=0;y[5]=0;y[6]=0;y[7]=0;
}

/// set edge parametrization for LineSite-LineSite edge (parallel case)
void EdgeProps::set_ll_para_parameters(Site* s1, Site* s2) {
    assert( s1->isLine() && s2->isLine() );
    type = PARA_LINELINE;
    
    // find a point (x1,y1) on the line s1
    // ax+by+c=0
    double x1(0),y1(0);
    if ( fabs( s1->a() ) >  fabs( s1->b() ) ) {
        y1=0;
        x1=-s1->c()/s1->a();
    } else {
        x1=0;
        y1=-s1->c()/s1->b();
    }
    
    // find a point (x2,y2) on the line s2
    // ax+by+c=0
    double x2(0),y2(0);
    if ( fabs( s2->a() ) >  fabs( s2->b() ) ) {
        y2=0;
        x2=-s2->c()/s2->a();
    } else {
        x2=0;
        y2=-s2->c()/s2->b();
    }
    
    // now e.g. the s2 line is given by
    // p = (x2,y2) + t*(-b2, a)
    // and we can find the projection of (x1,y1) onto s2 as
    // p1 = p2 = p0 + t*v
    Point p1(x1,y1);
    Point p2(x2,y2);
    Point v(-s2->b(),s2->a() );
    double t = (p1-p2).dot(v) / v.dot(v);
    Point p1_proj = p2+t*v;

    assert( (p1-p1_proj).norm() > 0 );
    
    // from this point, go a distance d/2 in the direction of the normal
    // to find a point through which the bisector passes
    x1 = x1 + (p1_proj-p1).x / 2;
    y1 = y1 + (p1_proj-p1).y / 2;
    // the tangent of the bisector (as well as the two line-sites) is a vector
    // (-b , a)

    x[0]=  x1;
    x[1]= -s1->b();
    y[0]= y1;
    y[1]= s1->a();
    
    x[2]=0;x[3]=0;x[4]=0;x[5]=0;x[6]=0;x[7]=0;
    y[2]=0;y[3]=0;y[4]=0;y[5]=0;y[6]=0;y[7]=0;
}

/// set edge parametrization for LineSite-LineSite edge
void EdgeProps::set_ll_parameters(Site* s1, Site* s2) {  // Held thesis p96
    assert( s1->isLine() && s2->isLine() );
    type = LINELINE;
    double delta = s1->a()*s2->b() - s1->b()*s2->a() ;

    // (numerically) parallel line segments - the generic LLL solver
    // is numerically unstable for parallel cases
    if (std::abs(delta) <= 1024.0*std::numeric_limits<double>::epsilon())
    {
        set_ll_para_parameters(s1,s2);
        return;
    }
   
    assert( delta != 0 );       
    double alfa1 = ( s1->b()*s2->c()-s2->b()*s1->c() ) / delta;
    double alfa2 = ( s2->a()*s1->c()-s1->a()*s2->c() ) / delta;
    double alfa3 = -( s2->b()-s1->b() ) / delta;
    double alfa4 = -( s1->a()-s2->a() ) / delta;
    
    // point (alfa1,alfa2) is the intersection point between the line-segments
    // vector (-alfa3,-alfa4) is the direction/tangent of the bisector
    x[0]=  alfa1;  
    x[2]= -alfa3; 
    y[0]=  alfa2;         
    y[2]= -alfa4;  

    x[1]=0;x[3]=0;x[4]=0;x[5]=0;x[6]=0;x[7]=0;
    y[1]=0;y[3]=0;y[4]=0;y[5]=0;y[6]=0;y[7]=0;
}

/// set edge parameters when s1 is PointSite and s2 is ArcSite
void EdgeProps::set_pa_parameters(Site* s1, Site* s2) {
    assert( s1->isPoint() && s2->isArc() );
    //std::cout << "set_pa_parameters()\n";
    
    type = HYPERBOLA; // hyperbola or ellipse?
    double lamb2(1.0);
    //if (s2->cw())
    //    lamb2 = +1.0;
    //else
    
    //sign=!sign;
    // distance between centers
    double d = sqrt( (s1->x() - s2->x())*(s1->x() - s2->x()) + (s1->y()-s2->y())*(s1->y()-s2->y()) );
    assert( d > 0 );
    if (d<=s2->r()) {
        lamb2=-1.0;
        sign=!sign;
    }
        
    //if (s2->cw()) {
    //    lamb2 = +1.0;
    //} else {
    //    lamb2 = -1.0;
    //}
    
    double alfa1 = ( s2->x() - s1->x() ) / d;
    double alfa2 = ( s2->y() - s1->y() ) / d;
    double alfa3 = ( s2->r()*s2->r() -  d*d) / (2*d);
    double alfa4 = ( lamb2 * s2->r()  ) / d;
    x[0] = s1->x();
    x[1] = alfa1*alfa3;
    x[2] = alfa1*alfa4;
    x[3] = alfa2;
    x[4] = 0; //r1;  PointSite has zero radius
    x[5] = +1; //lamb1; allways outward offset from PointSite
    x[6] = alfa3;
    x[7] = alfa4;
    
    y[0] = s1->y();
    y[1] = alfa2*alfa3;
    y[2] = alfa2*alfa4;
    y[3] = alfa1;
    y[4] = 0; //r1;     PointSite has zero radius
    y[5] = +1; //lamb1; allways outward offset from PointSite
    y[6] = alfa3;
    y[7] = alfa4;
    //print_params();
}


/// set edge parameters when s1 is ArcSite and s2 is LineSite
void EdgeProps::set_la_parameters(Site* s1, Site* s2) { 
    assert( s1->isLine() && s2->isArc() );
    std::cout << "set_la_parameters() sign= " << sign << " cw= " << s2->cw() << "\n";
    type = PARABOLA;
    double lamb2;
    if (s2->cw())
        lamb2 = +1.0;
    else
        lamb2 = -1.0;
    double alfa1 = s1->a(); //a2
    double alfa2 = s1->b(); //b2
    double alfa3 = ( s1->a()*s2->x() + s1->b()*s2->y() + s1->c() );
    double alfa4 = s2->r();
    double kk = +1; // # positive line-offset
    //if (alfa3 > 0) {
    //    kk = -1;
    //    assert(0);
    //}
    //sign = false;
    // figure out sign?
    
    x[0] = s2->x();
    x[1] = alfa1*alfa3;
    x[2] = alfa1*kk;
    x[3] = alfa2;
    x[4] = alfa4;
    x[5] = lamb2;
    x[6] = alfa3;
    x[7] = kk;
    
    y[0] = s2->y();
    y[1] = alfa2*alfa3;
    y[2] = alfa2*kk;
    y[3] = alfa1;
    y[4] = alfa4;
    y[5] = lamb2;
    y[6] = alfa3;
    y[7] = kk;
    print_params();
}


/// \return minumum t-value for this edge
/// this function dispatches to a helper-function based on the Site:s \a s1 and \a s2
// used only for positioning APEX vertices?
double EdgeProps::minimum_t( Site* s1, Site* s2) {
    if (s1->isPoint() && s2->isPoint())        // PP
        return minimum_pp_t(s1,s2);
    else if (s1->isPoint() && s2->isLine())    // PL
        return minimum_pl_t(s1,s2);
    else if (s2->isPoint() && s1->isLine())    // LP
        return minimum_pl_t(s2,s1);
    else if (s1->isLine() && s2->isLine())     // LL, non-parallel lines allways cross somewhere..
        return 0;
    else if (s1->isPoint() && s2->isArc() ) // PA
        return minimum_pa_t(s1,s2);
    else if (s2->isPoint() && s1->isArc() ) // AP
        return minimum_pa_t(s2,s1);
    else
        assert(0);
    // todo:  AP, AL, AA
    return -1;
}
/// minimum t-value for LINE edge between PointSite and PointSite
double EdgeProps::minimum_pp_t(Site* s1, Site* s2) {
    assert( s1->isPoint() && s2->isPoint() );
    double p1p2 = (s1->position() - s2->position()).norm() ;
    assert( p1p2 >=0 );
    return p1p2/2; // this splits point-point edges at APEX
}
/// minimum t-value for ::PARABOLA edge 
double EdgeProps::minimum_pl_t(Site* , Site* ) {
    double mint = - x[6]/(2.0*x[7]);
    assert( mint >=0 );
    return mint;
}
/// minimum t-value for edge between PointSite and ArcSite
double EdgeProps::minimum_pa_t(Site* s1, Site* s2) {
    assert( s1->isPoint() && s2->isArc() );
    double p1p2 = (s1->position() - s2->apex_point(s1->position()) ).norm(); // - s2->r() ;
    assert( p1p2 >=0 );
    return p1p2/2; // this splits point-point edges at APEX
}
/// print out edge parametrization
void EdgeProps::print_params() const {
    std::cout << "x-params: ";
    for (int m=0;m<8;m++)
        std::cout << x[m] << " ";
    std::cout << "sign= " << sign;
    std::cout << "\n";
    
    std::cout << "y-params: ";
    for (int m=0;m<8;m++)
        std::cout << y[m] << " ";
    std::cout << "\n";
}

} // end namespace
