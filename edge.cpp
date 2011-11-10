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

namespace ovd {

EdgeProps::EdgeProps() {
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;
    x[3] = 0;
    x[4] = 0;
    x[5] = 0;
    x[6] = 0;        
    x[7] = 0;
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    y[4] = 0;
    y[5] = 0;
    y[6] = 0;        
    y[7] = 0;
}

/*
* bisector formulas
* x = x1 - x2 - x3*t +/- x4 * sqrt( square(x5+x6*t) - square(x7+x8*t) )
* (same formula for y-coordinate)
* line (line/line)
* parabola (circle/line)
* hyperbola (circle/circle)
* ellipse (circle/circle)
*/
Point EdgeProps::point(double t) const {
    double discr1 = round( sq(x[4]+x[5]*t) - sq(x[6]+x[7]*t) );
    double discr2 = round( sq(y[4]+y[5]*t) - sq(y[6]+y[7]*t) );
    if ( (discr1 >= 0) && (discr2 >= 0) ) {
        double psig = sign ? +1 : -1;
        double nsig = sign ? -1 : +1;
        double xc = x[0] - x[1] - x[2]*t + psig * x[3] * sqrt( discr1 );
        double yc = y[0] - y[1] - y[2]*t + nsig * y[3] * sqrt( discr2 );
        if (xc!=xc) { // test for NaN!
        std::cout << xc << " , " << yc << " t=" << t << "\n";
        print_params();
        assert(0);
        }
        return Point(xc,yc);
    } else {
        std::cout << " warning bisector sqrt(-1) discr1=" << discr1 << " discr2=" << discr2 << "!\n";
        std::cout << " t= " << t << "\n";
        assert(0);
        return Point(0,0);
    }
}

void EdgeProps::set_parameters(Site* s1, Site* s2, bool sig) {
    sign = sig;
    if (s1->isPoint() && s2->isPoint())        // PP
        set_pp_parameters(s1,s2);
    else if (s1->isPoint() && s2->isLine())    // PL
        set_pl_parameters(s1,s2);
    else if (s2->isPoint() && s1->isLine())  {  // LP
        set_pl_parameters(s2,s1);
        sign = !sign;
    } else if (s1->isLine() && s2->isLine())     // LL
        set_ll_parameters(s2,s1);
    else
        assert(0);
        // AP
        // PA
        // AA
        // AL
        // LA
}

void EdgeProps::copy_parameters(EdgeProps& other) {
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
}   

// called for point(s1)-point(s2) edges
void EdgeProps::set_pp_parameters(Site* s1, Site* s2) {
    //std::cout << "set_pp_parameters()\n";
    // x = x1 - x2 - x3*t +/- x4 * sqrt( square(x5+x6*t) - square(x7+x8*t) )
    assert( s1->isPoint() && s2->isPoint() );
    double d = (s1->position() - s2->position()).norm(); //sqrt( sq(xc1-xc2) + sq(yc1-yc2) )
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
    x[7]=0;           // sqrt(  (0+t)^2 - (d/2 +0*t)^2 )
    y[0]=s1->y();     
    y[1]=alfa2*alfa3; 
    y[2]=0; 
    y[3]=-alfa1;       
    y[4]=0;           
    y[5]=+1;          
    y[6]=alfa3;       
    y[7]=0;
}

// called for point(s1)-line(s2) edges
void EdgeProps::set_pl_parameters(Site* s1, Site* s2) {
    assert( s1->isPoint() && s2->isLine() );
    
    type = PARABOLA;
    double alfa3 = s2->a()*s1->x() + s2->b()*s1->y() + s2->c();
    // figure out k, i.e. offset-direction for LineSite
    double kk = 1.0;
    if (alfa3>0.0) {
        kk = -1.0;
    } else {
        sign = !sign;
    }
    
    x[0]=s1->x();       // xc1
    x[1]=s2->a()*alfa3; // alfa1*alfa3
    x[2]=s2->a()*kk;      // -alfa1 = - a2 * k2?
    x[3]=s2->b();       // alfa2 = b2
    x[4]=0;             // alfa4 = r1 
    x[5]=+1;            // lambda1 (allways positive offset from PointSite?)
    x[6]=alfa3;        // alfa3= a2*xc1+b2*yc1+d2?
    x[7]=kk;            // -1 = k2 side of line??

    y[0]=s1->y();       // yc1
    y[1]=s2->b()*alfa3; // alfa2*alfa3
    y[2]=s2->b()*kk;      // -alfa2 = -b2
    y[3]=s2->a();       // alfa1 = a2
    y[4]=0;             // alfa4 = r1
    y[5]=+1;            // lambda1 (allways positive offset from PointSite?)
    y[6]=alfa3;         // alfa3
    y[7]=kk;            // -1 = k2 side of line??
}

// line(s1)-line(s2) edge
void EdgeProps::set_ll_parameters(Site* s1, Site* s2) {  // Held thesis p96
    assert( s1->isLine() && s2->isLine() );
    
    type = LINELINE;
    //sign = true; //sign does not matter because there is no sqrt()?        
    
    double delta =  s1->a()*s2->b() - s1->b()*s2->a() ;
    assert( delta != 0 );
       
    double alfa1 = ( s1->b()*s2->c()-s2->b()*s1->c() ) / delta;
    double alfa2 = ( s2->a()*s1->c()-s1->a()*s2->c() ) / delta;
    double alfa3 = -( s2->b()-s1->b() ) / delta;
    double alfa4 = -( s1->a()-s2->a() ) / delta;
    x[0]= alfa1;  
    x[2]= -alfa3; 
    y[0]= alfa2;         
    y[2]= -alfa4;  

    x[1]=0;
    x[3]=0;  
    x[4]=0;
    x[5]=0;
    x[6]=0;
    x[7]=0;
    y[1]= 0;
    y[3]=0; 
    y[4]=0;
    y[5]=0;
    y[6]=0;
    y[7]=0;
}

double EdgeProps::minimum_t( Site* s1, Site* s2) {
    if (s1->isPoint() && s2->isPoint())        // PP
        return minimum_pp_t(s1,s2);
    else if (s1->isPoint() && s2->isLine())    // PL
        return minimum_pl_t(s1,s2);
    else if (s2->isPoint() && s1->isLine())    // LP
        return minimum_pl_t(s2,s1);
    else if (s1->isLine() && s2->isLine())     // LL
        return 0;
    else
        assert(0);
    
    return -1;
}

double EdgeProps::minimum_pp_t(Site* s1, Site* s2) {
    assert( s1->isPoint() && s2->isPoint() );
    double p1p2 = (s1->position() - s2->position()).norm() ;
    assert( p1p2 >=0 );
    return p1p2/2; // this splits point-point edges at APEX
}

double EdgeProps::minimum_pl_t(Site* s1, Site* s2) {
    double mint = - x[6]/(2.0*x[7]);
    assert( mint >=0 );
    return mint;
}

void EdgeProps::print_params() const {
    std::cout << "x-params: ";
    for (int m=0;m<8;m++)
        std::cout << x[m] << " ";
    std::cout << "\n";
    std::cout << "y-params: ";
    for (int m=0;m<8;m++)
        std::cout << y[m] << " ";
    std::cout << "\n";
}

} // end namespace
