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

#include <string>
#include <iostream>

#include "voronoidiagram.hpp"
#include "common/point.hpp"

#include <boost/foreach.hpp>

#include "simple_svg_1.0.0.hpp"

#include <cmath>

#define PI 3.1415926535897932384626433832795
#define CIRCLE_FUZZ 1.e-9

inline ovd::Point scale(ovd::Point p) {
    double s = 500;
    return s*p+s*ovd::Point(1,1);
}

inline double scale(double d) {
    double s = 500;
    return s*d;
}

inline svg::Color get_edge_color(ovd::HEGraph& g, ovd::HEEdge e) {
    if ( g[e].type == ovd::LINESITE )
        return svg::Color::Yellow;
    if ( g[e].type == ovd::PARABOLA )
        return svg::Color::Cyan;
    if ( g[e].type == ovd::SEPARATOR )
        return svg::Color::Magenta;
    if ( g[e].type == ovd::LINELINE )
        return svg::Color::Green;
    if ( g[e].type == ovd::PARA_LINELINE )
        return svg::Color::Lime;
    if ( g[e].type == ovd::OUTEDGE)
         return svg::Color::Orange;
    return svg::Color::Blue;
}

inline void write_line_to_svg(ovd::HEGraph& g, svg::Document& doc, ovd::Point src, ovd::Point trg, svg::Color col) {
    ovd::Point src_p = scale( src );
    ovd::Point trg_p = scale( trg );
    
    svg::Polyline polyline( svg::Stroke(1, col) );
    polyline << svg::Point( src_p.x, src_p.y ) << svg::Point( trg_p.x, trg_p.y );
    doc << polyline;
}

inline void write_arc_to_svg(ovd::HEGraph& g, svg::Document& doc, ovd::Point src, ovd::Point trg, double r, ovd::Point ctr, bool cw, svg::Color col) {
    ovd::Point src_p = scale( src );
    ovd::Point trg_p = scale( trg );
    double radius = scale( r );
    ovd::Point ctr_p = scale( ctr );

    // determine angle theta
    ovd::Point start( src - ctr );
    ovd::Point end( trg - ctr );
    double theta1 = atan2( start.x, start.y );
    double theta2 = atan2( end.x, end.y );
    if ( !cw ) {
        while( (theta2 - theta1) > -CIRCLE_FUZZ )
            theta2 -= 2*PI;
    } else {
        while( (theta2 - theta1) < CIRCLE_FUZZ )
            theta2 += 2*PI;
    }
    double theta = theta2-theta1;

    double x_radius(radius), y_radius(radius), x_axis_rotation(0.);
    bool large_arc_flag( PI <= theta );
    bool sweep_flag( cw );
    svg::EllipticalArc arc(
        svg::Point( src_p.x, src_p.y ),
        x_radius, y_radius, x_axis_rotation,
        large_arc_flag, sweep_flag,
        svg::Point( trg_p.x, trg_p.y ),
        svg::Stroke(1, col)
    );
    doc << arc;
}

inline void write_edge_to_svg(ovd::HEGraph& g, svg::Document& doc, ovd::HEEdge e) {
    ovd::HEVertex src = g.source(e);
    ovd::HEVertex trg = g.target(e);
    ovd::Point src_p = scale( g[src].position );
    ovd::Point trg_p = scale( g[trg].position );
    
    svg::Color col = get_edge_color(g,e);
    svg::Polyline polyline( svg::Stroke(1, col) );
    if ( (g[e].type == ovd::SEPARATOR) || (g[e].type == ovd::LINE) || 
                     (g[e].type == ovd::LINESITE) || (g[e].type == ovd::OUTEDGE) || 
                     (g[e].type == ovd::LINELINE)  || (g[e].type == ovd::PARA_LINELINE)
        ) {
        // edge drawn as two points
        polyline << svg::Point( src_p.x, src_p.y) << svg::Point( trg_p.x, trg_p.y );
    } else if ( g[e].type == ovd::PARABOLA ) { 
        double t_src = g[src].dist();
        double t_trg = g[trg].dist();
        double t_min = std::min(t_src,t_trg);
        double t_max = std::max(t_src,t_trg);
        int nmax=40;
        for (int n=0;n<nmax;n++) {
            double t = t_min + ((t_max-t_min)/((nmax-1)*(nmax-1)))*n*n;
            ovd::Point pt = scale( g[e].point(t) );
            polyline <<  svg::Point(pt.x, pt.y) ;
        }
    }
    doc << polyline;
}

inline void write_pointsite_to_svg(ovd::HEGraph& g, svg::Document& doc, ovd::HEVertex v) {
    if ( g[v].type == ovd::POINTSITE ) {
        ovd::Point p = scale( g[v].position );
        doc << svg::Circle( svg::Point(p.x, p.y), 0.1, 
                              svg::Fill( svg::Color(100, 200, 120)), svg::Stroke(0.01, svg::Color(200, 250, 150) ) 
                            );
    }
}

inline void vd2svg(std::string filename, ovd::VoronoiDiagram* vd) {
    svg::Dimensions dimensions(1024, 1024);
    svg::Document doc(filename, svg::Layout(dimensions, svg::Layout::BottomLeft));
    
    ovd::HEGraph& g = vd->get_graph_reference();
    BOOST_FOREACH( ovd::HEEdge e, g.edges() ) {
        write_edge_to_svg(g,doc,e);
    }
    BOOST_FOREACH( ovd::HEVertex v, g.vertices() ) {
        write_pointsite_to_svg(g,doc,v);
    }
    doc.save();
}
