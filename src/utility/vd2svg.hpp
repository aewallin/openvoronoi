
// OpenVoronoi polygon example

#include <string>
#include <iostream>

#include <openvoronoi/voronoidiagram.hpp>
#include <openvoronoi/common/point.hpp>

#include <boost/foreach.hpp>

#include "simple_svg_1.0.0.hpp"

ovd::Point scale(ovd::Point p) {
    double s = 500;
    return s*p+s*ovd::Point(1,1);
}

svg::Color get_edge_color(ovd::HEGraph& g, ovd::HEEdge e) {
    if ( g[e].type == ovd::LINESITE )
        return svg::Color::Yellow;
    if ( g[e].type == ovd::PARABOLA )
        return svg::Color::Green;
    return svg::Color::Blue;
}

void write_edge_to_svd(ovd::HEGraph& g, svg::Document& doc, ovd::HEEdge e) {
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

void vd2svg(std::string filename, ovd::VoronoiDiagram* vd) {
    svg::Dimensions dimensions(1024, 1024);
    svg::Document doc(filename, svg::Layout(dimensions, svg::Layout::BottomLeft));
    
    ovd::HEGraph& g = vd->get_graph_reference();
    BOOST_FOREACH( ovd::HEEdge e, g.edges() ) {
        write_edge_to_svd(g,doc,e);
    }
    
    doc.save();
}
