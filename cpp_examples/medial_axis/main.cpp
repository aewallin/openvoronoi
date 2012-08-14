
#include <string>
#include <iostream>

#include <openvoronoi/medial_axis_filter.hpp>
#include <openvoronoi/medial_axis_walk.hpp>
#include <openvoronoi/voronoidiagram.hpp>
#include <openvoronoi/polygon_interior_filter.hpp>
#include <openvoronoi/utility/vd2svg.hpp>
#include <openvoronoi/version.hpp>

// OpenVoronoi example program. Uses MedialAxis filter to filter the complete Voronoi diagram
// down to the medial axis.
// then uses MedialAxisWalk to walk along the medial axis and draw clearance-disks
int main() {
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,100); // (r, bins)
    // double r: radius of circle within which all input geometry must fall. use 1 (unit-circle). Scale geometry if necessary.
    // int bins:  bins for face-grid search. roughly sqrt(n), where n is the number of sites is good according to Held.
     
    std::cout << ovd::version() << "\n"; // the git revision-string
    ovd::Point p0(-0.1,-0.2);
    ovd::Point p1(0.2,0.1);
    ovd::Point p2(0.4,0.2);
    ovd::Point p3(0.6,0.6);
    ovd::Point p4(-0.6,0.3);

    int id0 = vd->insert_point_site(p0);
    int id1 = vd->insert_point_site(p1);
    int id2 = vd->insert_point_site(p2);
    int id3 = vd->insert_point_site(p3);
    int id4 = vd->insert_point_site(p4);

    vd->insert_line_site(id0, id1);
    vd->insert_line_site(id1, id2);
    vd->insert_line_site(id2, id3);
    vd->insert_line_site(id3, id4);
    vd->insert_line_site(id4, id0);
    vd->check();

    
    // try commenting-out the line below; massive exterior clearance discs will
    // be drawn!
    ovd::polygon_interior_filter pi(true);
    vd->filter(&pi);
    ovd::medial_axis_filter ma;
    vd->filter(&ma);

    // save drawing to svg file.
    svg::Dimensions dimensions(1024, 1024);
    svg::Document doc("medial_axis.svg", svg::Layout(dimensions, svg::Layout::BottomLeft));
    ovd::HEGraph& g = vd->get_graph_reference();
    BOOST_FOREACH( ovd::HEEdge e, g.edges() ) {
        if( g[e].valid ) write_edge_to_svg(g,doc,e);
    }

    // walk the medial axis.
    ovd::MedialAxisWalk maw(g);
    ovd::MedialChainList chain_list = maw.walk();
    std::cout << "walking " << int(chain_list.size()) << " medial axis chains." << std::endl;
    BOOST_FOREACH( ovd::MedialChain chain, chain_list ) { // loop through each chain
        std::cout << "new chain length:" << int(chain.size()) << std::endl;
        BOOST_FOREACH( ovd::MedialPointList pt_list, chain ) { //loop through each medial-point list
            std::cout << "new point list length:" << int(pt_list.size()) << std::endl;
            BOOST_FOREACH( ovd::MedialPoint pt, pt_list ) { //loop through each medial-point
                std::cout << "pt:p:" << pt.p << ", clearance_radius:" << pt.clearance_radius << std::endl;
                if (pt.clearance_radius < 0.001) {
                    std::cout << "(the clearance radius is so small that the rendered circle will be too tiny to see.)" << std::endl;
                }
                ovd::Point ctr( scale( pt.p ) );
                double radius( scale( pt.clearance_radius ) );
                // draw a circle, centered on the medial-axis, with radius = clearance-disc
                svg::Circle clearance_disc( svg::Point( ctr.x, ctr.y ), 2*radius, svg::Fill(), svg::Stroke( 1, svg::Color::Red ) );
                doc << clearance_disc;
            }
        }
    }

    doc.save();
    std::cout << vd->print();
    delete vd;

    return 0;
}
