
// OpenVoronoi polygon example

#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include "voronoidiagram.hpp"
#include "version.hpp"
#include "common/point.hpp"
#include "utility/vd2svg.hpp"

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc,char *argv[]) {
    // Declare the supported options.
    po::options_description desc("This program calculates the voronoi diagram for n random non-intersecting line-segments\nAllowed options");
    desc.add_options()
        ("help", "produce help message")
        //("n", po::value<int>(), "set number of line-segments")
        ("d",  "run in debug-mode")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,10);
    
    if (vm.count("d")) {
        std::cout << "running in debug mode!\n";
        vd->debug_on();
    }
    
    std::cout << "OpenVoronoi version: " << ovd::version() << "\n";
    
    // create a voronoi diagram for a simple polygon
    //
    //   v1 ------ v2
    //   |         |
    //   |         |
    //   v6        v3
    //   |         |
    //   v5 ------v4

    std::vector<int> vertex_ids;
    std::vector<ovd::Point> vertices;
    vertices.push_back( ovd::Point(   0, 0.5) );
    vertices.push_back( ovd::Point( 0.5, 0.5) );
    vertices.push_back( ovd::Point( 0.7, 0.0) );
    vertices.push_back( ovd::Point( 0.4,-0.3) );
    vertices.push_back( ovd::Point(-0.2,-0.2) );
    vertices.push_back( ovd::Point(   0,   0) );
    
    // point-sites must be inserted first.
    // insert_point_site() returns an int-handle that is used when inserting line-segments
    BOOST_FOREACH(ovd::Point p, vertices ) {
        vertex_ids.push_back( vd->insert_point_site(p) );   
    }
    
    // now we insert line-segments
    for (unsigned int n=0;n<vertex_ids.size();n++) {
        int next=n+1;
        if (n==(vertex_ids.size()-1))
            next=0;
        vd->insert_line_site( vertex_ids[n], vertex_ids[next]);
    }
    
    std::cout << " Correctness-check: " << vd->check() << "\n";
    std::cout << vd->print();
    vd2svg("polygon.svg", vd);
    delete vd;
    return 0;
}

