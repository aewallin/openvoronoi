
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
        ("r", " ccw arc ")
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

    std::vector<int> id;
    std::vector<ovd::Point> vertices;
    /*       -----------------
     *      /      arc        \
     *     1                   0
     *     |                   |
     *     l3                  l1
     *     |                   |
     *     3-------l2----------2
     * */
    vertices.push_back( ovd::Point(0.1,0.1) );  // 0
    vertices.push_back( ovd::Point(-0.1,0.1) ); // 1
    vertices.push_back( ovd::Point(0.1,-0.1) ); // 2
    vertices.push_back( ovd::Point(-0.1,-0.1) ); // 3
    
    //vertices.push_back( ovd::Point(-0.03,-0.03) ); // 4?
    
    // point-sites must be inserted first.
    // insert_point_site() returns an int-handle that is used when inserting line-segments
    BOOST_FOREACH(ovd::Point p, vertices ) {
        id.push_back( vd->insert_point_site(p) );   
    }
    if (!vd->check()) return -1;
    
    // now we insert line-segments
    vd->insert_line_site(id[0],id[2]); // l1
    if (!vd->check()) return -1;
    vd->insert_line_site(id[3],id[2]); // l2
    if (!vd->check()) return -1;
    vd->insert_line_site(id[3],id[1]);
    if (!vd->check()) return -1;
    
    // then the arc
    ovd::Point center(0,0);
    if (vm.count("r")) // command line option --r
        vd->insert_arc_site( id[1], id[0] , center, false ); // ccw arc
    else 
        vd->insert_arc_site( id[0], id[1] , center, true ); // cw arc 
    if (!vd->check()) return -1;
        
    std::cout << " Correctness-check: " << vd->check() << "\n";
    std::cout << vd->print();
    std::string filename = "output.svg";
    if (vm.count("r"))
        filename = "output-r.svg";
    vd2svg(filename, vd);
    delete vd;
    return 0;
}

