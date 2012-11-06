// OpenVoronoi random polygon example
// uses random polygon generator "rpg" from https://github.com/aewallin/CGAL_RPG
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include "voronoidiagram.hpp"
#include "checker.hpp"
#include "version.hpp"
#include "common/point.hpp"
#include "utility/vd2svg.hpp"

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <randompolygon/rpg.hpp>

namespace po = boost::program_options;

typedef std::pair<double,double> point;

int main(int argc,char *argv[]) {
    // Declare the supported options.
    po::options_description desc("This program calculates the voronoi diagram for n random polygon\nAllowed options");
    desc.add_options()
        ("help", "produce help message")
        ("n", po::value<int>(), "set number of points")
        ("s", po::value<int>(), "seed for random number generator")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    unsigned int nmax = 100;
    unsigned int seed = 42;
    if (vm.count("n")) 
        nmax = vm["n"].as<int>();
    
    if (vm.count("s")) 
        seed = vm["s"].as<int>();

    std::cout << "Number of vertices in random polygon: " << nmax << "\n";
    int bins = (int)sqrt(nmax);
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,10*bins);
    std::cout << "OpenVoronoi version: " << ovd::version() << "\n";
    
    boost::timer tmr;
    std::vector< point > pt_list = rpg(nmax,seed);
    std::cout << "rpg done in " << tmr.elapsed() << " seconds\n" << std::flush;

    std::vector< int > point_id ;
    
    tmr.restart();
    BOOST_FOREACH(point p, pt_list ) {
        point_id.push_back( vd->insert_point_site( ovd::Point(p.first,p.second)) ); 
    }
    double t_points = tmr.elapsed();
    
    // now we insert line-segments
    tmr.restart();
    for( unsigned int n =0; n<point_id.size();n++) {
        int next = n+1;
        if (n==point_id.size()-1)
            next = 0;
        vd->insert_line_site( point_id[n], point_id[next]); 
    }
    double t_lines = tmr.elapsed();
    
    // sanity-check
    ovd::VoronoiDiagramChecker check(vd->get_graph_reference());
    if (!check.is_valid())
        return -1;
    else
        std::cout << "ovd::VoronoiDiagramChecker OK!\n";
    
    std::cout << "Points: " << t_points << " seconds \n";
    std::cout << "Lines: " << t_lines << " seconds \n";
    double norm = 2*nmax*log(2*(double)nmax)/log(2.0);
    std::cout << "Points: " << 1e6*t_points/norm << " us * n*log2(n)\n";
    std::cout << "Lines: " << 1e6*t_lines/norm << " us * n*log2(n)\n";
    std::cout << vd->print();
    vd2svg("random_polygon.svg", vd);
    delete vd;
    return 0;
}

