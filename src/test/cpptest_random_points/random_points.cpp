// OpenVoronoi random points example
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include "voronoidiagram.hpp"
#include "version.hpp"
#include "utility/vd2svg.hpp"

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

/// \test random PointSite test
int main(int argc,char *argv[]) {
    // Declare the supported options.
    po::options_description desc("This program calculates a poisson voronoi diagram\n Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("n", po::value<int>(), "set number of points")
        ("b", po::value<int>(), "set bin-count multiplier")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    
    unsigned int nmax = 100;
    unsigned int binmult = 1;
    if (vm.count("n")) 
        nmax = vm["n"].as<int>();
    if (vm.count("b")) 
        binmult = vm["b"].as<int>();
        
    std::cout << "Number of points: " << nmax << ".\n";
    int bins = std::min( 200, (int)sqrt(nmax) ); // more than 200x200 bins causes a lot of memory use, so limit.
    std::cout << "Using " << binmult*bins << " bins \n";
    
    // try to optimize number of bins by changing binmult and seing what is fast.
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,binmult*bins);
    
    std::cout << "version: " << ovd::version() << "\n";
    
    boost::mt19937 rng(42); // mersenne-twister random number generator
    boost::uniform_01<boost::mt19937> rnd(rng);
    
    std::vector<ovd::Point> pts;
    for (unsigned int m=0;m<nmax;m++) {
        double x = 0.8*rnd()-0.5; // x,y points safely within the unti circle
        double y = 0.8*rnd()-0.5;
        ovd::Point p(x,y);
        pts.push_back(p);
    }
    boost::timer tmr;
    BOOST_FOREACH(ovd::Point p, pts ) {
        vd->insert_point_site(p); // insert each point. This returns an int-handle which we do not use here.
    }
    double t = tmr.elapsed();
    std::cout << t << " seconds \n";
    double norm = nmax*log((double)nmax)/log(2.0);
    std::cout << 1e6*t/norm << " us * n*log2(n)\n";
    std::cout << vd->print();
    vd2svg("random_points.svg", vd);
    delete vd;
    return 0;
}

