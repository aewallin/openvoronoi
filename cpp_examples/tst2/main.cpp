
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#include <openvoronoi/voronoidiagram.hpp>

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

// random points
int main(int argc,char *argv[]) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("n", po::value<int>(), "set number of points")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    
    unsigned int nmax = 100;

    if (vm.count("n")) 
        nmax = vm["n"].as<int>();
    
    std::cout << "Number of points: " << nmax << ".\n";
    int bins = (int)sqrt(nmax);
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,10*bins);
    
    std::cout << "version: " << vd->version() << "\n";
    
    boost::mt19937 rng(42);
    boost::uniform_01<boost::mt19937> rnd(rng);
    
    std::vector<ovd::Point> pts;
    for (unsigned int m=0;m<nmax;m++) {
        double x = 0.8*rnd()-0.5;
        double y = 0.8*rnd()-0.5;
        ovd::Point p(x,y);
        pts.push_back(p);
    }
    boost::timer tmr;
    BOOST_FOREACH(ovd::Point p, pts ) {
        vd->insert_point_site(p);
    }
    //tmr.stop();
    double t = tmr.elapsed();
    std::cout << t << " seconds \n";
    double norm = nmax*log((double)nmax)/log(2.0);
    std::cout << 1e6*t/norm << " us * n*log2(n)\n";
    std::cout << vd->print();
    return 0;
}

