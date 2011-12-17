
#include <string>
#include <iostream>

#include <openvoronoi/voronoidiagram.hpp>

#include <boost/random.hpp>

// random points
int main() {
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,100);
    
    std::cout << vd->version() << "\n";
    
    boost::mt19937 rng(42);
    boost::uniform_01<boost::mt19937> rnd(rng);
    
    unsigned int nmax = 100;
    for (unsigned int m=0;m<nmax;m++) {
        double x = rnd();
        double y = rnd();
        ovd::Point p(x,y);
        vd->insert_point_site(p);
    }
    
    std::cout << vd->print();
    return 0;
}

