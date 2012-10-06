#include <string>
#include <iostream>

#include "voronoidiagram.hpp"
#include "version.hpp"

// very simple OpenVoronoi example program
int main(int argc, char **argv) {
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,100); // (r, bins)
    // double r: radius of circle within which all input geometry must fall. use 1 (unit-circle). Scale geometry if necessary.
    // int bins: bins for face-grid search. roughly sqrt(n), where n is the number of sites is good according to Held.
     
    std::cout << "OpenVoronoi version " << ovd::version() << "\n"; // the git revision-string
    std::cout << "build-type " <<ovd::build_type() << "\n"; // the build tybe (Release, Debug, Profile, etc)
    std::cout << "compiler " <<ovd::compiler() << "\n";
    std::cout << "system " <<ovd::system() << "\n";
    std::cout << "processor " <<ovd::processor() << "\n";
    
    ovd::Point p(0,0);
    vd->insert_point_site(p); // this returns an int-handle to the point-site, but we do not use it here.
    
    std::cout << vd->print();

    for(int i=0; i<argc; i++){

        std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
    }
    delete vd;
    return 0;
}
