
#include <string>
#include <iostream>

#include <openvoronoi/voronoidiagram.hpp>
#include <openvoronoi/version.hpp>

// very simple OpenVoronoi example program
int main() {
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,100); // (r, bins)
    // double r: radius of circle within which all input geometry must fall. use 1 (unit-circle). Scale geometry if necessary.
    // int bins:  bins for face-grid search. roughly sqrt(n), where n is the number of sites is good according to Held.
     
    std::cout << "Version: " << ovd::version() << "\n"; // the git revision-string
    std::cout << "Build type: " << ovd::build_type() << "\n";
    std::cout << "Compiler: " << ovd::compiler() << "\n";
    std::cout << "System: " << ovd::system() << "\n";
    std::cout << "Processor: " << ovd::processor() << "\n";
    
    ovd::Point p(0,0);
    vd->insert_point_site(p); // this returns an int-handle to the point-site, but we do not use it here.
    
    std::cout << vd->print();
    delete vd;
    return 0;
}
