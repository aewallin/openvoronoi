
#include <string>
#include <iostream>

#include <openvoronoi/voronoidiagram.hpp>

// very simple OpenVoronoi example program
int main() {
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,100); // (r, bins)
    // double r: radius of circle within which all input geometry must fall. use 1 (unit-circle). Scale geometry if necessary.
    // int bins:  bins for face-grid search. roughly sqrt(n), where n is the number of sites is good according to Held.
     
    std::cout << vd->version() << "\n"; // the git revision-string
    ovd::Point p(0,0);
    vd->insert_point_site(p); // this returns an int-handle to the point-site, but we do not use it here.
    
    std::cout << vd->print();

    return 0;
}
