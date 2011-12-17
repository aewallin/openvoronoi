
#include <string>
#include <iostream>

#include <openvoronoi/voronoidiagram.hpp>

// very simple OpenVoronoi example program
int main() {
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,100);
    std::cout << vd->version() << "\n";
    ovd::Point p(0,0);
    vd->insert_point_site(p);
    std::cout << vd->print();

    return 0;
}
