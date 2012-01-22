
#include <string>
#include <iostream>

#include <openvoronoi/offset.hpp>
#include <openvoronoi/voronoidiagram.hpp>

// very simple OpenVoronoi example program
int main() {
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,100); // (r, bins)
    // double r: radius of circle within which all input geometry must fall. use 1 (unit-circle). Scale geometry if necessary.
    // int bins:  bins for face-grid search. roughly sqrt(n), where n is the number of sites is good according to Held.
     
    std::cout << vd->version() << "\n"; // the git revision-string
    ovd::Point p0(-0.1,-0.2);
    ovd::Point p1(0.2,0.1);
    ovd::Point p2(0.4,0.2);
    ovd::Point p3(0.6,0.6);
    ovd::Point p4(-0.6,0.3);

    int id0 = vd->insert_point_site(p0);
    int id1 = vd->insert_point_site(p1);
    int id2 = vd->insert_point_site(p2);
    int id3 = vd->insert_point_site(p3);
    int id4 = vd->insert_point_site(p4);

    vd->insert_line_site(id0, id1);
    vd->insert_line_site(id1, id2);
    vd->insert_line_site(id2, id3);
    vd->insert_line_site(id3, id4);
    vd->insert_line_site(id4, id0);
    vd->check();

    // draw four offsets.
    ovd::HEGraph& g = vd->get_graph_reference();
    ovd::Offset offset(g);
    for (int i=1; i<5; i++) {
      ovd::OffsetLoops offset_list = offset.offset(i*0.008);
      BOOST_FOREACH( ovd::OffsetLoop loop, offset_list ) { // loop through each loop
        bool first = true;
        BOOST_FOREACH( ovd::OffsetVertex lpt, loop ) { // loop through each line/arc
          if (first) {
            first = false;
            std::cout << "first offset:p:" << lpt.p << std::endl;
          } else {
            std::cout << "offset:p:" << lpt.p << ",r:" << lpt.r << ",c:" << lpt.c << ",cw:" << lpt.cw << std::endl;
          }
        }
      }
    }

    std::cout << vd->print();

    return 0;
}
