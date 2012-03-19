// OpenVoronoi example. Build VD for a single glyph from true-type-tracer.
// true-type-tracer is from https://github.com/aewallin/truetype-tracer

#include <string>
#include <iostream>

#include <boost/foreach.hpp>
#include <boost/timer.hpp>

#include <openvoronoi/voronoidiagram.hpp>
#include <openvoronoi/checker.hpp>
#include <openvoronoi/version.hpp>
#include <openvoronoi/common/point.hpp>
#include <openvoronoi/utility/vd2svg.hpp>

#include <ttt/ttt.hpp>
#include <ttt/segment_writer.hpp>

#define TTFONT "/usr/share/fonts/truetype/freefont/FreeSerif.ttf"

Loops get_ttt_loops() {
    SEG_Writer my_writer; // this Writer writes G-code
    my_writer.set_scale(  2.5095362377e-4 );
    Ttt t( &my_writer, "P", false , TTFONT ); // ( Writer*, text, unicode?, path-to-font )
    
    Loops all_loops = my_writer.get_loops();
    
    extents ext = my_writer.get_extents();
    std::cout << " x-range : " << ext.minx << " - " << ext.maxx << "\n";
    std::cout << " y-range : " << ext.miny << " - " << ext.maxy << "\n";  
    
    Loops mod_loops;
    // ttt returns loops with duplicate start/endpoints
    // to avoid duplicates, remove the first point from each loop
    BOOST_FOREACH(Loop l, all_loops) {
        l.erase(l.begin()); // remove first element of loop
        mod_loops.push_back(l);
    }

    return mod_loops;
}

// ttt example program
int main() {
    
    Loops all_loops = get_ttt_loops();
    
    // print out the points
    int nloop =0;
    BOOST_FOREACH(Loop l, all_loops) {
        std::cout << "Loop " << nloop << " has " << l.size() << " points:\n";
        int npoint=0;
        BOOST_FOREACH(Point pt, l) {
            std::cout << " Point " << npoint << " x= " << pt.first << " y= " << pt.second << "\n";
            npoint++;
        }
        std::cout << "End Loop " << nloop << "\n";
        nloop++;
    } 
    
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,50);
    std::cout << "OpenVoronoi version: " << ovd::version() << "\n";
    
    // store the verted IDs here, for later inserting line-segments
    std::vector< std::vector<int> > loops_id;
    typedef std::pair<double,double> ttt_pt;
    boost::timer tmr;
    int n_points=0;
    BOOST_FOREACH(Loop loop, all_loops) {
        std::vector<int> loop_id;
        BOOST_FOREACH(ttt_pt pt, loop) {
            loop_id.push_back( vd->insert_point_site( ovd::Point(pt.first,pt.second)) ); 
            n_points++;
        }
        loops_id.push_back(loop_id);
    }
    std::cout << n_points << " PointSite:s inserted in  " << tmr.elapsed() << " seconds\n" << std::flush;
    
    // now we insert LineSite:s
    tmr.restart();
    int n_linesites=0;
    BOOST_FOREACH( std::vector<int> loop_id, loops_id) {
        for( unsigned int n=0; n<loop_id.size();n++) {
            int next = n+1;
            if (n==loop_id.size()-1)
                next = 0;
            vd->insert_line_site( loop_id[n], loop_id[next]); 
            n_linesites++;
        }
    }
    std::cout << n_linesites << " LineSites:s inserted in  " << tmr.elapsed() << " seconds\n" << std::flush;
    
    // write vd to svg file
    std::string filename = "ttt_glypoh.svg";
    vd2svg(filename, vd);
    std::cout << "Output written to SVG file: " << filename << "\n"; //tmr.elapsed() << " seconds\n" << std::flush;
    
    // sanity-check
    ovd::VoronoiDiagramChecker check(vd->get_graph_reference());
    if (!check.is_valid())
        return -1;
    else
        std::cout << "ovd::VoronoiDiagramChecker OK!\n";
    
    delete vd;
    return 0;
}
