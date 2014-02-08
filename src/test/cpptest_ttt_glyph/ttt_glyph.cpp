// OpenVoronoi example. Build VD for a single glyph from true-type-tracer.
// true-type-tracer is from https://github.com/aewallin/truetype-tracer

#include <string>
#include <iostream>
#include <algorithm>

#include <boost/foreach.hpp>
#include <boost/timer.hpp>
#include <boost/program_options.hpp>

#include "voronoidiagram.hpp"
#include "checker.hpp"
#include "version.hpp"
#include "common/point.hpp"
#include "utility/vd2svg.hpp"

#include <truetypetracer/ttt.hpp>
#include <truetypetracer/segment_writer.hpp>

#define TTFONT "/usr/share/fonts/truetype/freefont/FreeSerif.ttf"

Loops get_ttt_loops(int char_index=0, double ttt_scale=1e-4, int subdivision=50) {
    SEG_Writer my_writer; // this Writer writes G-code
    my_writer.set_scale( ttt_scale  ); // was 1e-4
    my_writer.set_conic_line_subdiv( subdivision );
    
    std::string all_glyphs = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    char c = all_glyphs.at(char_index);
    std::string glyph(&c,1);
    std::cout << " ttt glyph is " << glyph << "\n";
    
    Ttt t( &my_writer, glyph, false , TTFONT ); // Ttt( Writer*, text, unicode?, path-to-font )
    
    Loops all_loops = my_writer.get_loops();
    
    extents ext = my_writer.get_extents();
    std::cout << " x-range : " << ext.minx << " - " << ext.maxx << "\n";
    std::cout << " y-range : " << ext.miny << " - " << ext.maxy << "\n";  
    
    Loops mod_loops;
    // ttt returns loops with duplicate start/endpoints
    // to avoid duplicates, remove the first point from each loop
    BOOST_FOREACH(Loop l, all_loops) {
        std::reverse(l.begin(), l.end());
        l.erase(l.begin()); // remove first element of loop
        mod_loops.push_back(l);
    }

    return mod_loops;
}

namespace po = boost::program_options;



// ttt example program
int main(int argc,char *argv[]) {
    // Declare the supported options.
    po::options_description desc("This program calculates the voronoi diagram for a glyph from truetype-tracer.\nAllowed options");
    desc.add_options()
        ("help", "produce help message")
        ("c", po::value<int>(), "set character, an int in [0,51] ")
        ("s", po::value<int>(), "set scale (int), [0,5]. ")
        ("d", po::value<int>(), "subdivision (int), 50 ")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    
    int char_index = 0;

    if (vm.count("c")) {
        int c = vm["c"].as<int>();
        if ( 0<=c && c<= 51 )
            char_index = c;
        else {
            std::cout << "--c option out of range!\n";
            return 1;
        }
    }
    
    int subdivision = 50;
    if (vm.count("d")) {
        int d = vm["d"].as<int>();
        if ( 1<=d && d<= 100 )
            subdivision = d;
        else {
            std::cout << "--d option out of range!\n";
            return 1;
        }
    }
    
    // the different scales
    double scales[5] = {2e-4, 1e-4, 2e-5, 1e-5, 2e-6};
    const std::vector<double> scalevector(scales, scales+5);

    double scale = scalevector[0];
    if (vm.count("s")) {
        unsigned int s = vm["s"].as<int>();
        if ( 0<=s && s< scalevector.size() )
            scale = scalevector[s];
        else {
            std::cout << "--s option out of range!\n";
            return 1;
        }
    }
    
    Loops all_loops = get_ttt_loops(char_index, scale, subdivision);
    
    // print out the points
    int nloop =0;
    BOOST_FOREACH(Loop l, all_loops) {
        std::cout << "Loop " << nloop << " has " << l.size() << " points:\n";
        int npoint=0;
        BOOST_FOREACH(Point pt, l) {
            std::cout << " Point " << npoint << " x= " << pt.x << " y= " << pt.y << "\n";
            npoint++;
        }
        std::cout << "End Loop " << nloop << "\n";
        nloop++;
    } 
    
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,50);
    std::cout << "OpenVoronoi version: " << ovd::version() << "\n";
    
    // store the verted IDs in loop_id, for later inserting line-segments
    std::vector< std::vector<int> > loops_id;
    boost::timer tmr;
    int n_points=0;
    BOOST_FOREACH(Loop loop, all_loops) {
        std::vector<int> loop_id;
        BOOST_FOREACH( Point pt, loop) {
            loop_id.push_back( vd->insert_point_site( ovd::Point(pt.x,pt.y)) ); 
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
            std::cout << "Inserting LineSite: " << loop_id[n] << " - " << loop_id[next] << "\n" << std::flush;
            vd->insert_line_site( loop_id[n], loop_id[next]); 
            n_linesites++;
        }
    }
    std::cout << n_linesites << " LineSites:s inserted in  " << tmr.elapsed() << " seconds\n" << std::flush;
    
    // write vd to svg file
    std::string filename = "ttt_glyph.svg";
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
