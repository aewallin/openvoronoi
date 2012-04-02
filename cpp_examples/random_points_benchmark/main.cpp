#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <numeric>

#include <openvoronoi/voronoidiagram.hpp>
#include <openvoronoi/version.hpp>

#include <boost/random.hpp>
#include <boost/timer.hpp>

const char *BENCHMARK_FILE = "voronoi_benchmark.txt";

boost::mt19937 gen(static_cast<unsigned int>(time(NULL)));
boost::timer timer;

typedef ovd::Point point_type;
typedef double coordinate_type;

void construct_voronoi_points(std::vector<point_type>& points) {
    int bins =  (int)(sqrt(2)*sqrt(points.size()));  // number of bins for nearest-neighbor bin-search
    ovd::VoronoiDiagram* vd = new ovd::VoronoiDiagram(1,bins);
    for(unsigned int n=0;n<points.size();++n)
        vd->insert_point_site( points[n] ); // insert each point. This returns an int-handle which we do not use here.
    
    assert( vd->check() ); // this runs a sanity-check on the diagram. This is slow, so run only in debug mode.
    
    delete vd;
}

int main(int argc,char *argv[]) {

    coordinate_type minimum_coordinate = std::numeric_limits<uint32_t>::min();
    coordinate_type maximum_coordinate = std::numeric_limits<uint32_t>::max();
    
    std::cout << "OpenVoronoi " << ovd::version() << " " << ovd::build_type() << "\n";
    std::cout << "| Number of points | Number of tests | Time per one test |  usec/n*log2(n)   |" << std::endl;
    
    std::ofstream bench_file(BENCHMARK_FILE, std::ios_base::out | std::ios_base::app);
    bench_file << "Voronoi Benchmark Test (time in seconds):" << std::endl;
    bench_file << "| Number of points | Number of tests | Time per one test |" << std::endl;
    bench_file << std::setiosflags(std::ios::right | std::ios::fixed) << std::setprecision(6);

    #ifdef NDEBUG
        int max_points = 100000;
    #else
        int max_points = 10000;
    #endif
    // OpenVoronoi timings on i7-2600K (roughly)
    // 1k    0.1s
    // 10k   1.75s
    // 100k  19s
    
    std::vector<point_type> points;
    coordinate_type x, y;
    for (int num_points = 10; num_points <= max_points; num_points *= 10 ) {
        points.resize(num_points);
        
        int num_tests = max_points / num_points;
        
        std::vector< std::vector<point_type> > test_points;
        test_points.resize(num_tests);
        for (int cur = 0; cur < num_tests; cur++) {
            for (int cur_point = 0; cur_point < num_points; cur_point++) {
                x = static_cast<coordinate_type>(gen());
                y = static_cast<coordinate_type>(gen());

                // For OpenVoronoi we scale coordinates so they fit within a unit circle.
                // e.g. a box centered at (0,0) with side-length 0.8
                x -= (maximum_coordinate-minimum_coordinate); // center around x=0
                y -= (maximum_coordinate-minimum_coordinate); // center around y=0
                x *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
                y *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
                assert( (-0.4 <= x) && ( x <= +0.4 ) );
                assert( (-0.4 <= y) && ( y <= +0.4 ) );
                
                points[cur_point] = point_type(x, y);
            }
            test_points[cur] = points;
        }
        
        timer.restart();
        for (int cur = 0; cur < num_tests; cur++) {
            construct_voronoi_points( test_points[cur] );
        }
        double elapsed_time = timer.elapsed();
        double time_per_test = elapsed_time / num_tests;
        std::cout << "| " << std::setw(16) << num_points << " ";
        std::cout << "| " << std::setw(15) << num_tests << " ";
        std::cout << "| " << std::setw(17) << time_per_test << " ";
        std::cout << "| " << std::setw(17) << 1e6*time_per_test/(num_points*log(num_points)/log(2)) << " ";
        std::cout << "| " << std::endl << std::flush;
        
        bench_file << "| " << std::setw(16) << num_points << " ";
        bench_file << "| " << std::setw(15) << num_tests << " ";
        bench_file << "| " << std::setw(17) << time_per_test << " ";
        bench_file << "| " << std::endl;
    }
    bench_file.close();
    return 0;
}

