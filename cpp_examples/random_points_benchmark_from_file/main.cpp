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
const char *POINT_INPUT_FILE = "voronoi_point.txt";
const int POINT_RUNS = 10;

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
    std::vector<point_type> points;
    {
        std::ifstream input_file(POINT_INPUT_FILE);
        if (input_file.fail()) {
            std::cout << "Failed to open input file!\n";
            exit(-1);
        }
        int num_points;
        coordinate_type x, y;
        input_file >> num_points;
        points.reserve(num_points);
        for (int i = 0; i < num_points; ++i) {
            input_file >> x >> y;
            coordinate_type minimum_coordinate = std::numeric_limits<uint32_t>::min();
            coordinate_type maximum_coordinate = std::numeric_limits<uint32_t>::max();
            // For OpenVoronoi we scale coordinates so they fit within a unit circle.
            // e.g. a box centered at (0,0) with side-length 0.8
            x -= (maximum_coordinate-minimum_coordinate); // center around x=0
            y -= (maximum_coordinate-minimum_coordinate); // center around y=0
            x *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
            y *= 0.4/(maximum_coordinate-minimum_coordinate); // scale to fit within [-0.4, 0.4]
            assert( (-0.4 <= x) && ( x <= +0.4 ) );
            assert( (-0.4 <= y) && ( y <= +0.4 ) );
            
            points.push_back(point_type(x, y));
        }
        input_file.close();
    }
    
    std::vector<double> periods;
    {
        for (int i = 0; i < POINT_RUNS; ++i) {
            timer.restart();
            construct_voronoi_points(points);
            periods.push_back(timer.elapsed());
            std::cout << " run " << i << " / " << POINT_RUNS-1 << " done in " << periods[i] << " s \n" << std::flush;
        }
    }
    std::sort(periods.begin(), periods.end());
    // Using olympic system to evaluate average time.
    double elapsed_time = std::accumulate(periods.begin() + 2, periods.end() - 2, 0.0);
    elapsed_time /= (periods.size() - 4);
    std::cout << "Static test of " << points.size() << " points: " << elapsed_time << std::endl;
    
    std::ofstream bench_file(BENCHMARK_FILE, std::ios_base::out | std::ios_base::app);
    bench_file << std::setiosflags(std::ios::right | std::ios::fixed) << std::setprecision(6);
    bench_file << "Static test of " << points.size() << " points: " << elapsed_time << std::endl;
    bench_file.close();
    
    return 0;
}

