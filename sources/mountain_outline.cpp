#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
//#include <map>
//#include <set>
#include <unordered_set>
//#include <unordered_map>
#include <cmath>
#include <algorithm>

// #include <gnuplot-iostream.h>
#include <assert.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "vector3d.hpp"
#include "swiss_to_lat_lon.hpp"
#include "height_map.hpp"

const int horizon_angles = 360;
const int tile_size = 40;


namespace po = boost::program_options;

struct pos_hoz {
    vector3d pos;
    vector3d norm;
    short elevation_angles[horizon_angles];
    pos_hoz() : pos(0,0,0), norm(0,0,0) {
        memset(elevation_angles, 0, sizeof(short) * horizon_angles);
    }
};

bool EQ_DBL(double a, double b){
    return std::abs(a-b)<EPS_DBL;
}

int main(int ac, char** av) {
    
    bool test = true;
    //Variables to be assigned by program options
    double height_map_resolution = 25.0;;
    double north_bound = 1e100;
    double south_bound = -1e100;
    double east_bound = 1e100;
    double west_bound = -1e100;
    
    
    std::string input_file;
    std::string output_file;
    bool verbose = false;
    
    /*
    // Declare the supported options.
    po::options_description op_desc("Allowed options");
    
    op_desc.add_options()
    ("help", "print options table")
    ("input-file,i", po::value<std::string>(&input_file), "File containing height map data in x<space>y<space>z<newline> swiss coordinates format")
    ("output-dir,o", po::value<std::string>(&output_file), "File for output data (default data_out.hoz)")
    ("resolution,R", po::value<double>(&height_map_resolution)->default_value(25.0), "resolution of data (default: 25.0)")
    ("nmax",po::value<double>(&north_bound)->default_value(1e100), "maximum north coordinate to be treated (default 1e100)")
    ("nmin",po::value<double>(&south_bound)->default_value(-1e100), "minimum north coordinate to be treated (default -1e100)")
    ("emax",po::value<double>(&east_bound)->default_value(1e100), "maximum east coordinate to be treated (default 1e100)")
    ("emin",po::value<double>(&west_bound)->default_value(-1e100), "minimum east coordinate to be treated (default -1e100)")
    ("verbose, v", "Verbose: output lots of text")
    ;
    
    po::positional_options_description pd;
    pd.add("input-file", 1).add("output-file", 1);
    
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, op_desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout <<"MountainOutline [options] [input file] [output base]"<<std::endl<< op_desc << std::endl;
        return 1;
    }
    
    if (vm.count("verbose")){
        verbose = true;
    }
    if(!vm.count("input-file")){
        std::cout << "Input file must be specified" << std::endl;
        exit(255);
    }
    if(!vm.count("output-dir")){
        output_file = std::string("data_out.hoz");
        std::cout << "Output directory will be: " << std::endl;
    }
    */
    input_file = std::string("/Users/alexxx/HeightMaps/precieous/DHM25.xyz");
    std::ifstream ifs(input_file);
    if (!ifs.is_open())
        throw std::runtime_error("could not open file : " + std::string(input_file));
    
    height_map grid_points(ifs, south_bound, north_bound, east_bound, west_bound);
    
    
    //grid_points is a height_map, an unordered_set of all points in bounding box with
    //the maximum and minimum x,y,h of all values stored.
    
    
    std::pair<double,double> NE = swiss_to_lat_lon(grid_points.xmax()+height_map_resolution/2.0, grid_points.ymax()+height_map_resolution/2.0);
    std::pair<double,double> SW = swiss_to_lat_lon(grid_points.xmin()-height_map_resolution/2.0, grid_points.ymin()-height_map_resolution/2.0);
    
    std::cout << "NE: " << grid_points.xmax() <<", " << grid_points.ymax() << " SW: "<< grid_points.xmin() << " , " << grid_points.ymin() <<  std::endl;
    std::cout  << "NE: " << NE.first <<", " << NE.second << " SW: "<< SW.first << ", " << SW.second <<  std::endl;
    std::cout << "Maximum height: " << grid_points.hmax() << std::endl;
    std::cout << "Number of points in dataset: " << grid_points.size() << std::endl;
    
    if(test){
        double y_test = height_map_resolution*floor(601430.0/height_map_resolution);
        double x_test = height_map_resolution*floor(126243.0/height_map_resolution);
        auto it_v = grid_points.find_point(vector3d(x_test,y_test,0));      //get the gridpoint from the set
        std::cout << "Test point: " << x_test <<" "<< y_test << std::endl;
        std::cout << "Resolution: " << height_map_resolution << ", tile size: "<< tile_size << std::endl;
        vector3d v = *it_v;
        double elevation_angles[360];
        int theta_start = 0;
        int theta_end = 360;
        
        //        time_t start_time  = time(NULL);
        clock_t t = clock();
        for(int theta = theta_start;theta<theta_end;theta++){
            double phi = atan(grid_points.compute_elevation_angle(v,theta, height_map_resolution, 360.0/horizon_angles));
            short phi_short = (short)((phi / M_PI_2) * std::numeric_limits<short>::max());
            elevation_angles[theta] = phi_short;
            
        }
        double run_time = ((double)(clock()-t))/CLOCKS_PER_SEC;
        //        double run_time = difftime(time(NULL), start_time);
        std::cout << "phi = [";
        for(int theta = theta_start;theta<theta_end;theta++){
            std::cout << elevation_angles[theta] << ";";
        }
        std::cout << "];" << std::endl;
        std::cout << "runtime test point: " << run_time << std::endl;
        
        return 0;
    }
    
    

    return 0;
}






