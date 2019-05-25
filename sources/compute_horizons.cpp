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
    
    bool test = false;
    //Variables to be assigned by program options
    double height_map_resolution;
    double north_bound;
    double south_bound;
    double east_bound;
    double west_bound;
    
    
    std::string input_file;
    std::string output_dir;
    bool verbose = false;
    
    // Declare the supported options.
    po::options_description op_desc("Allowed options");
    op_desc.add_options()
    ("help", "print options table")
    ("input-file,i", po::value<std::string>(&input_file), "File containing height map data in x<space>y<space>z<newline> swiss coordinates format")
    ("output-dir,o", po::value<std::string>(&output_dir), "Directory for output files (default data_out_<date_time>)")
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
        std::cout <<"ComputeHorizons [options] [input file] [output base]"<<std::endl<< op_desc << std::endl;
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
        std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string time_string(30, 0);
        std::strftime(&time_string[0], time_string.size(), "%Yx%mx%dT%Hx%Mx%S", std::localtime(&now));
        output_dir = std::string("hoz_out_")+time_string;
        std::cout << "Output directory will be: "+output_dir << std::endl;
    }
    
    //create output directory
    boost::filesystem::path dir(output_dir);
    if(!boost::filesystem::create_directory(dir)) {
        std::cout << "Failed to creat output directory" << "\n";
        exit(255);
    }
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
            elevation_angles[theta] = grid_points.compute_elevation_angle(v, theta, height_map_resolution, 1.0);
        }
        double run_time = ((double)(clock()-t))/CLOCKS_PER_SEC;
        //        double run_time = difftime(time(NULL), start_time);
        std::cout << "phi = [";
        for(int theta = theta_start;theta<theta_end;theta++){
            std::cout << elevation_angles[theta] << ";";
        }
        std::cout << "];" << std::endl;
        std::cout << "runtime test point: " << run_time << std::endl;
        std::cout << "estimated runtime all points: " << run_time*grid_points.size() << std::endl;
        
        return 0;
    }
    
    
    
    
    std::string file_name;
    std::ofstream ofs;
    
    time_t start_time  = time(NULL);
    int N=0;
    int NT=0;
    
    int first_tile_y = (int) (grid_points.ymin()+10.0*height_map_resolution);
    int first_tile_x = (int) (grid_points.xmin()+10.0*height_map_resolution);
    int N_tiles_x = (int) ceil((grid_points.xmax()-grid_points.xmin()-20.0*height_map_resolution)/tile_size);
    int N_tiles_y = (int) ceil((grid_points.ymax()-grid_points.ymin()-20.0*height_map_resolution)/tile_size);
    
    pos_hoz tile_points[tile_size*tile_size];
    bool tile_empty = true;
    
    for(int x_tile = 0;x_tile<N_tiles_x; x_tile++){
        for(int y_tile = 0; y_tile<N_tiles_y;y_tile++){
            double coord_x = (first_tile_x + height_map_resolution*tile_size*x_tile);
            double coord_x_end = coord_x+height_map_resolution*tile_size;
            int tile_point = 0;
            tile_empty = true;
            while(coord_x < coord_x_end){
                double coord_y = (first_tile_y + height_map_resolution*tile_size*y_tile);
                double coord_y_end = coord_y+height_map_resolution*tile_size;
                while(coord_y < coord_y_end){
                    
                    //get point and points to the north, west, south and east
                    vector3d v(coord_x, coord_y, 0.0);
                    vector3d vN(v.x+height_map_resolution,v.y,0);
                    vector3d vW(v.x,v.y+height_map_resolution,0);
                    vector3d vS(v.x-height_map_resolution,v.y,0);
                    vector3d vE(v.x,v.y-height_map_resolution,0);
                    
                    auto itv=grid_points.find_point(v);
                    auto itE=grid_points.find_point(vE);
                    auto itN=grid_points.find_point(vN);
                    auto itW=grid_points.find_point(vW);
                    auto itS=grid_points.find_point(vS);
                    
                    if (grid_points.is_end(itv) ||grid_points.is_end(itE) || grid_points.is_end(itN) || grid_points.is_end(itW) || grid_points.is_end(itS)){
                        coord_y = coord_y + height_map_resolution;
                        tile_points[tile_point].pos = vector3d(0,0,0);
                        tile_point++;
                        continue;
                        //if one of the neighboring points are not found we continue with the next point
                    }
                    v = *itv;
                    vE = *itE;
                    vN = *itN;
                    vW = *itW;
                    vS = *itS;
                    tile_points[tile_point].pos = v;
                    
                    vector3d Normal = ((((vE-v)^(vN-v))+((vW-v)^(vS-v)))/2).norm();       //normal is the crossproduct of two prependicular differences. Avgd.
                    
                    tile_points[tile_point].norm = Normal;
                    
                    for(int theta=0;theta<horizon_angles;theta++){
                        double phi = grid_points.compute_elevation_angle(v,theta, height_map_resolution, 360.0/horizon_angles);
                        tile_points[tile_point].elevation_angles[theta] = (short)((phi / M_PI_2) * std::numeric_limits<short>::max());
                    }
                    tile_empty = false;
                    N++;
                    coord_y = coord_y + height_map_resolution;
                    tile_point++;
                    if (N%40==0) {
                        double run_time = difftime(time(NULL), start_time);
                        int rdays = (int)floor(run_time/86400.0);
                        int rhours = (int)floor((run_time-rdays*86400)/3600.0);
                        int rminutes = (int)floor((run_time-rdays*86400.0-rhours*3600.0)/60.0);
                        int rseconds = (int)floor(run_time-rdays*86400.0-rhours*3600.0-rminutes*60.0);
                        double remaining_time = ((grid_points.size() - N)*1.0)*run_time/(N*1.0);
                        int days = (int)floor(remaining_time/86400.0);
                        int hours = (int)floor((remaining_time-days*86400)/3600.0);
                        int minutes = (int)floor((remaining_time-days*86400-hours*3600.0)/60.0);
                        int seconds = (int)floor(remaining_time-days*86400-hours*3600.0-minutes*60.0);
                        double progress = round(1000.0*100.0*(N*1.0)/(grid_points.size()*1.0))/1000.0;
                        if(N > 0){
                            std::cout <<"   Progress: "<< progress <<" %  ("<<N<<", "<<NT<<")  time: "<<rdays<<"d"<<rhours<<":"<<rminutes<<":"<<rseconds<<" time left: "<<days<<"d"<<hours<<":"<<minutes<<":"<<seconds<<"      \r";
                            std::cout.flush();
                        }
                    }
                    
                }
                coord_x = coord_x + height_map_resolution;
            }
            if(!tile_empty){
                std::string tile_name;
                int tile_x = (first_tile_x + height_map_resolution*tile_size*x_tile);
                int tile_y = (first_tile_y + height_map_resolution*tile_size*y_tile);
                tile_name = output_dir + std::string("/tile_") + std::to_string(tile_x) + std::string("_") + std::to_string(tile_y) + std::string(".hoz");
                //            std::cout << std::endl << "tile_name " << tile_name << std::endl;
                ofs.open(tile_name, std::iostream::out | std::iostream::binary);
                //            std::cout << "first point" << tile_points[0].pos << "norm: " << tile_points[0].norm<< std::endl;
                //open output file
                if (!ofs.is_open()){
                    std::cout << "Can't open output file " << tile_name << std::endl;
                }
                ofs.write((char*)&tile_points[0], sizeof(pos_hoz)*tile_size*tile_size);
                ofs.flush();
                ofs.close();
                NT++;
            }
            
        }
    }
    return 0;
}






