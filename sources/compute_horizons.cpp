#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <algorithm>

// #include <gnuplot-iostream.h>
#include <assert.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "vector3d.hpp"
#include "swiss_to_lat_lon.hpp"
#include "import_heightmap_horizon.hpp"



namespace po = boost::program_options;


class pos_hoz {
public:
    vector3d pos;
    vector3d norm;
    std::vector<signed short> elevation_angles;
    std::vector<unsigned char> distances;
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&pos, sizeof(vector3d));
        ofs.write((char*)&norm, sizeof(vector3d));
        ofs.write((char*)&elevation_angles[0], sizeof(short)*elevation_angles.size());
        ofs.write((char*)&distances[0], distances.size());
    }
};

bool EQ_DBL(double a, double b){
    return std::abs(a-b)<EPS_DBL;
}





int main(int ac, char** av) {
    //Variables to be assigned by program options
    bool test = false;
    double height_map_resolution;
    double north_bound;
    double south_bound;
    double east_bound;
    double west_bound;
    int horizon_angles;
    bool elevation_dependant_sun_intensity;
    double tile_size;
    
    std::string input_file;
    std::string output_dir;
    bool verbose = false;
    
    // Declare the supported options.
    po::options_description op_desc("Allowed options");
    op_desc.add_options()
    ("help", "print options table")
    ("input-file,i", po::value<std::string>(&input_file), "File containing height map data in x<space>y<space>z<newline> swiss coordinates format")
    ("output-dir,o", po::value<std::string>(&output_dir), "Directory for output files (default data_out_<date_time>)")
    ("resolution,R", po::value<double>(&height_map_resolution)->default_value(200.0), "resolution of data (default: 200.0)")
    ("tile-size, T", po::value<double>(&tile_size)->default_value(40), "Output files will have size tile_size x tile_size (default 40)")
    ("nmax",po::value<double>(&north_bound)->default_value(1e100), "maximum north coordinate to be treated (default 1e100)")
    ("nmin",po::value<double>(&south_bound)->default_value(-1e100), "minimum north coordinate to be treated (default -1e100)")
    ("emax",po::value<double>(&east_bound)->default_value(1e100), "maximum east coordinate to be treated (default 1e100)")
    ("emin",po::value<double>(&west_bound)->default_value(-1e100), "minimum east coordinate to be treated (default -1e100)")
    ("horizonangles,h", po::value<int>(&horizon_angles)->default_value(360), "number of angles to compute horizon for (default 360)")
    ("elevation,E", "Include effect of terrain elevation in the computation of sun intensity")
    ("verbose, v", "Verbose: output lots of text")
    ;
    
    po::positional_options_description pd;
    pd.add("input-file", 1).add("output-file", 1);
    
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, op_desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout <<"CrunchGeoData [options] [input file] [output base]"<<std::endl<< op_desc << std::endl;
        return 1;
    }
    
    if (vm.count("elevation")){
        elevation_dependant_sun_intensity = true;
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
    

    
    std::unordered_map<vector3d, vector3d , hash> grid_points;   //grid_points is unordered_map of all points in bounding box with
    //a vector to hold the normal
    
    double north_x_max;
    double east_y_max;
    double south_x_min;
    double west_y_min;
    
    import_heightmap(input_file, south_bound, north_bound, east_bound, west_bound,
                     horizon_angles, grid_points ,north_x_max, south_x_min, east_y_max, west_y_min);
    
    
    std::pair<double,double> NE = swiss_to_lat_lon(north_x_max+height_map_resolution/2.0, east_y_max+height_map_resolution/2.0);
    std::pair<double,double> SW = swiss_to_lat_lon(south_x_min-height_map_resolution/2.0, west_y_min-height_map_resolution/2.0);
    
    std::cout << "NE: " << north_x_max <<", " << east_y_max << " SW: "<< south_x_min << " , " << west_y_min <<  std::endl;
    std::cout /*<< std::setprecision(9)*/ << "NE: " << NE.first <<", " << NE.second << " SW: "<< SW.first << ", " << SW.second <<  std::endl;
    std::cout << "Number of points in dataset: " << grid_points.size() << std::endl;
    
//    double average_latitude=((NE.first+SW.first)/2.0)*M_PI/180.0;
    
    
    
//std::cout <<"   Progress: "<< 100*(N*1.0)/(grid_points.size()*1.0) <<" %    \n";
    
//    std::ofstream ofs(output_file);    //open output file
    
    int N=0;

    
    std::string file_name;
    std::ofstream ofs;
    time_t start_time  = time(NULL);
    
    for(auto grid_point : grid_points){
        vector3d v = grid_point.first;
        
        if(test){
            double y_test = height_map_resolution*floor(601430.0/height_map_resolution);
            double x_test = height_map_resolution*floor(126243.0/height_map_resolution);
            auto it_v = grid_points.find(vector3d(x_test,y_test,0));      //get the two gridpoint from the set
            std::cout << "Test point: " << x_test <<" "<< y_test << std::endl;
            std::cout << "Resolution: " << height_map_resolution << ", tile size: "<< tile_size << std::endl;
            if (it_v == grid_points.end()){
                std::cout << "Test point not found" << std::endl;
                exit(-3);
            }
            v = it_v->first;
        }
        
        
        N++;
        
        if (1) {
            double run_time = difftime(time(NULL), start_time);
            int rhours = (int)floor(run_time/3600.0);
            int rminutes = (int)floor((run_time-rhours*3600.0)/60.0);
            int rseconds = (int)floor(run_time-rhours*3600.0-rminutes*60.0);
            double remaining_time = ((grid_points.size() - N)*1.0)*run_time/(N*1.0);
            int hours = (int)floor(remaining_time/3600.0);
            int minutes = (int)floor((remaining_time-hours*3600.0)/60.0);
            int seconds = (int)floor(remaining_time-hours*3600.0-minutes*60.0);
            double progress = round(1000.0*100.0*(N*1.0)/(grid_points.size()*1.0))/1000.0;
            std::cout <<"   Progress: "<< progress <<" %  ("<<N<<")  time: "<<rhours<<":"<<rminutes<<":"<<rseconds<<" time left: "<<hours<<":"<<minutes<<":"<<seconds<<"      \r";
            std::cout.flush();
        }
        vector3d vNx(v.x+height_map_resolution*10,v.y,0);            //get points 10x resolution to the north, west, south and east
        vector3d vWx(v.x,v.y+height_map_resolution*10,0);
        vector3d vSx(v.x-height_map_resolution*10,v.y,0);
        vector3d vEx(v.x,v.y-height_map_resolution*10,0);
        auto itEx=grid_points.find(vEx);
        auto itNx=grid_points.find(vNx);
        auto itWx=grid_points.find(vWx);
        auto itSx=grid_points.find(vSx);
        if (itEx == grid_points.end()||itNx == grid_points.end()||itWx == grid_points.end()||itSx == grid_points.end()){
            continue;     //if we are to close to border, dont compute anything, continue with next point
        }
        
        vector3d vN(v.x+height_map_resolution,v.y,0);            //get points to the north, west, south and east
        vector3d vW(v.x,v.y+height_map_resolution,0);
        vector3d vS(v.x-height_map_resolution,v.y,0);
        vector3d vE(v.x,v.y-height_map_resolution,0);
        
        auto itE=grid_points.find(vE);
        auto itN=grid_points.find(vN);
        auto itW=grid_points.find(vW);
        auto itS=grid_points.find(vS);
        
        if (itE == grid_points.end()||itN == grid_points.end()||itW == grid_points.end()||itS == grid_points.end()){
            exit(-1);
        }
        
        vE=itE->first;
        vN=itN->first;
        vW=itW->first;
        vS=itS->first;
        
        vector3d Normal = ((((vE-v)^(vN-v))+((vW-v)^(vS-v)))/2).norm();       //normal is the crossproduct of two prependicular differences. Avgd.
        
        
        std::vector<double> horizon_elevation_angles;   //vector to hold elevation angles at theta
        horizon_elevation_angles.assign(horizon_angles, -M_PI_2);
        std::vector<double> dists(horizon_angles);     //distance to currently min elevation
        std::vector<vector3d> hozpoints(horizon_angles);     //distance to currently min elevation

        
        size_t edge_points = 0;
        size_t interior_points = 0;
        
        
        
        
        
        
        for(int theta = 0; theta<horizon_angles;theta++){        //compute for each angle theta (index of angles)
            double th = theta*M_PI*2.0/horizon_angles;           //angle in radians corresponding to theta index
            double sin_theta = sin(th);   //calculate sin of the angle in radians
            if(!EQ_DBL(sin_theta, 0)){
                double sin_sign = sin_theta>0?1.0:-1.0;
                for(double y=v.y+height_map_resolution*sin_sign;(y >= west_y_min && y <= east_y_max);y+=sin_sign*height_map_resolution){//increase/decrease if theta +/-
                    double x;
                    double x1;
                    double x2;
                    if(EQ_DBL(th, M_PI_2)||EQ_DBL(th, -M_PI_2)){
                        x = v.x;
                        x1 = x;
                        x2 = x;
                    }
                    else{
                        x = v.x+(y-v.y)/tan(th);                             //find the x that goes with y for this theta
                        x1 = height_map_resolution*floor(x/height_map_resolution);    //find the nearest lower gridpoint
                        x2 = x1 + height_map_resolution;                  //Add height_map_resolution to get nearest higher gridpoint
                    }
                    double dist = (v-vector3d(x,y,v.z)).length();
                    double phiMaxTheta = atan(5000.0/dist);      //maximal phi at this theta (with height 5000 m)
                    double phi;
                    double height;
                    
                    if(x1>=north_x_max||x1<south_x_min||x2>=north_x_max||x2<south_x_min||horizon_elevation_angles[theta]>phiMaxTheta){
                        break;
                    }
                    if((int)x%(int)height_map_resolution){                       //if x is not a grid point
                        
                        auto it_vec1 = grid_points.find(vector3d(x1,y,0));      //get the two gridpoints from the set
                        auto it_vec2 = grid_points.find(vector3d(x2,y,0));
                        if (it_vec1 == grid_points.end()||it_vec2 == grid_points.end()){
                            edge_points++;
                            continue;
                        }
                        vector3d vec1 = it_vec1->first;
                        vector3d vec2 = it_vec2->first;
                        height = vec1.z*(x2-x)/height_map_resolution + vec2.z*(x-x1)/height_map_resolution-v.z;       //compute height at x via linear interpolation
                        phi = atan(height/dist);
                        
                        
                    }
                    else{//if x is a gridpoint
                        auto it_vec = grid_points.find(vector3d(x,y,0));  //get vector
                        if (it_vec == grid_points.end()){
                            //exit(-1);
                            edge_points++;
                            continue;
                        }
                        vector3d vec = it_vec->first;
                        height = vec.z-v.z;                          //get height
                        phi = atan(height/dist);
                        }
                    if(phi>horizon_elevation_angles[theta]){//see if larger
                        horizon_elevation_angles[theta]=phi;
                        dists[theta] = dist/1000.0;
                        hozpoints[theta] = vector3d(x,y,height+v.z);
                    }

                    interior_points++;
                }
                
            }
            
            double cos_theta = cos(th);
            if(!EQ_DBL(cos_theta, 0)){
                double cos_sign = cos_theta>0?1.0:-1.0;
                for(double x=v.x+height_map_resolution*cos_sign;(x >= south_x_min && x <= north_x_max);x+=height_map_resolution*cos_sign){//increase/decrease if theta +/-

                    double y = v.y+(x-v.x)*tan(th);                             //find the y that goes with x for this theta
                    double y1 = height_map_resolution*floor(y/height_map_resolution);    //find the nearest lower gridpoint
                    double y2 = y1 + height_map_resolution;                  //Add height_map_resolution to get nearest higher gridpoint
                    double dist = (v-vector3d(x,y,v.z)).length();
                    double phiMaxTheta = atan(5000.0/dist);      //maximal phi at this theta (with height 5000 m)
                    double height;
                    
                    if(y1>=east_y_max||y1<=west_y_min||y2>=east_y_max||y2<=west_y_min||horizon_elevation_angles[theta]>phiMaxTheta){
                        break;
                    }
                    
                    double phi;
                    if((int)y%(int)height_map_resolution){                       //if y is not a grid point
                        auto it_vec1 = grid_points.find(vector3d(x,y1,0));      //get the two gridpoints from the set
                        auto it_vec2 = grid_points.find(vector3d(x,y2,0));
                        if (it_vec1 == grid_points.end()||it_vec2 == grid_points.end()){
                            //exit(-1);
                            edge_points++;
                            continue;
                        }
                        vector3d vec1 = it_vec1->first;
                        vector3d vec2 = it_vec2->first;
                        height = vec1.z*(y2-y)/height_map_resolution + vec2.z*(y-y1)/height_map_resolution-v.z;       //compute height at y via linear interpolation
                        phi = atan(height/dist);
                        
                    }
                    else{//if y is a gridpoint
                        auto it_vec = grid_points.find(vector3d(x,y,0));  //get vector
                        if (it_vec == grid_points.end()){
                            //exit(-1);
                            edge_points++;
                            continue;
                        }
                        vector3d vec = it_vec->first;
                        height = vec.z-v.z;                          //get height
                        phi = atan(height/dist);
                        
                    }
                    if(phi>horizon_elevation_angles[theta]){//see if larger
                        horizon_elevation_angles[theta]=phi;
                        dists[theta] = dist/1000.0;
                        hozpoints[theta] = vector3d(x,y,height+v.z);
                    }

                    
                }
                interior_points++;
            }
            
            
            
            
        }
        std::string tile_name;
        int tile_x = floor(v.x/(tile_size*height_map_resolution))*tile_size*height_map_resolution;
        int tile_y = floor(v.y/(tile_size*height_map_resolution))*tile_size*height_map_resolution;
        tile_name = output_dir + std::string("/tile_") + std::to_string(tile_x) + std::string("_") + std::to_string(tile_y) + std::string(".hoz");
        std::cout << "tile_name " << tile_name << std::endl;
        //std::string("./") +
        
        
        
        if(file_name.compare(tile_name)!=0){
            if(ofs.is_open())
                ofs.close();
            file_name = tile_name;
            std::cout << "file_name " << file_name << std::endl;
            ofs.open(file_name, std::iostream::out | std::iostream::app | std::iostream::binary);
   //open output file
            if (!ofs.is_open()){
                std::cout << "Can't open output file " << file_name << std::endl;
                exit(-2);
            }
        }
        pos_hoz outputs;
        outputs.pos = v;
        outputs.norm = Normal;
//        ofs << v.x << v.y << v.z;
//        ofs << Normal;
        for(int k=0;k<horizon_angles;k++){
            outputs.elevation_angles.push_back(horizon_elevation_angles[k]*32767/(M_PI/2));
            char dst = (char) (dists[k]>255?255:dists[k]);
            outputs.distances.push_back(dst);
 //           ofs << horizon_elevation_angles[k];
  //          ofs << dists[k];
            //ofs << hozpoints[k] << std::endl;
        }
        outputs.write(ofs);
        //ofs.write((char*)outputs, (horizon_angles*2+6)*sizeof(double));
        if(test){
            break;
        }
        ofs.flush();
        //ofs << std::endl;
        
    }
    std::cout << grid_points.size() << std::endl;
    return 0;
}
