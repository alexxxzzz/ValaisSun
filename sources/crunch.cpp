#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <gnuplot-iostream.h>
#include <assert.h>
#include <boost/program_options.hpp>

#include "vector3d.hpp"
#include "swiss_to_lat_lon.hpp"
#include "import_heightmap.hpp"

namespace po = boost::program_options;



bool EQ_DBL(double a, double b){
    return std::abs(a-b)<EPS_DBL;
}





int main(int ac, char** av) {
    //Variables to be assigned by program options
    double height_map_resolution;
    double north_bound;
    double south_bound;
    double east_bound;
    double west_bound;
    int sun_angles;
    int times_per_year;
    int start_day;
    double summer_angle_panel;
    double winter_angle_panel;
    bool use_terrain_normals = true;
    bool compute_best_fixed_angle = false;
    bool compute_best_summer_winter_angle = false;
    bool compute_horizon = false;

    
    std::string input_file;
    std::string output_file;
    bool verbose = false;
    bool elevation_dependant_sun_intensity = false;
    
    // Declare the supported options.
    po::options_description op_desc("Allowed options");
    op_desc.add_options()
    ("help", "print options table")
    ("input-file,i", po::value<std::string>(&input_file), "File containing height map data in x<space>y<space>z<newline> swiss coordinates format")
    ("output-file-base,o", po::value<std::string>(&output_file)->default_value("data_out"), "Base filname for output files (default data_out)")
    ("resolution,R", po::value<double>(&height_map_resolution)->default_value(200.0), "resolution of data (default: 200.0)")
    ("nmax",po::value<double>(&north_bound)->default_value(1e100), "maximum north coordinate to be treated (default 1e100)")
    ("nmin",po::value<double>(&south_bound)->default_value(-1e100), "minimum north coordinate to be treated (default -1e100)")
    ("emax",po::value<double>(&east_bound)->default_value(1e100), "maximum east coordinate to be treated (default 1e100)")
    ("emin",po::value<double>(&west_bound)->default_value(-1e100), "minimum east coordinate to be treated (default -1e100)")
    ("elevation,E", "Include effect of terrain elevation in the computation of sun intensity")
    ("sunangles,s", po::value<int>(&sun_angles)->default_value(360), "number of angles to compute sunlight from (default 360)")
    ("times,t", po::value<int>(&times_per_year)->default_value(12), "number of times per year to calculate at (default 12)")
    ("start-day, D", po::value<int>(&start_day)->default_value(20), "Day of the year to output first data. (default 20)")
    ("fixed-angle, F", po::value<double>(&summer_angle_panel), "Use fixed angle for solar panel inclination to compute sun intensity.")
    ("summer-angle", po::value<double>(&summer_angle_panel), "Fixed angle for solar panel inclination during summer to compute sun intensity..")
    ("winter-angle", po::value<double>(&winter_angle_panel), "Fixed angle for solar panel inclination during winter to compute sun intensity..")
    ("best-fixed-angle, B", "Use best fixed angle for this latitude to compute sun intensity.")
    ("best-two-season-angles, BSW", "Use best summer and winter angle for this latitude to compute sun intensity.")
    ("terrain-normals, T", "Use terrain normal for each point to compute sun intensity, default option if none of the solar panel angles options are specified.")
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
    if(vm.count("fixed-angle")){
        winter_angle_panel = summer_angle_panel;
        use_terrain_normals = false;
    }
    if(vm.count("summer-angle") && !vm.count("winter-angle")){
        std::cout << "Please specify also a winter angle, or use option fixed-angle" << std::endl;
        exit(255);
    }
    if(vm.count("best-fixed-angle")){
        use_terrain_normals = false;
        compute_best_fixed_angle = true;
    }
    if(vm.count("best-two-season-angles")){
        use_terrain_normals = false;
        compute_best_summer_winter_angle = true;
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
    
    
    std::unordered_map<vector3d, std::vector<double>, hash> grid_points;   //grid_points is unordered_map of all points in bounding box with
                                                                            //a vector to hold average sun power for the days sun is computed
    
    double north_x_max;
    double east_y_max;
    double south_x_min;
    double west_y_min;
    
    import_heightmap(input_file, south_bound, north_bound, east_bound, west_bound,
                     times_per_year, grid_points ,north_x_max, south_x_min, east_y_max, west_y_min);
    
    
    std::pair<double,double> NE = swiss_to_lat_lon(north_x_max+height_map_resolution/2.0, east_y_max+height_map_resolution/2.0);
    std::pair<double,double> SW = swiss_to_lat_lon(south_x_min-height_map_resolution/2.0, west_y_min-height_map_resolution/2.0);
    
    std::cout << "NE: " << north_x_max <<", " << east_y_max << " SW: "<< south_x_min << " , " << west_y_min <<  std::endl;
    std::cout << std::setprecision(9) << "NE: " << NE.first <<", " << NE.second << " SW: "<< SW.first << ", " << SW.second <<  std::endl;
    std::cout << "Number of points in dataset: " << grid_points.size() << std::endl;
    
    double average_latitude=((NE.first+SW.first)/2.0)*M_PI/180.0;
    vector3d summer_normal_panel;
    vector3d winter_normal_panel;

    if(!use_terrain_normals){
        if(compute_best_fixed_angle){
        }
        if(compute_best_summer_winter_angle){
        }
        if(average_latitude > 0){  //if latitude is greater than zero, point south
            summer_normal_panel.x = -1*sin(summer_angle_panel*M_PI/180);
            summer_normal_panel.y = 0;
            summer_normal_panel.z = 1*cos(summer_angle_panel*M_PI/180);

            winter_normal_panel.x = -1*sin(winter_angle_panel*M_PI/180);
            winter_normal_panel.y = 0;
            winter_normal_panel.z = 1*cos(winter_angle_panel*M_PI/180);
        }
        else{  //if latitude is less than zero, point north
            summer_normal_panel.x = 1*sin(summer_angle_panel*M_PI/180);
            summer_normal_panel.y = 0;
            summer_normal_panel.z = 1*cos(summer_angle_panel*M_PI/180);
            
            winter_normal_panel.x = 1*sin(winter_angle_panel*M_PI/180);
            winter_normal_panel.y = 0;
            winter_normal_panel.z = 1*cos(winter_angle_panel*M_PI/180);
        }

    }
    
    
    std::vector<std::vector<double> > sun_elevation_angle(times_per_year,std::vector<double>(sun_angles));
    //Elevation angle PHI of sun for a day and hour
    
    std::vector<std::vector<double> > sun_intensity_day_angle(times_per_year,std::vector<double>(sun_angles));
    
    std::vector<std::vector <vector3d> > sun_vec_day_angle(times_per_year,std::vector<vector3d>(sun_angles));
    //Unit vectors for sun's direction in day, hour.
    
    
    for(int day = 0; day < times_per_year; day++){
        double N=(365.0*day)/times_per_year+start_day;  //We are computing values for the Nth day of the year
        double phi_axis=-asin(0.39779*cos((0.98565*(N+10)+1.914*sin(0.98565*(N-2)*M_PI/180))*M_PI/180));  //Earth axis inclination for day N
        for(int thH=0;thH<sun_angles;thH++){//Angle theta for sun; 0 is North, +90 is West.
            int indTH= thH;
            double thHrad = (180-thH)*M_PI/180.0;  //turn around for sun formula angle
            double phiElv=asin(sin(average_latitude)*sin(phi_axis)+cos(average_latitude)*cos(thHrad)*cos(phi_axis));  //Sun elevation, imperical fomula, see wikipedia
            
            thHrad = M_PI-thHrad;          //In our coordinate system (North x+, West y+) we must inverse the sun theta. Theta corresponds to sun ray direction.
            
            vector3d sun_vec(cos(phiElv) * cos(thHrad) , cos(phiElv) * sin(thHrad) , sin(phiElv));  //Unit vector for direction sun ray come from
            sun_elevation_angle[(int)day][indTH]=phiElv;  //stick sun elevation in its datastructure
            
            sun_vec_day_angle[(int)day][indTH]=sun_vec;  //stick sun vector in its datastructure
            
            double air_mass_coefficient =sqrt(708.0*708.0*sin(phiElv)*sin(phiElv)+2.0*708.0+1.0)-708.0*sin(phiElv);  //wikipedia
            sun_intensity_day_angle[(int)day][indTH]=pow(0.7,pow(air_mass_coefficient,0.678))/0.7;  //stick sun intensity in its datastructure
        }
    }
    
    
    
    
    std::vector<double> max_sun_intensity;
    max_sun_intensity.assign(times_per_year, 0);
    std::vector<double> min_sun_intensity;
    min_sun_intensity.assign(times_per_year, 1e12);
    
    
    int N=0;
    
    std::cout <<"   Progress: "<< 100*(N*1.0)/(grid_points.size()*1.0) <<" %    \n";
    
    for(auto grid_point : grid_points){
        vector3d v = grid_point.first;
        N++;
        
        if (1) {
            std::cout <<"   Progress: "<< 100.0*(N*1.0)/(grid_points.size()*1.0) <<" %  ("<<N<<")    \r";
            std::cout.flush();
        }
        vector3d vNx(v.x+height_map_resolution*10,v.y,0);            //get points 10xresolution to the north, west, south and east
        vector3d vWx(v.x,v.y+height_map_resolution*10,0);
        vector3d vSx(v.x-height_map_resolution*10,v.y,0);
        vector3d vEx(v.x,v.y-height_map_resolution*10,0);
        auto itEx=grid_points.find(vEx);
        auto itNx=grid_points.find(vNx);
        auto itWx=grid_points.find(vWx);
        auto itSx=grid_points.find(vSx);
        if (itEx == grid_points.end()||itNx == grid_points.end()||itWx == grid_points.end()||itSx == grid_points.end()){
            std::vector<double> L;
            L.assign(times_per_year,-1.0);  //border points get value -1.0 to indicate we have coputed no value for them
            grid_points[v]=L;
            continue;
        }
        
        //v = *it_v;
        
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
        horizon_elevation_angles.assign(sun_angles, 0);
        //    std::vector<vector3d> horizon(sun_angles);
        std::vector<double> dists(sun_angles);
        
        size_t edge_points = 0;
        size_t interior_points = 0;
        
        for(int theta = 0; theta<sun_angles;theta++){        //compute for each angle theta
            double th = theta*M_PI/180;
            {
                double sin_theta = sin(th);   //calculate sin of the angle in radians
                double yend = EQ_DBL(copysign(1,sin_theta),1)?east_y_max:0;     //if the sinus is positive endvalue is east_y_max, if negative 0
                for(double y=v.y+copysign(height_map_resolution,sin_theta);(y >= 0 && y <= east_y_max);y+=copysign(height_map_resolution,sin_theta)){//increase/decrease if theta +/-
                    double x = v.x+(y-v.y)/tan(th);                             //find the x that goes with y for this theta
                    double x1 = floor(x-((int)x%(int)height_map_resolution));       //find the nearest lower gridpoint by subtracting remainder according to height_map_resolution
                    double x2 = x1 + height_map_resolution;                  //Add height_map_resolution to get nearest higher gridpoint
                    double phi = horizon_elevation_angles[theta];
                    double dist = (v-vector3d(x,y,v.z)).length();
                    double phiMaxTheta = atan(5000.0/dist);      //maximal phi at this theta (with height 5000 m)
                    
                    
                    if(x1>=north_x_max||x1<0||x2>=north_x_max||x2<0||phi>phiMaxTheta){
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
                        double height = vec1.z*(x2-x)/height_map_resolution + vec2.z*(x-x1)/height_map_resolution-v.z;       //compute height at x via linear interpolation
                        dist=(v-vector3d(x,y,v.z)).length();
                        phi = atan(height/dist);

                        if(phi>horizon_elevation_angles[theta]){//see if larger
                            horizon_elevation_angles[theta]=phi;
                            dists[theta] = dist/1000;
                        }
                        
                    }
                    else{//if x is a gridpoint
                        auto it_vec = grid_points.find(vector3d(x,y,0));  //get vector
                        if (it_vec == grid_points.end()){
                            //exit(-1);
                            edge_points++;
                            continue;
                        }
                        vector3d vec = it_vec->first;
                        double height = vec.z-v.z;                          //get height
                        double dist=(v-vector3d(x,y,v.z)).length();
                        double phi = atan(height/dist);
                        if(phi>horizon_elevation_angles[theta]){//see if larger
                            horizon_elevation_angles[theta]=phi;
                            dists[theta] = dist/1000;
                        }
                        
                    }
                    interior_points++;
                }
            }
            double cos_theta = cos(th);
            double xend = EQ_DBL(copysign(1,cos_theta),1)?north_x_max:0;     //if the cosinus is positive endvalue is north_x_max, if negative 0
            
            for(double x=v.x+copysign(height_map_resolution,cos_theta);(x >= 0 && x <= north_x_max);x+=copysign(height_map_resolution,cos_theta)){//increase/decrease if theta +/-
                double y = v.y+(x-v.x)*tan(th);                             //find the y that goes with x for this theta
                double y1 = floor(y-((int)y%(int)height_map_resolution));       //find the nearest lower gridpoint by subtracting remainder according to height_map_resolution
                double y2 = y1 + height_map_resolution;                  //Add height_map_resolution to get nearest higher gridpoint
                double phi = horizon_elevation_angles[theta];
                double dist = sqrt((v.x-x)*(v.x-x) + (v.y-y)*(v.y-y));
                double phiMaxTheta = atan(5000.0/dist);      //maximal phi at this theta (with height 5000 m)
                
                if(y1>=north_x_max||y1<0||y2>=north_x_max||y2<0||phi>phiMaxTheta){
                    break;
                }
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
                    double height = vec1.z*(y2-y)/height_map_resolution + vec2.z*(y-y1)/height_map_resolution-v.z;       //compute height at y via linear interpolation
                    double dist=(v-vector3d(x,y,v.z)).length();
                    double phi = atan(height/dist);
                    if(phi>horizon_elevation_angles[theta]){//see if larger
                        horizon_elevation_angles[theta]=phi;
                        dists[theta] = dist/1000;
                    }
                    
                }
                else{//if y is a gridpoint
                    auto it_vec = grid_points.find(vector3d(x,y,0));  //get vector
                    if (it_vec == grid_points.end()){
                        //exit(-1);
                        edge_points++;
                        continue;
                    }
                    vector3d vec = it_vec->first;
                    double height = vec.z-v.z;                          //get height
                    double dist=(v-vector3d(x,y,v.z)).length();
                    double phi = atan(height/dist);
                    if(phi>horizon_elevation_angles[theta]){//see if larger
                        horizon_elevation_angles[theta]=phi;
                        dists[theta] = dist/1000;
                    }
                    
                }
                interior_points++;
            }
            
            
            
            
        }
        
        
        
        
        int k = 0;
        
        std::vector<double> average_sun_intensity;
        average_sun_intensity.assign(times_per_year, 0.0);
        double sun_intensity = 0;
        
        for(int j=0;j<times_per_year;j++){
            for (k=0;k<sun_angles;k++) {
                if(horizon_elevation_angles[k]>sun_elevation_angle[j][k]){
                    sun_intensity=0;
                }
                else{
                    double height = 0;
                    if(elevation_dependant_sun_intensity){
                        height = v.z;
                    }
                    //spherical shell approximation for air mass attenuation, see wikipedia
                    // http://en.wikipedia.org/wiki/Air_mass_%28solar_energy%29
                    double phi = sun_elevation_angle[j][k];
                    double R_E = 6371000.0;      //Earth radius, meters
                    double y_atm = 9000.0;       //Athmospheric thickness, meters
                    double r = R_E/y_atm;
                    double c = height/y_atm;
                    double air_mass_coefficient = sqrt((r+c)*(r+c)*cos(phi)*cos(phi) + (2*r+1+c)*(1-c))-(r+c)*cos(phi);
                    double I_0 = 1.353; // kW/m^2
                    double I = 1.1 * I_0 * pow(0.7, pow(air_mass_coefficient, 0.678));
                    sun_intensity = I * (sun_vec_day_angle[j][k] * Normal);
                    
                    average_sun_intensity[j]+=sun_intensity;
                    
                }
                average_sun_intensity[j] = average_sun_intensity[j] / sun_angles;
                grid_point.second[j]=average_sun_intensity[j];
                max_sun_intensity[j] = (average_sun_intensity[j] > max_sun_intensity[j])?average_sun_intensity[j]:max_sun_intensity[j];
                min_sun_intensity[j] = (average_sun_intensity[j] < min_sun_intensity[j])?average_sun_intensity[j]:min_sun_intensity[j];
            }
            assert(average_sun_intensity.size() == times_per_year);
            grid_points[v]=average_sun_intensity;
            
            assert(grid_point.second.size() == times_per_year);
            
            
        }
    }
        std::cout << grid_points.size() << std::endl;
    
        for(int k = 0; k < times_per_year; ++k){
            std::ostringstream oss("");
            oss << k;
            std::ofstream ofs(output_file + oss.str() + ".xyz");
            if (!ofs.is_open()){
                exit(-2);
            }
            for(auto grid_point : grid_points){
                vector3d v = grid_point.first;
                if(grid_point.second.size()>k){
                    ofs << v.x << " " << v.y << " " << v.z << " ";
                    ofs << grid_points[v][k] << std::endl;
                }
            }
        }
        return 0;
    }
