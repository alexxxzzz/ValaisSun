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
#include "Swiss_To_LatLon.hpp"

namespace po = boost::program_options;




bool EQ_DBL(double a, double b){
        return std::abs(a-b)<EPS_DBL;
}



const double thLat=46.26*M_PI/180.0;


int main(int ac, char** av) {
    //Variables to be assigned by program options
    double Resolution;
    double Nmax;
    double Nmin;
    double Emax;
    double Emin;
    int SunAngles;
    int TimesPerYear;
    int StartDay;
    std::string InputFile;
    std::string OutputFile;

    // Declare the supported options.
    po::options_description op_desc("Allowed options");
    op_desc.add_options()
    ("help", "print options table")
    ("input-file,i", po::value<std::string>(&InputFile), "File containing height map data in x<space>y<space>z<newline> swiss coordinates format")
    ("output-file-base,o", po::value<std::string>(&OutputFile)->default_value("data_out"), "Base filname for output files (default data_out)")
    ("resolution,R", po::value<double>(&Resolution)->default_value(200.0), "resolution of data (default: 200.0)")
    ("nmax",po::value<double>(&Nmax)->default_value(1e100), "maximum north coordinate to be treated (default 1e100)")
    ("nmin",po::value<double>(&Nmin)->default_value(-1e100), "minimum north coordinate to be treated (default -1e100)")
    ("emax",po::value<double>(&Emax)->default_value(1e100), "maximum east coordinate to be treated (default 1e100)")
    ("emin",po::value<double>(&Emin)->default_value(-1e100), "minimum east coordinate to be treated (default -1e100)")
    ("sunangles,s", po::value<int>(&SunAngles)->default_value(360), "number of angles to compute sunlight from (default 360)")
    ("times,t", po::value<int>(&TimesPerYear)->default_value(12), "number of times per year to calculate at (default 12)")
    ("start-day, D", po::value<int>(&StartDay)->default_value(20), "Day of the year to output first data. (default 20)")
    ;
    
    po::positional_options_description pd;
    pd.add("input-file", 1).add("output-file", 1);
    
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, op_desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout <<"CrunchGeoData [options] [input file] [output base]"<<std::endl<< op_desc << "\n";
        return 1;
    }
    
    
    exit(0);
    
    assert(ac == 2);
    std::cout << "openning file : " << InputFile <<  std::endl;
    std::unordered_map<vector3d, std::vector<double>, hash> locations;   //Locations is unordered_set of all points in bounding box
    // std::set<vector3d> locations;
    std::ifstream ifs(av[1]);
    if (!ifs.is_open())
        throw std::runtime_error("could not open file : " + std::string(av[1]));
    
    
    double Xmax=0;
    double Ymax=0;
    double Xmin=1e12;
    double Ymin=1e12;
    
    while (!ifs.eof()) {
        double x;
        double y;
        double h;
        ifs >> y;  //import East coordinate in y
        ifs >> x;  //import North coordinate in x
        ifs >> h;  //import Height in h
        
        
        vector3d vec(x, y, h);
        if(vec.y<Emax && vec.y>Emin && vec.x<Nmax && vec.x>Nmin){ //Check that the imported coordinate is in our bounding rectangle
            Xmax = x>Xmax?x:Xmax;
            Ymax = y>Ymax?y:Ymax;
            Xmin = x<Xmin?x:Xmin;
            Ymin = y<Ymin?y:Ymin;

            vec.y=-vec.y;              //y component of vector is direction WEST
                                       //x component of vecotr is direction NORTH
            
            //THETA ANGLES are from North to West axis
            std::vector<double> t(TimesPerYear);
            //t.assign(TimesPerYear,3);
            locations.insert(std::make_pair(vec, t));         //Locations is unordered_set of all points in bounding box
            
        }
    }
    
    std::cout << "NE: " << Xmax <<", " << Ymax << " SW: "<< Xmin << " , " << Ymin <<  std::endl;

    std::pair<double,double> NE = Swiss_To_LatLon(Xmax+Resolution/2.0, Ymax+Resolution/2.0);
    std::pair<double,double> SW = Swiss_To_LatLon(Xmin-Resolution/2.0, Ymin-Resolution/2.0);

    std::cout << std::setprecision(9) << "NE: " << NE.first <<", " << NE.second << " SW: "<< SW.first << ", " << SW.second <<  std::endl;
    
    
    int Npoints = locations.size();
    std::cout << "Number of points: " << Npoints << std::endl;
    
    std::vector<std::vector<double> > phiAxMonthAngle(TimesPerYear,std::vector<double>(SunAngles));  //Elevation angle PHI of Earth's axis for a month(22nd) and hour
    std::vector<std::vector<double> > sunIntensityMonthAngle(TimesPerYear,std::vector<double>(SunAngles));
    std::vector<double> pts_theta(SunAngles);   //vector used for plotting the suns elevation, x-vector, hours
    std::vector<double> pts_thElv_y(SunAngles);  //vector used for plotting the suns elevation, y-vector, elevation angle
    std::vector<double> pts_dx(SunAngles);
    std::vector<double> pts_dy(SunAngles);
    
    std::vector<std::vector <vector3d> > Sun_Vec_Month_Angle(TimesPerYear,std::vector<vector3d>(SunAngles));  //Unit vectors for sun's direction in month(22nd), hour.
    
    //    int Plot_Month=5;  //Month we will plot the sun elevation for
    
    for(int Month = 0; Month < TimesPerYear; Month++){
        double N=(365.0*Month)/TimesPerYear+StartDay;  //We are computing values for the Nth day of the year
        double phiAxMonth=-asin(0.39779*cos((0.98565*(N+10)+1.914*sin(0.98565*(N-2)*M_PI/180))*M_PI/180));  //Earth axis inclination for day N
        for(int thH=0;thH<SunAngles;thH++){//Angle theta for sun; 0 is North, +90 is West.
            int indTH= thH;
            double thHrad = (180-thH)*M_PI/180.0;  //turn around for sun formula angle
            double phiElv=asin(sin(thLat)*sin(phiAxMonth)+cos(thLat)*cos(thHrad)*cos(phiAxMonth));  //Sun elevation, imperical fomula, see wikipedia
            
            thHrad = M_PI-thHrad;          //In our coordinate system (North x+, West y+) we must inverse the sun theta. Theta corresponds to sun ray direction.
            
            vector3d sun_vec(cos(phiElv) * cos(thHrad) , cos(phiElv) * sin(thHrad) , sin(phiElv));  //Unit vector for direction sun ray come from
            //std::cout << "Month: " << (int)Month<< "indH: "<< indH << std::endl;
            //            phiAxMonthAngle[(int)Month][indH]=thElv;
            phiAxMonthAngle[(int)Month][indTH]=phiElv;  //stick sun elevation in its datastructure
            
            Sun_Vec_Month_Angle[(int)Month][indTH]=sun_vec;  //stick sun vector in its datastructure
            
            double air_mass_coefficient =sqrt(708.0*708.0*sin(phiElv)*sin(phiElv)+2.0*708.0+1.0)-708.0*sin(phiElv);  //wikipedia
            sunIntensityMonthAngle[(int)Month][indTH]=pow(0.7,pow(air_mass_coefficient,0.678))/0.7;
            
            /*           if(Month==Plot_Month){//for plotting
             pts_theta[indTH]=thH;
             pts_thElv_y[indTH]=phiElv*180.0/M_PI;
             pts_dx[indTH] = cos(thHrad)*10.0;
             pts_dy[indTH] = sin(thHrad)*5.0;
             //std::cout <<"Hour: "<< hour <<" dx: "<< pts_dx[indH]<<" dy: "<< pts_dy[indH]<< std::endl;
             }*/
        }
    }
    
    
    
    
   // vector3d v(Ntest, Wtest,0);
   // auto it_v = locations.find(v);          //get our test vector
    
  //  if (it_v == locations.end()){
  //      std::cout <<"Testpoint not found"<< std::endl;
  //      exit(-1);
  //  }
    std::vector<double> Max_Sun;
    Max_Sun.assign(TimesPerYear, 0);
    std::vector<double> Min_Sun;
    Min_Sun.assign(TimesPerYear, 1e12);
    
    
    int N=0;
    
    std::cout <<"MAX: "<<Max_Sun[5]<<"   MIN: "<<Min_Sun[5]<<"   Progress: "<< 100*(N*1.0)/(Npoints*1.0) <<" %    \n";

    for(auto pair : locations){
        vector3d v = pair.first;
//        std::vector<double>& vL=pair.second;
        N++;
        
        if (1) {
            std::cout <<"MAX: "<<Max_Sun[5]<<"   MIN: "<<Min_Sun[5]<<"   Progress: "<< 100.0*(N*1.0)/(Npoints*1.0) <<" %  ("<<N<<")    \r";
            std::cout.flush();
        }
        vector3d vNx(v.x+Resolution*10,v.y,0);            //get points 10xresolution to the north, west, south and east
        vector3d vWx(v.x,v.y+Resolution*10,0);
        vector3d vSx(v.x-Resolution*10,v.y,0);
        vector3d vEx(v.x,v.y-Resolution*10,0);
        auto itEx=locations.find(vEx);
        auto itNx=locations.find(vNx);
        auto itWx=locations.find(vWx);
        auto itSx=locations.find(vSx);
        if (itEx == locations.end()||itNx == locations.end()||itWx == locations.end()||itSx == locations.end()){
            std::vector<double> L;
            L.assign(TimesPerYear,-1.0);
            locations[v]=L;
            continue;
        }
        
        //v = *it_v;
        
        vector3d vN(v.x+Resolution,v.y,0);            //get points to the north, west, south and east
        vector3d vW(v.x,v.y+Resolution,0);
        vector3d vS(v.x-Resolution,v.y,0);
        vector3d vE(v.x,v.y-Resolution,0);
        
        auto itE=locations.find(vE);
        auto itN=locations.find(vN);
        auto itW=locations.find(vW);
        auto itS=locations.find(vS);
        
        if (itE == locations.end()||itN == locations.end()||itW == locations.end()||itS == locations.end()){
            exit(-1);
        }
        
        vE=itE->first;
        vN=itN->first;
        vW=itW->first;
        vS=itS->first;
        
        vector3d Normal = ((((vE-v)^(vN-v))+((vW-v)^(vS-v)))/2).norm();       //normal is the crossproduct of two prependicular differences. Avgd.
        
        /*    std::vector<double> sun_intensity_with_normal(SunAngles);
         for(int k=0;k<SunAngles;k++){
         sun_intensity_with_normal[k]=std::abs(Sun_Vec_Month_Angle[Plot_Month][k]*Normal)*sunIntensityMonthAngle[Plot_Month][k];
         }*/
        
        /*
         Gnuplot gp;    //Use gnuplot
         gp << "set xrange [0:360]\nset yrange [-1:1]\n";  //set range for plotting
         gp << "plot '-' with lines, '-' with lines\n";                       //plot x,y with lines connecting them
         gp.send1d(boost::make_tuple(pts_theta,sunIntensityMonthAngle[Plot_Month]));
         gp.send1d(boost::make_tuple(pts_theta,sun_intensity_with_normal));
         */
        //    vector3d Normal1= ((vE-v)^(vN-v)).norm();
        //    vector3d Normal2= ((vW-v)^(vS-v)).norm();
        
        //   std::cout << "Norm: " << Normal << " Norm1: " << Normal1 << " Norm2: " << Normal2 << std::endl;
        
        std::vector<double> horizon_elevation_angles;   //vector to hold elevation angles at theta
        horizon_elevation_angles.assign(SunAngles, 0);
        //    std::vector<vector3d> horizon(SunAngles);
        std::vector<double> dists(SunAngles);
        
        size_t edge_points = 0;
        size_t interior_points = 0;
        
        for(int theta = 0; theta<SunAngles;theta++){        //compute for each angle theta
            double th = theta*M_PI/180;
            {
                double sin_theta = sin(th);   //calculate sin of the angle in radians
                double yend = EQ_DBL(copysign(1,sin_theta),1)?Ymax:0;     //if the sinus is positive endvalue is Ymax, if negative 0
                for(double y=v.y+copysign(Resolution,sin_theta);(y >= 0 && y <= Ymax);y+=copysign(Resolution,sin_theta)){//increase/decrease if theta +/-
                    double x = v.x+(y-v.y)/tan(th);                             //find the x that goes with y for this theta
                    double x1 = floor(x-((int)x%(int)Resolution));       //find the nearest lower gridpoint by subtracting remainder according to Resolution
                    double x2 = x1 + Resolution;                  //Add Resolution to get nearest higher gridpoint
                    double phi = horizon_elevation_angles[theta];
                    double dist = (v-vector3d(x,y,v.z)).length();
                    double phiMaxTheta = atan(5000.0/dist);      //maximal phi at this theta (with height 5000 m)
                    
                    
                    if(x1>=Xmax||x1<0||x2>=Xmax||x2<0||phi>phiMaxTheta){
                        break;
                    }
                    if((int)x%(int)Resolution){                       //if x is not a grid point
                        
                        //std::cout << vector3d(x1,y,0) << ", " << vector3d(x2,y,0) << " , " << x << std::endl;
                        auto it_vec1 = locations.find(vector3d(x1,y,0));      //get the two gridpoints from the set
                        auto it_vec2 = locations.find(vector3d(x2,y,0));
                        if (it_vec1 == locations.end()||it_vec2 == locations.end()){
                            //exit(-1);
                            
                            // std::cout <<(it_vec1 == locations.end())<<(it_vec2 == locations.end())<<"theta: "<<theta<<" y: "<<y<<", x1: "<<x1<<", x2: "<<x2<<" oops Xmax: "<<Xmax<< ", Ymax: " << Ymax << std::endl;
                            edge_points++;
                            continue;
                        }
                        vector3d vec1 = it_vec1->first;
                        vector3d vec2 = it_vec2->first;
                        double height = vec1.z*(x2-x)/Resolution + vec2.z*(x-x1)/Resolution-v.z;       //compute height at x via linear interpolation
                        dist=(v-vector3d(x,y,v.z)).length();
                        phi = atan(height/dist);
                        //std::cout <<height<< ", " << dist << " , " << phi << std::endl;
                        if(phi>horizon_elevation_angles[theta]){//see if larger
                            horizon_elevation_angles[theta]=phi;
                            dists[theta] = dist/1000;
                        }
                        
                    }
                    else{//if x is a gridpoint
                        auto it_vec = locations.find(vector3d(x,y,0));  //get vector
                        if (it_vec == locations.end()){
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
            double xend = EQ_DBL(copysign(1,cos_theta),1)?Xmax:0;     //if the cosinus is positive endvalue is Xmax, if negative 0
            
            for(double x=v.x+copysign(Resolution,cos_theta);(x >= 0 && x <= Xmax);x+=copysign(Resolution,cos_theta)){//increase/decrease if theta +/-
                double y = v.y+(x-v.x)*tan(th);                             //find the y that goes with x for this theta
                double y1 = floor(y-((int)y%(int)Resolution));       //find the nearest lower gridpoint by subtracting remainder according to Resolution
                double y2 = y1 + Resolution;                  //Add Resolution to get nearest higher gridpoint
                double phi = horizon_elevation_angles[theta];
                double dist = sqrt((v.x-x)*(v.x-x) + (v.y-y)*(v.y-y));
                double phiMaxTheta = atan(5000.0/dist);      //maximal phi at this theta (with height 5000 m)
                
                if(y1>=Xmax||y1<0||y2>=Xmax||y2<0||phi>phiMaxTheta){
                    break;
                }
                if((int)y%(int)Resolution){                       //if y is not a grid point
                    auto it_vec1 = locations.find(vector3d(x,y1,0));      //get the two gridpoints from the set
                    auto it_vec2 = locations.find(vector3d(x,y2,0));
                    if (it_vec1 == locations.end()||it_vec2 == locations.end()){
                        //exit(-1);
                        edge_points++;
                        continue;
                    }
                    vector3d vec1 = it_vec1->first;
                    vector3d vec2 = it_vec2->first;
                    double height = vec1.z*(y2-y)/Resolution + vec2.z*(y-y1)/Resolution-v.z;       //compute height at y via linear interpolation
                    double dist=(v-vector3d(x,y,v.z)).length();
                    double phi = atan(height/dist);
                    if(phi>horizon_elevation_angles[theta]){//see if larger
                        horizon_elevation_angles[theta]=phi;
                        dists[theta] = dist/1000;
                    }
                    
                }
                else{//if y is a gridpoint
                    auto it_vec = locations.find(vector3d(x,y,0));  //get vector
                    if (it_vec == locations.end()){
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
        //    std::vector<double> thetas(SunAngles);
        //    std::vector<double> phisD(SunAngles);
        //    std::vector<double> Sun_Intensity_All_Inclusive(SunAngles);
        
        
        std::vector<double> total_sun;
        total_sun.assign(TimesPerYear, 0.0);
        double sun_intensity = 0;
        
        for(int j=0;j<TimesPerYear;j++){
            for (k=0;k<SunAngles;k++) {
                //           thetas[k]=k;
                //           phisD[k]=horizon_elevation_angles[k]*180/M_PI;
                if(horizon_elevation_angles[k]>phiAxMonthAngle[j][k]){
                    sun_intensity=0;
                }
                else{
                    sun_intensity= std::abs(Sun_Vec_Month_Angle[j][k]*Normal)*sunIntensityMonthAngle[j][k];
                }
                total_sun[j]+=sun_intensity;
                
            }
            pair.second[j]=total_sun[j];
            Max_Sun[j] = (total_sun[j] > Max_Sun[j])?total_sun[j]:Max_Sun[j];
            Min_Sun[j] = (total_sun[j] < Min_Sun[j])?total_sun[j]:Min_Sun[j];
            //std::cout <<"total_sun["<<j <<"]: "<<total_sun[j]<<" pair.second["<<j <<"]: "<<pair.second[j]<<std::endl;
        }
        assert(total_sun.size() == TimesPerYear);
        locations[v]=total_sun;
//        std::cout <<"total_sun[0]: "<<total_sun[0]<<" locations[v][0]: "<<locations[v][0]<<std::endl;
        //vL = total_sun;
        assert(pair.second.size() == TimesPerYear);
        
        
        //std::cout <<"interior_points: "<<interior_points <<" edge_points: "<<edge_points<<std::endl;
        
        
        /*    Gnuplot gp2;
         gp2 << "set xrange [0:360]\nset yrange [-60:90]\n";
         gp2 << "plot '-' with lines, '-' with vectors\n";
         gp2.send1d(boost::make_tuple(thetas, phisD));
         gp2.send1d(boost::make_tuple(pts_theta, pts_thElv_y, pts_dx, pts_dy));
         
         Gnuplot gp3;
         gp3 << "set xrange [0:360]\nset yrange [-1:1]\n";
         gp3 << "plot '-' with lines\n";
         gp3.send1d(boost::make_tuple(thetas, Sun_Intensity_All_Inclusive));
         */
        
        
        //    gp2.send1d(boost::make_tuple(thetas, dists));
        
        //    gp2.send1d(boost::make_tuple(pts_Hour_theta, pts_thElv_y));
    }
    std::cout << locations.size() << std::endl;
    std::cout <<"Max_Sun: "<<Max_Sun[5]<<"Min_Sun: "<<Min_Sun[5]<<std::endl;
    
    for(int k = 0; k < TimesPerYear; ++k){
        std::ostringstream oss("");
        oss << k;
        std::ofstream ofs(OutputFile + oss.str() + ".xyz");
        if (!ofs.is_open()){
            exit(-2);
        }
        for(auto pair : locations){
            vector3d v = pair.first;
//            std::vector<double>& vL=pair.second;
//            std::cout <<"vL.size(): "<<vL.size()<<"pair.second.size(): "<<pair.second.size()<<std::endl;
            if(pair.second.size()>k){
                ofs << v.x << " " << v.y << " " << v.z << " ";
                ofs << (locations[v][k]-Min_Sun[k])/(Max_Sun[k]-Min_Sun[k]) << std::endl;
//                ofs << locations[v][k] << std::endl;
            }
        }
    }
    return 0;
}
