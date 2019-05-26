//
//  height_map.cpp
//  ComputeHorizons
//
//  Created by alexxx on 25/05/2019.
//

#include <iostream>
#include <fstream>
#include "height_map.hpp"



height_map::height_map (std::ifstream& ifs, double south_bound, double north_bound, double east_bound, double west_bound){
    
    x_max = -1e100;
    x_min = 1e100;
    y_max = -1e100;
    y_min = 1e100;
    h_max = -1e100;
    h_min = 1e100;
    
    while (!ifs.eof()) {
        double x;
        double y;
        double h;
        ifs >> y;  //import East coordinate in y
        ifs >> x;  //import North coordinate in x
        ifs >> h;  //import Height in h
        
        vector3d vec(x, y, h);
        if(vec.y<east_bound && vec.y>west_bound && vec.x<north_bound && vec.x>south_bound){ //Check that the imported coordinate is in our bounding rectangle
            x_max = x>x_max?x:x_max;
            y_max = y>y_max?y:y_max;
            h_max = h>h_max?h:h_max;
            x_min = x<x_min?x:x_min;
            y_min = y<y_min?y:y_min;
            h_min = h>h_min?h:h_min;
            //y component of vector is increasing direction EAST
            //x component of vector is increasing direction NORTH
            
            vector3d normal(0,0,0);
            locations.insert(vec);         //Locations is unordered_set of all points in bounding box
            
        }
    }
}

std::unordered_set<vector3d, hash>::const_iterator height_map::find_point(vector3d point){
    return locations.find(point);
}

double height_map::compute_elevation_angle(vector3d point, double theta, double height_map_resolution, double angular_resolution){
    double resolution_distance = height_map_resolution / tan(angular_resolution*M_PI/180.0);
    double sin_theta = sin(theta*M_PI/180.0);
    double cos_theta = cos(theta*M_PI/180.0);
    //double tan_theta = tan(theta*M_PI/180.0);
    int x_step = (abs(cos_theta)<=EPS_DBL?0:round(cos_theta/abs(cos_theta)))*height_map_resolution;
    //x_step is +/-height_map_resolution according to sign of cos_theta. Or zero when cos_theta < EPS_DBL
    int y_step = (abs(sin_theta)<=EPS_DBL?0:round(sin_theta/abs(sin_theta)))*height_map_resolution;
    sin_theta = abs(sin_theta) < EPS_DBL ? EPS_DBL : sin_theta;  //sin_theta should not be zero
    cos_theta = abs(cos_theta) < EPS_DBL ? EPS_DBL : cos_theta;  //cos_theta should not be zero
    double tan_theta = sin_theta / cos_theta;
    double cot_theta = cos_theta / sin_theta;
    auto found_point = locations.find(point); //find the point that was passed to us. necessary because
    //caller might have only passed x and y coord, not z
    if(found_point == locations.end()){
        throw std::runtime_error("point not found"); //error if not found
    }
    point = *found_point;
    
    double phi = -acos(6378000/(6378000+point.z));   //start our output elevation angle at maximum negative
    
    double x_coord = x_step;
    double y_coord = y_step;
    double max_dist = 1e50;
    double dist = 0;
    double phi2;
    while(dist<max_dist){
        double res_mult = (floor(dist/resolution_distance) + 1);
     //   std::cout << "theta = " << theta << " phi = " << phi << " dist = "<< dist << " max_dist = " << max_dist << std::endl;
        if(abs(x_coord/cos_theta-y_coord/sin_theta)<EPS_DBL || x_coord/cos_theta < EPS_DBL ||y_coord/sin_theta < EPS_DBL){
            auto v_itr = locations.find(vector3d(point.x + x_coord, point.y + y_coord,0));
            if(v_itr == locations.end())
                break;
            vector3d v = *v_itr;
            dist = point.distxy(v);
            phi2 = asin((v.z-point.z)/(dist));
            x_coord = x_coord + x_step*res_mult;
            y_coord = y_coord + y_step*res_mult;
            
        }
        
        else if(x_coord/cos_theta < y_coord/sin_theta){
            double y_coord_x = x_coord*tan_theta;
            double y_1 = trunc(y_coord_x/height_map_resolution)*height_map_resolution;
            auto v1_itr = locations.find(vector3d(point.x + x_coord, point.y + y_1,0));
            if(v1_itr == locations.end())
                break;
            auto v2_itr = locations.find(vector3d(point.x + x_coord, point.y + y_1 + y_step,0));
            if(v2_itr == locations.end())
                break;
            vector3d v1 = *v1_itr;
            vector3d v2 = *v2_itr;
            double height_difference = (abs(y_coord_x-y_1)/height_map_resolution)*v2.z+(abs(y_coord_x-(y_1+y_step))/height_map_resolution)*v1.z-point.z;
            dist = point.distxy(vector3d(point.x + x_coord, point.y + y_coord_x,0));
            phi2 = asin(height_difference/dist);
            x_coord = x_coord + x_step*res_mult;
        }
        else if(x_coord/cos_theta > y_coord/sin_theta){
            double x_coord_y = y_coord*cot_theta;
            double x_1 = trunc(x_coord_y/height_map_resolution)*height_map_resolution;
            auto v1_itr = locations.find(vector3d(point.x + x_1, point.y + y_coord,0));
            if(v1_itr == locations.end())
                break;
            auto v2_itr = locations.find(vector3d(point.x + x_1 + x_step, point.y + y_coord,0));
            if(v2_itr == locations.end())
                break;
            vector3d v1 = *v1_itr;
            vector3d v2 = *v2_itr;
            double height_difference = (abs(x_coord_y-x_1)/height_map_resolution)*v2.z+(abs(x_coord_y-(x_1+x_step))/height_map_resolution)*v1.z-point.z;
            dist = point.distxy(vector3d(point.x + x_coord_y, point.y + y_coord,0));
            phi2 = asin(height_difference/dist);
            y_coord = y_coord + y_step*res_mult;
        }
        if(phi2 > phi){
            phi = phi2;
            if(phi > EPS_DBL)
                max_dist = (h_max-point.z)/tan(phi);
        }

    }
    return phi;
}


