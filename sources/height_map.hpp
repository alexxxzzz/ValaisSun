//
//  height_map.hpp
//  ComputeHorizons
//
//  Created by alexxx on 25/05/2019.
//

#ifndef height_map_hpp
#define height_map_hpp

#include <unordered_set>
#include "vector3d.hpp"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class height_map {
    std::unordered_set<vector3d, hash> locations;
    double x_max, x_min, y_max, y_min, h_max, h_min;
public:
    height_map (std::ifstream&, double , double , double , double);
    std::unordered_set<vector3d, hash>::const_iterator find_point(vector3d);
    double compute_elevation_angle(vector3d, double, double, double);
    int size(){
        return locations.size();
    }
    double xmax(){
        return x_max;
    }
    double xmin(){
        return x_min;
    }
    double ymax(){
        return y_max;
    }
    double ymin(){
        return y_min;
    }
    double hmax(){
        return h_max;
    }
    double hmin(){
        return h_min;
    }
    bool is_end(std::unordered_set<vector3d, hash>::const_iterator point){
        return (point == locations.end());
    }
};

#endif /* height_map_hpp */