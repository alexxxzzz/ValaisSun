//
//  vectorHM.hpp
//  ComputeHorizons
//
//  Created by alexxx on 26/05/2019.
//

#ifndef vectorHM_hpp
#define vectorHM_hpp

struct vectorHM{
    int x;
    int y;
    float z;
    
    vectorHM(int x_ = 0, int y_ = 
             0, float z_ = 0.0) : x(x_), y(y_), z(z_) {}
    
    double distxy(const vectorHM& vec) const{
        return sqrt((double)(x-vec.x)*(double)(x-vec.x)+(double)(y-vec.y)*(double)(y-vec.y));
    }
    vector3d operator + (const vector3d& vec) const{
        return vector3d(x+vec.x,y+vec.y,z+vec.z);
    }
    vector3d operator - (const vector3d& vec) const{
        return vector3d(x-vec.x,y-vec.y,z-vec.z);
    }
    
    bool operator < (const vectorHM& vec) const{
        if(x<vec.x){
            return true;
        }
        else if((x == vec.x) && (y<vec.y)){
            return true;
        }
        else{
            return false;
        }
    }
    bool operator == (const vector3d& vec) const{
        return ((x == vec.x) && (y == vec.y));
    }
};

inline std::ostream& operator << (std::ostream& os, const vector3d& vec){
    os << "<" << vec.x << ", " << vec.y << ", " << vec.z << ">";
    return os;
}

struct hash {
    size_t operator() (const vector3d& vec) const {
        return std::hash<int>() ((int)(vec.y/25))+(((int)(vec.x/25)) << 16);
    }
};

#endif /* vectorHM_hpp */
