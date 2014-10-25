#ifndef VECTOR3D_HEADER_DEFINED
#define VECTOR3D_HEADER_DEFINED


const double EPS_DBL = 1e-12;

struct vector3d{
    double x;
    double y;
    double z;
    
    vector3d(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0) : x(x_), y(y_), z(z_) {}
    
    double operator *(const vector3d& vec) const{
        return x * vec.x + y * vec.y + z * vec.z;
    }
    vector3d operator / (double d) const {
        return vector3d(x / d, y / d, z / d);
    }
    vector3d norm() const{
        return *this / sqrt(*this * *this);
    }
    double length() const{
        return sqrt(*this * *this);
    }
    vector3d operator ^ (const vector3d& vec) const{
        return vector3d(y*vec.z-z*vec.y, z*vec.x-x*vec.z,x*vec.y-y*vec.x);
    }
    vector3d operator + (const vector3d& vec) const{
        return vector3d(x+vec.x,y+vec.y,z+vec.z);
    }
    vector3d operator - (const vector3d& vec) const{
        return vector3d(x-vec.x,y-vec.y,z-vec.z);
    }
    
    bool operator < (const vector3d& vec) const{
        if(x<vec.x){
            return true;
        }
        else if(((x-vec.x)<EPS_DBL) && (y<vec.y)){
            return true;
        }
        else{
            return false;
        }
    }
    bool operator == (const vector3d& vec) const{
        return (std::abs(x-vec.x)<EPS_DBL) && (std::abs(y-vec.y)<EPS_DBL);
    }
    vector3d interpLin(const vector3d& vec1, const vector3d& vec2){
        double dist1 = (this->x-vec1.x) * (this->x-vec1.x) + (this->y-vec1.y) * (this->y-vec1.y);
        double dist2 = (this->x-vec2.x) * (this->x-vec2.x) + (this->y-vec2.y) * (this->y-vec2.y);
        double h = vec1.z * dist1 / (dist1 + dist2) + vec2.z * dist2 / (dist1 + dist2);
        return vector3d(this->x,this->y,h);
    }
};

std::ostream& operator << (std::ostream& os, const vector3d& vec){
    os << "<" << vec.x << ", " << vec.y << ", " << vec.z << ">";
    return os;
}

struct hash {
    size_t operator() (const vector3d& vec) const {
        return (size_t)vec.y+((size_t)vec.x << 32);
    }
};

#endif // VECTOR3D_HEADER_DEFINED