
void import_heightmap(std::string input_file, double south_bound, double north_bound, double east_bound, double west_bound, int horizon_angles,
                      std::unordered_map<vector3d, vector3d, hash>& grid_points ,double& x_max, double& x_min, double& y_max, double& y_min){
    
    x_max = -1e100;
    x_min = 1e100;
    y_max = -1e100;
    y_min = 1e100;
    
    std::ifstream ifs(input_file);
    if (!ifs.is_open())
        throw std::runtime_error("could not open file : " + std::string(input_file));
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
            x_min = x<x_min?x:x_min;
            y_min = y<y_min?y:y_min;
            
                                       //y component of vector is direction EAST
                                       //x component of vector is direction NORTH
                                       //THETA ANGLES are from North to West axis

            
            vector3d normal(0,0,0);
            grid_points.insert(std::make_pair(vec, normal));         //Locations is unordered_set of all points in bounding box
            
        }
    }
    

    
}