std::pair<double,double> Swiss_To_LatLon(double north, double east){
    east -= 600000.0;                                             // Convert origin to "civil" system, where Bern has coordinates 0,0.
    north -= 200000.0;
    
    east /= 1E6;                                                // Express distances in 1000km units.
    north /= 1E6;
    
    double lon = 2.6779094;                                        // Calculate longitude in 10000" units.
    lon += 4.728982 * east;
    lon += 0.791484 * east * north;
    lon += 0.1306 * east * north * north;
    lon -= 0.0436 * east * east * east;
    
    double lat = 16.9023892;                                       // Calculate latitude in 10000" units.
    lat += 3.238272 * north;
    lat -= 0.270978 * east * east;
    lat -= 0.002528 * north * north;
    lat -= 0.0447 * east * east * north;
    lat -= 0.0140 * north * north * north;
    
    lon *= 100.0 / 36.0;                                            // Convert longitude and latitude back in degrees.
    lat *= 100.0 / 36.0;
    return std::make_pair(lat,lon);
    
}
