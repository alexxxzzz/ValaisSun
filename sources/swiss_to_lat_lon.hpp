#pragma once

std::pair<double,double> swiss_to_lat_lon(double north, double east){
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

std::pair<double,double> lat_lon_to_swiss(double lat, double lon){
    lat *= 3600;                                                // Convert latitude and longitude in seconds.
    lon *= 3600;
    
    lat -= 169028.66;                                           // Shift the origin in Bern.
    lon -= 26782.5;
    
    lat /= 10000;                                               // Convert latitude and longitude in 10000" units.
    lon /= 10000;
    
    double east = 600072.37;                                       // Calculate easting [m].
    east += 211455.93 * lon;
    east -= 10938.51 * lon * lat;
    east -= 0.36 * lon * lat * lat;
    east -= 44.54 * lon * lon * lon;
    
    double north = 200147.07;                                      // Calculate northing [m].
    north += 308807.95 * lat;
    north += 3745.25 * lon * lon;
    north += 76.63 * lat * lat;
    north -= 194.56 * lon * lon * lat;
    north += 119.79 * lat * lat * lat;
    
    
    return std::make_pair(north,east);
}
