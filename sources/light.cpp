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
#include <opencv2/opencv.hpp>

#include "vector3d.hpp"
#include "swiss_to_lat_lon.hpp"

//const double Emin=542200;   //Valais
//const double Emax=680200;
//const double Nmin=79600;
//const double Nmax=167800;

const double Emin = 484650;
const double Emax = 843369;    //Suisse
const double Nmin = 74049;
const double Nmax = 299128;

bool EQ_DBL(double a, double b){
    return std::abs(a-b)<EPS_DBL;
}


int main(int ac, char** av) {
    assert(ac >= 2);
    std::cout << "openning file : " << av[1] <<  std::endl;
    std::map<vector3d, double> locations;   //Locations is unordered_set of all points in bounding box
    // std::set<vector3d> locations;
    std::ifstream ifs(av[1]);
    if (!ifs.is_open())
        throw std::runtime_error("could not open file : " + std::string(av[1]));

    double MaxX=0;
    double MaxY=0;
    double MinX=1e12;
    double MinY=1e12;
    
    while (!ifs.eof()) {
        double x;
        double y;
        double h;
        double L;
        ifs >> y;  //import North coordinate in y
        ifs >> x;  //import East coordinate in x
        ifs >> h;  //import Height in h
        ifs >> L;  //import Light in L
        vector3d vec(x, y, h);
        locations.insert(std::make_pair(vec, L));         //Locations is map of all points
        
        MaxX = x>MaxX?x:MaxX;
        MaxY = y>MaxY?y:MaxY;
        MinX = x<MinX?x:MinX;
        MinY = y<MinY?y:MinY;

        
    }
    std::cout << "X: " << MinX <<", " << MaxX << " Y: "<< MinY << " , " << MaxY <<  std::endl;

    double East = Emax - MinX;
    double West = Emax - MaxX;
    double North = Nmin + MaxY;
    double South = Nmin + MinY;
    std::pair<double,double> NE = Swiss_To_LatLon(North, East);
    std::pair<double,double> SW = Swiss_To_LatLon(South, West);

    std::cout << "NE: " << North <<", " << East << " SW: "<< South << " , " << West <<  std::endl;
    std::cout << std::setprecision(9)  << "NE: " << NE.first <<", " << NE.second << " SW: "<< SW.first << " , " << SW.second <<  std::endl;

    double BernNorth = 199937 - Nmin;
    double BernEast = Emax - 599960;
    
    double JetNorth = 118118 - Nmin;
    double JetEast = Emax - 501018;
    
    double ZurichNorth = 246608 - Nmin;
    double ZurichEast = Emax - 683719;
    
    
    std::cout << "Bern: " << BernNorth <<", " << BernEast << std::endl;
    
    //    y = Emax - E ;
//    x = N - Nmin;

    std::cout << "points:" << locations.size() <<  std::endl;
    //cv::Mat img((int)(MaxX*1.5/200), (int)(MaxY*1.5/200), CV_8UC3);
    cv::Mat img((int)((MaxY-MinY)*1.0/200+1),(int)((MaxX-MinX)*1.0/200+1),  CV_8UC4);
    
    std::cout << "imgsize: " << img.cols << " , "<< img.rows << " MaxX-MinX: " << (MaxX-MinX)/200.0 << " MaxY-MinY: " << (MaxY-MinY)/200.0 << std::endl;
    
    //exit(0);
    
    img = cv::Scalar(255,255,255,0);
    
     for(auto pair : locations){
         int x = (0.0*(MaxX-MinX)/200.0+(MaxX-pair.first.x)/200.0);
         int y = (0.0*(MaxY-MinX)/200.0+(MaxY-pair.first.y)/200.0);
         double L = pair.second;
         double hue;
         double  sat;
         double val;
         int alpha;
         int red;
         int green;
         int blue;
         if(L<0){
             hue = 0;
             sat = 0;
             val = 0;
             alpha = 0;
         }
         else{
             hue = (120.0-120.0*L);
             sat = 255.0;
             val = 255.0;
             alpha = 120;
         }
         
         if(abs(pair.first.y-BernNorth)<100&&abs(pair.first.x-BernEast)<100){
             sat = 0;
             val = 0;
             alpha = 255;
             std::cout << "Bern!" << std::endl;
         }
         if(abs(pair.first.y-JetNorth)<100&&abs(pair.first.x-JetEast)<100){
             sat = 0;
             val = 0;
             alpha = 255;
             std::cout << "Jet d'Eau!" << std::endl;
         }
         if(abs(pair.first.y-ZurichNorth)<100&&abs(pair.first.x-ZurichEast)<100){
             sat = 0;
             val = 0;
             alpha = 255;
             std::cout << "Zurich!" << std::endl;
         }

         
         
         double C = (val/255.0)*(sat/255.0);
         double X = C*(1-std::abs((int)(hue/30.0)%2+(hue/30.0)-floor(hue/30.0)-1));
         double m = val/255.0-C;
         C = C + m;
         X = X + m;
         if(hue<30){
             red = C*255.0;
             green = X*255.0;
             blue = m*255.0;
         }
         else if(hue<60){
             red = X*255.0;
             green = C*255.0;
             blue = m*255.0;
         }
         else if(hue<90){
             red = m*255.0;
             green = C*255.0;
             blue = X*255.0;
         }
         else if(hue<120){
             red = m*255.0;
             green = X*255.0;
             blue = C*255.0;
         }
         else if(hue<150){
             red = X*255.0;
             green = m*255.0;
             blue = C*255.0;
         }
         else{
             red = C*255.0;
             green = m*255.0;
             blue = X*255.0;
         }
         
         img.at<cv::Vec4b>(y,x)[0]=blue;
         img.at<cv::Vec4b>(y,x)[1]=green;
         img.at<cv::Vec4b>(y,x)[2]=red;
         img.at<cv::Vec4b>(y,x)[3]=alpha;
//         std::cout << "colors "<<red <<" , "<<green <<", "<<blue <<", "<< m << std::endl;
     }
    
//    cv::cvtColor(img, img, CV_HSV2BGR);
//    cv::namedWindow("xyz");
//    cv::imshow("xyz",img);
//    sleep(100.0);
    
    std::cout << "fileout "<< av[2] << std::endl;
    cv::imwrite(av[2], img);
    return 0;
}
