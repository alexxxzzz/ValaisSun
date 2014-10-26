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

bool EQ_DBL(double a, double b){
    return std::abs(a-b)<EPS_DBL;
}

int main(int ac, char** av) {
    assert(ac == 2);
    std::cout << "openning file : " << av[1] <<  std::endl;
    std::map<vector3d, double> locations;   //Locations is unordered_set of all points in bounding box
    // std::set<vector3d> locations;
    std::ifstream ifs(av[1]);
    if (!ifs.is_open())
        throw std::runtime_error("could not open file : " + std::string(av[1]));

    double MaxX=0;
    double MaxY=0;
    
    while (!ifs.eof()) {
        double x;
        double y;
        double h;
        double L;
        ifs >> x;  //import East coordinate in x
        ifs >> y;  //import North coordinate in y
        ifs >> h;  //import Height in h
        ifs >> L;  //import Light in L
        vector3d vec(x, y, h);
        locations.insert(std::make_pair(vec, L));         //Locations is map of all points
        
        MaxX = x>MaxX?x:MaxX;
        MaxY = y>MaxY?y:MaxY;
        
    }
    std::cout << "points:" << locations.size() <<  std::endl;
    //cv::Mat img((int)(MaxX*1.5/200), (int)(MaxY*1.5/200), CV_8UC3);
    cv::Mat img((int)(MaxX*1.0/200), (int)(MaxY*1.0/200), CV_8UC4);
    
    std::cout << "imgsize: " << img.cols << " , "<< img.rows << std::endl;
    
    
    img = cv::Scalar(255,255,255,0);
    
     for(auto pair : locations){
//         int x = (0.25*MaxX/200.0+(MaxX-pair.first.x)/200.0);
//        int y = (0.25*MaxY/200.0+(MaxY-pair.first.y)/200.0);
         int x = (0.0*MaxX/200.0+(MaxX-pair.first.x)/200.0);
         int y = (0.0*MaxY/200.0+(MaxY-pair.first.y)/200.0);
//         int y =(0.25*MaxY/200.0+pair.first.y/200.0);
         double L = pair.second;
         int hue;
         int sat;
         int val;
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
             hue = (int)(120.0-120.0*L);
             sat = 255;
             val = 255;
             alpha = 90;
         }
         double C = (val/255.0)*(sat/255.0);
         double X = C*(1-std::abs(((int)(hue*2.0/60.0))%2-1));
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
         
         img.at<cv::Vec4b>(x,y)[0]=blue;
         img.at<cv::Vec4b>(x,y)[1]=green;
         img.at<cv::Vec4b>(x,y)[2]=red;
         img.at<cv::Vec4b>(x,y)[3]=alpha;
     }
    
//    cv::cvtColor(img, img, CV_HSV2BGR);
//    cv::namedWindow("xyz");
//    cv::imshow("xyz",img);
//    sleep(100.0);
    
    std::cout << "fileout "<< av[2] << std::endl;
    cv::imwrite(av[2], img);
    return 0;
}
