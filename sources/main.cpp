#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>

struct position {
    double x;
    double y;
    double h;
};

const double Emin=542277.0;
const double Emax=680234.0;
const double Nmin=79595.0;
const double Nmax=167709.0;

double Etest=596589.0;
double Ntest=112180.0;

int main(int ac, char** av) {
    std::vector<position> map;
    std::ifstream ifs(av[1]);
    while (!ifs.eof()) {
        position p;
        ifs >> p.x;
        ifs >> p.y;
        ifs >> p.h;
        if(p.x<Emax&&p.x>Emin&&p.y<Nmax&&p.y>Nmin){
            map.push_back(p);
        }
        std::vector <double> phis(360);
    }
//    for (auto it : map) {
//        std::cout << it.x << ", " << it.y << " : " << it.h << std::endl;
//    }
    
    
    std::cout << map.size() << std::endl;
    return 0;
}