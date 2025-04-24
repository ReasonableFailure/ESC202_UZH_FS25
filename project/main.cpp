#include <iostream>
#include <cmath>
#include "pedestrian.h"


int main(){
    std::vector<double> ini = {1.0,3.5,2.0};
    std::vector<double> zwei = {0.0,0.0,0.0};
    Pedestrian marco = Pedestrian(zwei);
    marco.setPosition(ini);
    return 0;
}