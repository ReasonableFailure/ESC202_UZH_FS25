#include "pedestrian.h"
#include <vector>
#include <cmath>
#include <iostream>

Pedestrian::Pedestrian(const std::vector<double>& location) : position(location){}

void Pedestrian::setPosition(const std::vector<double>& location){
    this->position = location;
    for(int i = 0; i < location.size(); ++i){
        std::cout << location[i] << std::endl;
        std::cout << this->position[i] << std::endl;
    }

}

