#include <iostream>
#include <vector>

class Pedestrian{
    private:
        std::vector<double> position;
        std::vector<double> velocity;
        std::vector<double> acceleration;
    protected:

    public:
        Pedestrian(const std::vector<double> & location);
        // ~Pedestrian();
        void setPosition(const std::vector<double>& location);
        // void setVelocity(const std::vector<double>& speed);
        // void setAcceleration(const std::vector<double>& acc);
        // double getAbsoluteSpeed();
      

};