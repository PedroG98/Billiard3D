#include <iostream>
#include "Billiard.h"
#include <fstream>
#include <cmath>

// TODO
// Optimización: Ainv manual

constexpr double a = 1.0, b = 1.0;
Vector<double> ellipsoid(double s, double t)
{
    using namespace std;
    return Vector<double>{b * cos(t) * sin(s), b * sin(t) * sin(s), a * cos(s)};
}

double ellipsoidDA(double s, double t)
{
    using namespace std;
    return b*sin(s) * sqrt(a*a*sin(s)*sin(s) + b*b*cos(s)*cos(s));
}


int main()
{
    PointMass pointMass;
    pointMass.m = 1;
    pointMass.pos = Vector<double>{0.0, 0.2, 0.5};
    pointMass.vel = Vector<double>{0.5, 0.5, 0.5};

    Surface S(ellipsoid, ellipsoidDA, 10, 0, M_PI, 0, 2*M_PI);

    // I = MR^2 = M
    // Masa y momento de inercia correspondiente

    // Definicion de sistema
    DynamicBilliard system(pointMass, S, {0,0,0}, {-0.05,-0.05,-0.05}, {0,0.5,0}, {0.,0.,0.}, 50, 0.3, 50);

    std::ofstream out("path.dat");
    std::ofstream conservation("conservation.dat");
    std::ofstream pointcare("pointcare.dat");
    int count = 0;
    while (system.getTime() < 1000)
    {
        system.step();
        count++;

        if (count % 100 == 0)
        {
            out << system.getTime() << '\t' << system.getParticlePos() << system.getRotation().flatten() << system.getCM();
            out << std::endl;
            conservation << system.getTime() << '\t' << system.getE() << '\t' << system.getP() << '\t' << system.getL() << std::endl;
        }
        if (count % 5000 == 0)
        {
            std::cout << system.getTime();
            system.flushCutsToFile(pointcare);
        }
    }
    out.close();
    conservation.close();
    system.flushCutsToFile(pointcare);
    pointcare.close();

    return 0;
}
