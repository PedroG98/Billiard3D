#include "NumericalIntegrator.h"
#include <cmath>





// PEFRL -------------------------------------------------------------------------------------------------
// Valores constantes para el integrador
constexpr double zeta   = 0.1786178958448091;
constexpr double lambda = -0.2123418310626054;
constexpr double xi     = -0.6626458266981849e-1;

PEFRLIntegrator::PEFRLIntegrator(Vector<double> (*F)(const Vector<double>&, double), Vector<double>& r0, Vector<double>& v0, double mass, double dt)
    : F(F), t(0.0), dt(dt), r(r0), v(v0), m(mass)
{}


// Paso
void PEFRLIntegrator::step()
{
    Vector<double> y1  = r + zeta * v * dt;
    Vector<double> y1p = v + (1 - 2 * lambda) / (2*m) * F(y1, t) * dt;
    Vector<double> y2  = y1 + xi * y1p * dt;
    Vector<double> y2p = y1p + lambda / m * F(y2, t) * dt;
    Vector<double> y3  = y2 + (1 - 2 * (xi + zeta)) * y2p * dt;
    Vector<double> y3p = y2p + lambda / m * F(y3, t) * dt;
    Vector<double> y4  = y3 + xi * y3p * dt;

    v = y3p + (1 - 2 * lambda) / (2*m) * F(y4, t) * dt;
    r = y4 + xi * v * dt;

    t += dt;
}


// Getters
Vector<double> PEFRLIntegrator::getR() const
{
    return r;
}
Vector<double> PEFRLIntegrator::getV() const
{
    return v;
}
double PEFRLIntegrator::getTime() const
{
    return t;
}



/*
std::pair<Vector<double>, Vector<double>> verletStep(Vector<double> (*a)(const Vector<double>&, const double),
                                                     const Vector<double>& xn1, const Vector<double>& xn0, double t, double dt)
{
    Vector<double> xn2 = 2.0*xn1 - xn0 + a(xn1, t) * dt*dt;
    Vector<double> vn1 = (xn2 - xn0) / (2 * dt);
    return std::make_pair(xn2, vn1);
}


std::pair<Vector<double>, Vector<double>> velVerletStep(Vector<double> (*a)(const Vector<double>&, const double),
                                                        const Vector<double>& xn, const Vector<double>& vn, double t, double dt)
{
    Vector<double> xn1 = xn + vn*dt + a(xn, t)*dt*dt/2.0;
    Vector<double> vn1 = vn + dt/2 * (a(xn, t) + a(xn1, t + dt));   // dt/2 (a(n) + a(n+1))
    return std::make_pair(xn1, vn1);
}
*/
