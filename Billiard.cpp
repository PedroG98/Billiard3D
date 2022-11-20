#include "Billiard.h"
#include <cmath>

#include <algorithm>



// Masa puntual y superficie
double PointMass::getE() const
{
    return m * dot(vel,vel) / 2;
}


Surface::Surface(Vector<double> (*gamma)(double s, double t), double (*da)(double s, double t), double mass, double si, double sf, double ti, double tf) :
    r(gamma), da(da), mass(mass), si(si), sf(sf), ti(ti), tf(tf),
    I(momentOfInertia()), Iinv(I.inverse()) // RAII
{}


Vector<double> Surface::operator()(double s, double t) const
{
    return r(s, t);
}

double Surface::dA(double s, double t) const
{
    return da(s, t);
}

Matrix<double> Surface::momentOfInertia() const
{
    // Alta precisión pues es sólo 1 vez
    double Ix = 0.0, Iy = 0.0, Iz = 0.0;
    double surfaceArea = 0.0;        // Para comparar con masa
    double ds = (sf - si) / 128;
    double dt = (tf - ti) / 128;
    for (unsigned int i = 0; i < 128; i++)
        for (unsigned int j = 0; j < 128; j++)
        {
            // Aux
            double s = si + i*ds;
            double t = ti + j*dt;
            Vector<double> r_ = r(s, t);
            double dA_ = dA(s, t);

            Ix += (r_[1]*r_[1] + r_[2]*r_[2]) * dA_ * ds*dt;
            Iy += (r_[0]*r_[0] + r_[2]*r_[2]) * dA_ * ds*dt;
            Iz += (r_[1]*r_[1] + r_[0]*r_[0]) * dA_ * ds*dt;
            surfaceArea += dA_ * ds * dt;
        }
    Matrix<double> ret = Matrix<double>(3, 3, {Ix, 0, 0,
                                               0, Iy, 0,
                                               0, 0, Iz});
    return ret * mass / surfaceArea;
}


// Aux

Matrix<double> rot(unsigned int axis, double angle)
{
    using namespace std;    //cmath

    Matrix<double> I = Identity<double>(3);
    // Generar ángulos respectivos
    if (axis == 0)  // x
    {
        I(1,1) = cos(angle);    I(1,2) = -sin(angle);
        I(2,1) = sin(angle);    I(2,2) = cos(angle);
    }
    else if (axis == 1)  // y
    {
        throw "Not implemented";
    }
    else if (axis == 2)  // z
    {
        I(0,0) = cos(angle);    I(0,1) = -sin(angle);
        I(1,0) = sin(angle);    I(1,1) = cos(angle);
    }
    else
        throw std::invalid_argument("Ejes de rotacion desde 0 a 2");

    return I;
}
// Matriz de rotación 3D derivada. dlambda/dtheta
/*
Matrix<double> Drot(unsigned int axis, double angle)
{
    using namespace std;    //cmath

    Matrix<double> I(3,3, 0.0);
    // Generar ángulos respectivos
    if (axis == 0)  // x
    {
        I(1,1) = -sin(angle);    I(1,2) = -cos(angle);
        I(2,1) = cos(angle);   I(2,2) = -sin(angle);
    }
    else if (axis == 1)  // y
    {
        throw "Not implemented";
    }
    else if (axis == 2)  // z
    {
        I(0,0) = -sin(angle);    I(0,1) = -cos(angle);
        I(1,0) = cos(angle);   I(1,1) = -sin(angle);
    }
    else
        throw std::invalid_argument("Ejes de rotacion desde 0 a 2");

    return I;
}*/



// Billiard -------------------------------------------------------------------------------------------------
// Constructor
Billiard::Billiard(const PointMass& particle, const Surface& surface,
                   const Vector<double>& R0, const Vector<double>& V0, const Vector<double>& theta0, const Vector<double>& thetap0,
                   double height, double elasticity, unsigned int mesh_partition) :
    S(surface), P(particle), pos(R0), vel(V0), theta(theta0), thetap(thetap0), U(height), alpha(elasticity),
    integrator(Rk4Integrator<Vector<double> (Billiard::*)(const Vector<double>&,double), Billiard>())
{
    Rt = lambdaT(theta0);
    R = Rt.transposed();
    calculateOmega();

    // Generar mesh de detección
    partitionSizeT = (S.tf-S.ti) / mesh_partition;
    partitionSizeS = (S.sf-S.si) / mesh_partition;
    for (unsigned int t = 0; t < mesh_partition; t++)
        for (unsigned int s = 0; s < mesh_partition; s++)
        {
            double t_ = S.ti + partitionSizeT * t;
            double s_ = S.si + partitionSizeS * s;
            detectionMesh.push_back(std::make_pair(t_,s_));
        }

    // Integrador rk4
    Vector<double> q0 = {P.pos, pos, theta, P.vel, vel, thetap};
    integrator = Rk4Integrator(&Billiard::F, *this, q0, 0.01, true, 1e-3, 0.01);
}


// Revisar apuntes para derivación. xi(r) = phi'(r)/r
double Billiard::Xi(double r) const
{
    using namespace std;
    if (r > alpha)
        return 0;
    /* Cuadrática
    else
        return 2*U/alpha * (1/alpha - 1/r);//-U / (alpha * r);
    */
    // Bump
    else
    {
        double xma2 = r*r - alpha*alpha;
        return -2*U*M_E*alpha*alpha * exp(alpha*alpha / xma2) / (xma2*xma2);
    }
}


// r = P.r - Rcm - lambda gamma
Vector<double> Billiard::gamma(double s, double t, const Vector<double>& Rcm) const
{
    // S está en sistema S', Rt para pasar a sistema lab S
    return Rcm + Rt * S(s,t);
}


// Matriz de rotación de ángulos de Euler
// Covariante: Para los vectores unitarios. La inversa para los componentes
Matrix<double> Billiard::lambdaT(const Vector<double>& theta_) const
{
    // Consiste en la triple multiplicación matricial de rotaciones en z,x,z
    // Testeado en Wolfram
    using namespace std;    //cmath


    // Esto se escribe explícitamente abajo
    // return rot(2, theta_[2]) * rot(0, theta_[1]) * rot(2, theta_[0]);
    // Por mathematica
    double c1 = cos(theta_[0]), c2 = cos(theta_[1]), c3 = cos(theta_[2]);
    double s1 = sin(theta_[0]), s2 = sin(theta_[1]), s3 = sin(theta_[2]);
    return Matrix<double>(3, 3, {c1*c3 - c2*s1*s3 , -c3*s1 - c1*c2*s3, s2*s3,
                                 c2*c3*s1 + c1*s3 , c1*c2*c3 - s1*s3 , -c3*s2,
                                 s1*s2            , c1*s2            , c2});
}




// Modifica miembro integrationMesh
void Billiard::generateIntegrationMesh(const Vector<double>& r_,
                                       const Vector<double>& R_)
{
    // Itera sobre los puntos de detección y encuentra los que hacen contacto con la partícula
    std::vector<std::pair<double,double>> in_contact;
    for (const auto& detector : detectionMesh)
    {
        double d = (r_ - gamma(detector.second, detector.first, R_)).norm();
        if (d <= alpha)
            in_contact.push_back(detector);
    }

    // Genera mesh sobre cuadrado en torno a cada punto detectado, sin overlap, 10x10 puntos
    integrationMesh.clear();
    for (const auto& detected : in_contact)
        for (unsigned int i = 0; i < 10; i++)
            for (unsigned int j = 0; j < 10; j++)
            {
                double t = detected.first + partitionSizeT * (-.5 + i/10.);
                double s = detected.second + partitionSizeS * (-.5 + j/10.);
                integrationMesh.push_back(std::make_pair(t, s));
            }
}



// Integración rectangular simple en 2 variables
Vector<double> Billiard::calculateSxigamma(const Vector<double>& r_, const Vector<double>& R_) const
{
    // Tamaño de partición
    double hs = partitionSizeS / 10;
    double ht = partitionSizeT / 10;

    Vector<double> ret = Vector<double>(3);
    for (const auto& point : integrationMesh)
    {
        double s = point.second, t = point.first;
        Vector<double> gamma_ = gamma(s, t, R_);
        double dA = S.dA(s, t);
        ret += dA * hs*ht * Xi((r_ - gamma_).norm()) * gamma_;
    }

    return ret;
}


// Idem, sólo xi
double Billiard::calculateSxi(const Vector<double>& r_, const Vector<double>& R_) const
{
    double hs = partitionSizeS / 10;
    double ht = partitionSizeT / 10;

    double ret = 0.0;
    for (const auto& point : integrationMesh)
    {
        double s = point.second, t = point.first;
        Vector<double> gamma_ = gamma(s, t, R_);
        double dA = S.dA(s, t);
        ret += dA * hs*ht * Xi((r_ - gamma_).norm());
    }
    return ret;
}

void Billiard::calculateOmega()
{
    omega = A(theta) * thetap;
}


// Matriz de rotación A
Matrix<double> Billiard::A(const Vector<double>& theta_) const
{
    using namespace std;
    double c1 = cos(theta_[0]), c2 = cos(theta_[1]),
           s1 = sin(theta_[0]), s2 = sin(theta_[1]);
    return Matrix<double>(3,3,{0, c1 , s1*s2,
                               0, -s1, c1*s2,
                               1, 0  , c2    });
}


Matrix<double> Billiard::Aprime(const Vector<double>& theta_, const Vector<double>& thetap_) const
{
    double c1 = cos(theta_[0]), c2 = cos(theta_[1]),
           s1 = sin(theta_[0]), s2 = sin(theta_[1]);
    Matrix<double> A1 = Matrix<double>(3,3,{0, -s1 , c1*s2 ,
                                            0, -c1 , -s1*s2,
                                            0, 0   , 0   }),
                   A2 = Matrix<double>(3,3,{0, 0 , s1*c2,
                                            0, 0 , c1*c2,
                                            0, 0 , -s2});
                // A3 = 0
    return A1 * thetap_[0] + A2 * thetap_[1];
}


// Derivadas
Vector<double> Billiard::Dtheta(const Vector<double>& thetap_) const
{
    // Trivial
    return thetap_;
}


Vector<double> Billiard::Dv(const Vector<double>& r) const
{
    return 1./P.m * (Sxigamma - Sxi * r);
}

Vector<double> Billiard::DV(const Vector<double>& r,
                            const Vector<double>& Rcm,
                            const Vector<double>& theta,
                            const Vector<double>& Rp,
                            const Vector<double>& omega) const
{
    throw "Use direct call to Dv";
}

Vector<double> Billiard::Dthetap(const Vector<double>& r_,
                                 const Vector<double>& Rcm_,
                                 const Vector<double>& theta_,
                                 const Vector<double>& thetap_) const
{
    // Se calcula el torque
    Vector<double> tau = cross(Rcm_ - r_, Sxigamma) + cross(r_, Rcm_) * Sxi;
    // Matrices relevantes
    Matrix<double> A_ = A(theta_);
    Matrix<double> Agrad = Aprime(theta_, thetap_);
    Vector<double> omega = A_ * thetap_;

    using namespace std;
    double cosec2 = 1./sin(theta_[1]), sin1 = sin(theta_[0]);
    double cot2   = 1./tan(theta_[1]), cos1 = cos(theta_[0]);
    Matrix<double> Ainv(3,3, {-sin1*cot2 , -cos1*cot2 , 1,
                              cos1       , -sin1      , 0,
                              sin1*cosec2, cos1*cosec2, 0});
    // Fórmula larga
    return Ainv * (S.Iinv * (R * tau - cross(omega, S.I*omega)) - Agrad * thetap_);
}

Vector<double> Billiard::F(const Vector<double>& q, double t)
{
    Vector<double> r_ = {q[0],q[1],q[2]};
    Vector<double> R_ = {q[3],q[4],q[5]};
    Vector<double> theta_ = {q[6],q[7],q[8]};
    Vector<double> rp_ = {q[9],q[10],q[11]};
    Vector<double> Rp_ = {q[12],q[13],q[14]};
    Vector<double> thetap_ = {q[15],q[16],q[17]};

    // Obtención de parámetros generales. Algunos guardados como miembros en la clase
    Rt = lambdaT(theta_);
    R = Rt.transposed();
    generateIntegrationMesh(r_, R_);
    Sxi = calculateSxi(r_, R_);
    Sxigamma = calculateSxigamma(r_, R_);

    Vector<double> Dr_ = rp_;
    Vector<double> DR_ = Rp_;
    Vector<double> Dtheta_ = Dtheta(thetap_);
    Vector<double> Drp_ = Dv(r_);
    Vector<double> DRp_ = -P.m / S.mass * Drp_;
    Vector<double> Dthetap_ = Dthetap(r_, R_, theta_, thetap_);
    return Vector<double>{Dr_, DR_, Dtheta_, Drp_, DRp_, Dthetap_};
}



void Billiard::step()
{
    // Integración por rk4
    integrator.step();
    // Copiar valores a variables internas. Muy optimizable probablemente
    Vector<double> new_vals = integrator.getR();
    P.pos  = new_vals(0, 3);
    pos    = new_vals(3, 6);
    theta  = new_vals(6, 9);
    P.vel  = new_vals(9, 12);
    vel    = new_vals(12, 15);
    thetap = new_vals(15, 18);

    calculateOmega();
}



// Getters de variables
Vector<double> Billiard::getParticlePos() const
{
    return P.pos;
}
Matrix<double> Billiard::getRotation() const
{
    return Rt;
}
Vector<double> Billiard::getCM() const
{
    return pos;
}
double Billiard::getTime() const
{
    return integrator.getTime();
}


// Propiedades dinámicas
double Billiard::getE() const
{
    return P.getE() + S.mass * dot(vel,vel) / 2.0 + dot(omega, S.I * omega) / 2.0;
}


Vector<double> Billiard::getP() const
{
    return P.m * P.vel + S.mass * vel;
}


Vector<double> Billiard::getL() const
{
    return P.m * cross(P.pos, P.vel) + S.mass * cross(pos, vel) + Rt * S.I * omega;
}



// Cortes de Pointcare -----------------------------------------------------------------------------
DynamicBilliard::DynamicBilliard(const PointMass& particle, const Surface& surface,
                                 const Vector<double>& R0, const Vector<double>& V0, const Vector<double>& theta0, const Vector<double>& omega0,
                                 double height, double elasticity, unsigned int mesh_partition) :
    Billiard(particle, surface, R0, V0, theta0, omega0, height, elasticity, mesh_partition)
{
    particleAboveXYPlane = (Rt * particle.pos)[2] >= 0;
}


// Cálculo de cortes
void DynamicBilliard::step()
{
    Billiard::step();

    // Proyectar posición de partícula sobre sistema principal. Trasladar antes
    Vector<double> rp = R * (P.pos - pos);      // Ineficiente, da lo mismo supongo
    // Si hay cambio de signo en z, hay corte y se guarda
    if ( ! ((rp[2] >= 0) == particleAboveXYPlane) )
    {
        particleAboveXYPlane = rp[2] >= 0;
        pointcareCutsQueue.push(Vector<double>{rp[0], rp[1], (double)particleAboveXYPlane});
    }

}

// Flush de cortes a archivo
void DynamicBilliard::flushCutsToFile(std::ofstream& ofs)
{
    while (!pointcareCutsQueue.empty())
    {
        ofs << pointcareCutsQueue.front() << std::endl;
        pointcareCutsQueue.pop();
    }
}
