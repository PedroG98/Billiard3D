#ifndef BILLIARD_H
#define BILLIARD_H

#include "Matrix.h"
#include "NumericalIntegrator.h"

#include <queue>
#include <fstream>
#include <vector>
#include <utility>


// Masa puntual
struct PointMass {
    double m;
    Vector<double> pos;
    Vector<double> vel;

    // Getter
    double getE() const;
};



// Superficie paramétrica de la forma r(s,t)
class Surface
{
    private:
        // Función paramétrica
        Vector<double> (*r)(double s, double t);
        double (*da)(double s, double t);

        // Calcular momento de inercia. ASUME EJES PRINCIPALES DADOS POR AHORA
        Matrix<double> momentOfInertia() const;

    public:
        // Llamando curva y rango
        Surface(Vector<double> (*gamma)(double s, double t), double (*da)(double s, double t), double mass, double si, double sf, double ti, double tf);

        // Llamar
        Vector<double> operator()(double s, double t) const;
        double dA(double s, double t) const;

        // Rangos en s,t
        const double si,sf;
        const double ti,tf;

        // Momento de inercia
        double mass;
        const Matrix<double> I;
        const Matrix<double> Iinv;
};


// Clase para billiard generalizado
class Billiard
{
    protected:
        // Caja:
        // Superficie descriptiva
        Surface S;
        // Mesh de detección
        double partitionSizeT, partitionSizeS;
        std::vector<std::pair<double,double>> detectionMesh; // Estática
        std::vector<std::pair<double,double>> integrationMesh;         // Generada dinámicamente cada llamada a F

        // Propiedades dinámicas
        PointMass P;
        // Dentro de surface

        // Variables generalizadas DENTRO DE INTEGRADOR

        Vector<double> pos, vel;        // Posicion y velocidad de CM
        Vector<double> theta;           // Ángulos de Euler
        Vector<double> thetap;
        Vector<double> omega;           // Velocidad angular, utilitario

        // Integrador
        Rk4Integrator<Vector<double> (Billiard::*)(const Vector<double>&,double), Billiard> integrator;

        // Potencial
        double alpha;   // Elasticidad
        double U;       // Altura


        // Funciones auxiliares
        // Matriz de rotación
        Matrix<double> R, Rt;
        // Calcula gamma (curva) en punto s,t sobre la superficie (Sistema Lab)
        Vector<double> gamma(double s, double t, const Vector<double>& R) const;
        // Matriz de rotación de ángulos de Euler
        Matrix<double> lambdaT(const Vector<double>& theta_) const;

        // Integrales sobre superficie
        // Detecta posición y genera mesh en la cercanía, para integrar
        void generateIntegrationMesh(const Vector<double>& r,
                                     const Vector<double>& R);
        double Xi(double r) const;
        // Guardados en variables para evitar recalcular
        double Sxi;
        Vector<double> Sxigamma;
        double calculateSxi(const Vector<double>& r, const Vector<double>& R_) const;
        Vector<double> calculateSxigamma(const Vector<double>& r, const Vector<double>& R_) const;
        // Generar matriz generadora de omega w = A(theta)thetap
        Matrix<double> A(const Vector<double>& theta_) const;
        // Derivada temporal de A
        Matrix<double> Aprime(const Vector<double>& theta_, const Vector<double>& thetap_) const;
        // Genera omega. No interactuar con F por favor
        void calculateOmega();

        // Derivadas
        Vector<double> Dtheta(const Vector<double>& thetap_) const;
        Vector<double> Dv(const Vector<double>& r) const;
        Vector<double> DV(const Vector<double>& r,
                          const Vector<double>& R,
                          const Vector<double>& theta,
                          const Vector<double>& Rp,
                          const Vector<double>& omega) const;
        Vector<double> Dthetap(const Vector<double>& r,
                              const Vector<double>& R,
                              const Vector<double>& theta_,
                              const Vector<double>& thetap_) const;
        // qp = F(q) para los 18 componentes
        Vector<double> F(const Vector<double>& q, double t);

    public:
        Billiard(const PointMass& particle, const Surface& surface,
                 const Vector<double>& R0, const Vector<double>& V0, const Vector<double>& theta0, const Vector<double>& omega0,
                 double height, double elasticity, unsigned int mesh_partition);

        // Integración
        virtual void step();

        // Getters de variables
        Vector<double> getParticlePos() const;
        Matrix<double> getRotation() const;
        Vector<double> getCM() const;
        double getTime() const;

        // Getters de propiedades dinámicas
        double getE() const;
        Vector<double> getP() const;
        Vector<double> getL() const;

        // Para permitir la llamada de métodos privados
        friend class Rk4Integrator<Vector<double> (Billiard::*)(const Vector<double>&,double), Billiard>;
};


// Contiene corte de Pointcare en plano xy del sistema solidario
class DynamicBilliard : public Billiard
{
    protected:
        bool particleAboveXYPlane;

        // Cola con cortes de pointcare a flushear a output (Vectores 2D)
        std::queue<Vector<double>> pointcareCutsQueue;

    public:
        DynamicBilliard(const PointMass& particle, const Surface& surface,
        const Vector<double>& R0, const Vector<double>& V0, const Vector<double>& theta0, const Vector<double>& omega0,
        double height, double elasticity, unsigned int mesh_partition);

        // Agrega cálculo de cortes
        virtual void step();

        // Flush de cortes a archivo
        void flushCutsToFile(std::ofstream& ofs);
};



#endif // BILLIARD_H
