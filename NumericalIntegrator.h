/* Clases que aplican esquemas de integración numéricos a una ecuación diferencial ordinaria
 * Reciben parámetros, función y valor inicial
 */
#ifndef NUMERICALINTEGRATOR_H
#define NUMERICALINTEGRATOR_H

#include "Vector.h"
#include <utility>      // Para std::pair
#include <type_traits>  // Function template check



// Integrador rk4 para una EDO de la forma x' = F(x, t)
// sobre member Parent
template <class Function, class Parent>
class Rk4Integrator
{
    private:
        // Función integrada
        Function F;
        Parent *P;

        // Tiempo actual y paso
        double t;
        double dt;
        double maxDt;

        // Posicion actual
        Vector<double> q;

        // Para paso variable
        bool adaptative;
        double tolerance;
        void adjustTimeStep();

        // Paso rk4 con paso explícito. Retorna valor de nuevo r, no modifica!
        Vector<double> rk4Step(double _dt) const;

    public:
        // Constructor toma función, posición inicial y paso inicial
        // Adicionalmente opción para paso variable con tolerancia dada
        Rk4Integrator();
        Rk4Integrator(Function F, Parent& P, const Vector<double>& r0, double dt,
                      bool adaptative=false, double tolerance = 0.0, double max_step = 0.0);

        // Dar paso de simulación
        void step();

        // Getters
        // Posición actual
        Vector<double> getR() const;
        Vector<double>& getR();
        // Tiempo actual
        double getTime() const;
};


// Integrador PEFRL para un sistema con fuerza F=F(r,t)
class PEFRLIntegrator
{
    private:
        // Fuerza
        Vector<double> (*F)(const Vector<double>&, double);

        // Tiempo actual y paso
        double t;
        double dt;

        // Posicion velocidad y masa
        Vector<double> r;
        Vector<double> v;
        double m;

    public:
        // Constructor toma función, posición y velocidad inicial, masa y paso
        PEFRLIntegrator(Vector<double> (*F)(const Vector<double>&, double), Vector<double>& r0, Vector<double>& v0, double mass, double dt);

        // Dar paso PEFRL
        void step();

        // Getters
        // Posición
        Vector<double> getR() const;
        // Velocidad
        Vector<double> getV() const;
        // Tiempo
        double getTime() const;
};






// Definiciones ----------------------------------------------------------------------------------------------
// Integrador Rk4 --------------------------------------------------------------------------------------
// Constructor

template <class Function, class Parent>
Rk4Integrator<Function, Parent>::Rk4Integrator()
{}

template <class Function, class Parent>
Rk4Integrator<Function, Parent>::Rk4Integrator(Function F, Parent& P, const Vector<double>& r0, double dt,
                                               bool adaptative, double tolerance, double max_step) :
    F(F), P(&P), t(0), dt(dt), q(r0), adaptative(adaptative), tolerance(tolerance), maxDt(max_step)
{
    // Revisar que Function tenga la signatura adecuada
    //static_assert(std::is_function<Function>());
    //std::result_of<Function(const Vector<double>&, double)> r;
    //std::cout << std::is_same<decltype(r), Vector<double>>::value << std::endl;
}


// Paso rk4, retorna nuevo r sin modificarlo
template <class Function, class Parent>
Vector<double> Rk4Integrator<Function, Parent>::rk4Step(double _dt) const
{
    Vector<double> k1 = P->F(q                , t);
    Vector<double> k2 = P->F(q + k1 * _dt/2.0 , t + _dt/2.0);
    Vector<double> k3 = P->F(q + k2 * _dt/2.0 , t + _dt/2.0);
    Vector<double> k4 = P->F(q + k3 * _dt     , t + _dt);

    return q + dt / 6 * (k1 + 2.0*k2 + 2.0*k3 + k4);
}


// Ajustar paso según tolerancia, según los valores para rk4 utilizados en clases
template <class Function, class Parent>
void Rk4Integrator<Function, Parent>::adjustTimeStep()
{
    // Se calcula la nueva posicion para paso completo y para medio paso
    Vector<double> fullStep = rk4Step(dt);
    // Para la segunda llamada de half step se ajusta el tiempo y posicion
    Vector<double> tmp = q;
    q = rk4Step(dt / 2);
    t += dt / 2;
    Vector<double> halfStep = rk4Step(dt / 2);
    // Deshace avance
    t -= dt / 2;
    q = tmp;

    // Error y parámetros S1, S2
    double S1 = 0.9, S2 = 4.0;
    double error = (fullStep - halfStep).norm();
    // t_test
    double dtTest = std::pow(tolerance / error, 0.2) * dt;

    // Nuevo paso con ajuste ideal
    dt = dtTest;
    if (dt > maxDt)
        dt = maxDt;
    if (S1 * dtTest < dt / S2)
        dt /= S2;
    else if (S1 * dtTest > dt * S2)
        dt *= S2;
    else
        dt = S1 * dtTest;
}




// Dar paso, incluye ajuste de ser necesario
template <class Function, class Parent>
void Rk4Integrator<Function, Parent>::step()
{
    // Si hay paso variable, ajustar
    if (adaptative)
        adjustTimeStep();

    // Llama a rk4 y modifica las variables
    q = rk4Step(dt);
    t += dt;
}

// Getters

template <class Function, class Parent>
Vector<double> Rk4Integrator<Function, Parent>::getR() const
{
    return q;
}
template <class Function, class Parent>
Vector<double>& Rk4Integrator<Function, Parent>::getR()
{
    return q;
}
template <class Function, class Parent>
double Rk4Integrator<Function, Parent>::getTime() const
{
    return t;
}


#endif // NUMERICALINTEGRATOR_H
