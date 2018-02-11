#ifndef CHOLERA_SIR_H
#define CHOLERA_SIR_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math.h>
#include "DiffEq_Sim.h"

using namespace std;

class Cholera_Sim : public DiffEq_Sim {

    private:
        const double b;
        const double beta;
        const double C;
        const double epsilon;
        const double gamma;
        const double mu;
        const double rho;

    public:
        Cholera_Sim() : b(0.0), beta(0.0), C(0.0), epsilon(0.0), gamma(0.0), mu(0.0), rho(0.0) { nbins=4; }
        Cholera_Sim(double _b, double _beta, double _C, double _epsilon, double _gamma, double _mu, double _rho): 
                    b(_b), beta(_beta), C(_C), epsilon(_epsilon), gamma(_gamma), mu(_mu), rho(_rho) { nbins=4; }
        ~Cholera_Sim() {};

        void initialize(double S, double I, double Y, double R) {
            x = new double[nbins];
            x[0] = S; 
            x[1] = I; 
            x[2] = Y;
            x[3] = R;
        }

        void derivative(double const x[], double dxdt[]) {
            const double S = x[0];
            const double I = x[1];
            const double Y = x[2];
            const double R = x[3];
            const double N = S+I+Y+R;
            const double _t = get_time();

            dxdt[0] = (b*N) - sin(4*_t/(2*M_PI))*(beta*S*(I+Y)/N) - (mu*S) + (rho*Y) + (epsilon*R);
            dxdt[1] = (sin(4*_t/(2*M_PI))*C*beta*S*(I+Y)/N) - (gamma*I) - (mu*I);
            dxdt[2] = (sin(4*_t/(2*M_PI))*(1-C)*beta*S*(I+Y)/N) - (rho*Y) - (mu*Y);
            dxdt[3] = (gamma*I) - (epsilon*R) - (mu*R);
        }

};

#endif
