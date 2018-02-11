#ifndef CHOLERA_SIR_H
#define CHOLERA_SIR_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math.h>
#include "DiffEq_Sim.h"
#include <array>

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
        double vecCounter = 0;
        double prevTime = 0;

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
    
    double getRain(double t){
        array<double,57> rain = {(267.029+45.824+30.495), 26.68579, 169.00367, 275.74283, 418.06171, 91.20199, 140.94853, 140.94853, 456.00995, 39.24, 98.1, 215.82, 307.38, 18.086, 280.333, 253.204, 352.677, 39.537, 237.222, 263.58, 777.561, 68.464, 136.928, 179.718, 470.69, 11.488, 183.808, 333.152, 620.352, 164.065, 94.985, 198.605, 405.845, 30.348, 106.218, 227.61, 394.524, 10.091, 171.547, 272.457, 544.914, 40.552, 172.346, 263.588, 537.314, 16.198, 121.485, 251.069, 421.148, 110.964, 206.076, 467.634, 515.19, 13.714, 143.997, 198.853, 322.279};
        
        if((int)floor(t)%3 and !(int)floor(prevTime)%3){
            vecCounter++;
        }
        prevTime = t;
        return rain[vecCounter]/3.0;//rainfall given for three month time increments -- converting to monthly
    }
    
    double getDelta(double t){
        double modTime = (int)floor(t)%12;
        array<double,12> meanRain = {74.892,8.407,14.214,28.9,46.386,76.279,43.7,60.8,94.15,102.729,207.964,196.021};//assuming rainfall starts in december (one month delay)
        
        return meanRain[modTime];
    }

        void derivative(double const x[], double dxdt[]) {
            const double S = x[0];
            const double I = x[1];
            const double Y = x[2];
            const double R = x[3];
            const double N = S+I+Y+R;
            const double _t = get_time();
            double waterTerm = getRain(_t)/(getRain(_t)+getDelta(_t));

            dxdt[0] = (b*N) - waterTerm*(beta*S*(I+Y)/N) - (mu*S) + (rho*Y) + (epsilon*R);
            dxdt[1] = (waterTerm*C*beta*S*(I+Y)/N) - (gamma*I) - (mu*I);
            dxdt[2] = (waterTerm*(1-C)*beta*S*(I+Y)/N) - (rho*Y) - (mu*Y);
            dxdt[3] = (gamma*I) - (epsilon*R) - (mu*R);
        }

};

#endif
