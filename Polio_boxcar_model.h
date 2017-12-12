//
//  SIR_H.h
//  
//
//  Created by Celeste on 10/1/17.
//
//

#ifndef Polio_boxcar_model_h
#define Polio_boxcar_model_h

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include "DIFFEQ_SIM.h"
#include <math.h>

using namespace std;

class Polio_boxcar_model : public DiffEq_Sim {
    
private:
    //birth/death rate
    const double mu;
    //effective contact rates
    const double b1;
    const double b2;
    const double b3;
    const double b4;
    //recovery rates
    const double g1;
    const double g2;
    const double g3;
    const double g4;
    //waning rates
    const double r1;
    const double r2;
    const double r3;
    
public:
    Polio_boxcar_model() : mu(0.0), b1(0.0), b2(0.0),b3(0.0),b4(0.0), g1(0.0), g2(0.0), g3(0.0), g4(0.0), r1(0.0), r2(0.0), r3(0.0) {nbins=9;}
    Polio_boxcar_model(double bd, double infect1, double infect2, double infect3, double infect4, double rec1, double rec2, double rec3, double rec4, double wane1, double wane2, double wane3): mu(bd), b1(infect1), b2(infect2), b3(infect3), b4(infect4), g1(rec1), g2(rec2), g3(rec3), g4(rec4), r1(wane1), r2(wane2), r3(wane3) {nbins=9;}
    ~Polio_boxcar_model() {};
    
    void initialize(double S, double I1, double I2, double I3, double I4, double R1, double R2, double R3, double R4) {
        x = new double[nbins];
        x[0] = S;
        x[1] = I1;
        x[2] = I2;
        x[3] = I3;
        x[4] = I4;
        x[5] = R1;
        x[6] = R2;
        x[7] = R3;
        x[8] = R4;
    }
    
    void derivative(double const x[], double dxdt[]) {
        //key:
        //x[0] -> S         x[5] -> R_1
        //x[1] -> I_1       x[6] -> R_2
        //x[2] -> I_2       x[7] -> R_3
        //x[3] -> I_3       x[8] -> R_4
        //x[4] -> I_4
        
        int totalPop = 0;
        int infectPop = 0;
        for(int i = 0; i < nbins; i++){
            totalPop += x[i];
            if(i >0 and i < 5){
                infectPop += x[i];
            }
        }
        dxdt[0] = mu*totalPop - (b1*infectPop/(double)totalPop + mu)*x[0];
        dxdt[1] = b1*x[0]*infectPop/(double)totalPop - (g1 + mu)*x[1];
        dxdt[2] = b2*x[5]*infectPop/(double)totalPop - (g2 + mu)*x[2];
        dxdt[3] = b3*x[6]*infectPop/(double)totalPop - (g3 + mu)*x[3];
        dxdt[4] = b4*x[7]*infectPop/(double)totalPop - (g4 + mu)*x[4];
        dxdt[5] = g1*x[1] - (b2*infectPop/(double)totalPop + mu)*x[5] + r1*x[6];
        dxdt[6] = g2*x[2] - (b3*infectPop/(double)totalPop + r1 + mu)*x[6] + r2*x[7];
        dxdt[7] = g3*x[3] - (b4*infectPop/(double)totalPop + r2 + mu)*x[7] + r3*x[8];
        dxdt[8] = g4*x[4] - (r3 + mu)*x[8];
        
    }
    void printX(){for(int i=0; i < nbins; i++) { cout << x[i] << " ";} cout << endl; }
    
};


#endif /* Polio_boxcar_model_h */
