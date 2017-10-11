//
//  PDE_Simulator.cpp
//  
//
//  Created by Celeste on 10/1/17.
//
//

#include <iostream>
#include <math.h>
#include <assert.h>
#include <vector>
#include <array>
//#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_math.h>
//#include "SIR_H.h"

using namespace std;

int main(int argc, char** argv){
   
    assert(argc==3);

//within-host parameters
    const double mu0 = 0.1841; //pathogen growth rate
    const double b0 = 1; //initial pathogen concentration
    const double y0 = 1; //initial antibody concentration
    const double y1 = pow(10,5);//peak antibody concentration
    const double t1 = 28; //duration of infection (days)
    const double mu = (1/t1)*log(y1/y0); //antibody growth rate
    const double c = b0*((mu-mu0)/(y0*(exp((mu-mu0)*t1)-1))); //pathogen clearance rate
    const double r = 1.69; //immunity shape parameter
    const double nu = 1.41; //immunity waning rate
    
//between-host parameters
    const double K1 = (100/(double)365); //I1 daily contact rate
    const double K2 = (50/(double)365); //Ir daily contact rate
    const double alpha1 = (13/(double)365); //I1 daily recovery parameter
    const double alpha2 = (31/(double)365); //Ir daily recovery parameter
    const double omega = (0.2/(double)365); //daily waning parameter
    const double C = 0.0001; //saturation constant
    const double delta = (0.02/(double)365); //daily birth/death rate
    //const double totalPop = 101; //total population size
    
//time parameters
    const double dt = 0.1; //time step
    const unsigned int N = atoi(argv[1]); //delta t * N = final time of interest T (days)
    const unsigned int M = atoi(argv[2]); //delta t * M = final time since infection tau (days)

//convergence parameters
    //const double epsilon = 0.001;
    //const unsigned int eq_interval = 10000;
    //pair<unsigned int, double> obs_min;
    //pair<unsigned int, double> obs_max;
    //const double T = dt*N;
    //const double tau = dt*M;
    
//vectors for linking functions
    vector<double> beta1(M,0.0);
    vector<double> beta2(M,0.0);
    vector<double> gamma1(M,0.0);
    vector<double> gamma2(M,0.0);
    vector<double> rho(M,0.0);
    
//vectors for within-host solution
//Construct linking functions
    for (unsigned int j = 0; j < M; ++j) {
        const double t = dt*j;
        double pathogen;
        double antibody;
        if(t<=t1){
            pathogen = b0*exp(mu0*t) - (c*y0*(exp(mu*t)-exp(mu0*t))/(mu-mu0));
            antibody = y0*exp(mu*t);
        }
        else{
            pathogen = 0;
            antibody = y1*pow((1+(r-1)*pow(y1,r-1)*nu*(t-t1)),-(1/(r-1)));
        }
        const double beta1Link = (K1*pathogen)/(pathogen+C);
        const double beta2Link = (K2*pathogen)/(pathogen+C);
        const double gamma1Link = (alpha1*antibody)/(pathogen+C);
        const double gamma2Link = (alpha2*antibody)/(pathogen+C);
        double rhoLink = 0;
        if(j>(t1/dt)){
            rhoLink = omega/(antibody+C);
        }
        beta1[j]  = beta1Link;
        beta2[j]  = beta2Link;
        gamma1[j] = gamma1Link;
        gamma2[j] = gamma2Link;
        rho[j]    = rhoLink;
    }
//vectors for between-host solution
    vector<double> S(N,0.0);
    S[0] = 100;
    vector<vector<double> > I1(N,vector<double>(M,0.0));
    vector<vector<double> > R(N,vector<double>(M,0.0));
    vector<vector<double> > Ir(N,vector<double>(M,0.0));
    vector<double> symptomaticIncidence (N,0.0);
//initialize between-host compartments
    I1[0][0] = 1;
    symptomaticIncidence[0] = 1.0/101;
    //obs_min = {0,symptomaticIncidence[0]};
    //obs_max = {0,symptomaticIncidence[0]};

//Finite Difference Method
    // use backward Euler difference quotient to approximate time derivatives
    // approximate integrals using right end point rule
    
    for(unsigned int k=1; k<N; k++){//looping through time (rows)
        //calculate the total population
        double I1pop = 0;
        double Rpop =0;
        double Irpop=0;
        for(unsigned int j = 0; j < M; j++){
            I1pop += I1[k-1][j];
            Rpop += R[k-1][j];
            Irpop += Ir[k-1][j];
        }
        double totalPop = S[k-1] + I1pop + Rpop + Irpop;
        assert(totalPop>0);
        double intSum = 0;
        double intSum1 = 0;
        for(unsigned int j=0; j<M; j++){//looping through time since infection (columns)
            //linearize to get k-1 index
            intSum  += dt * (beta1[j] * I1[k-1][j] + beta2[j] * Ir[k-1][j]);
            intSum1 += dt * (gamma1[j] * I1[k-1][j] + gamma2[j] * Ir[k-1][j]);
            //intSum1,intsum get too large when N,M>20
        }
        S[k] = (S[k-1] + dt*delta*totalPop)/(1+dt*intSum*(1.0/totalPop)+dt*delta);//
        double intSum2 =0;
        for(unsigned int j=0; j<M; j++){//columns
            if(j==0){
                R[k][j] = intSum1;//boundary condition
            } else {
                R[k][j] = R[k-1][j-1]/(1 + dt*((rho[j]/totalPop)*intSum+delta));
            }
            intSum2+= dt*rho[j]*R[k][j];
        }
        double intSum3=0;
        for(unsigned int j=0; j<M; j++){//columns
            if(j==0){
                //boundary conditions
                I1[k][j] = intSum*S[k]/totalPop;
                Ir[k][j] = (1.0/totalPop)*intSum2*intSum;
            } else {
                I1[k][j] = I1[k-1][j-1]/(1+dt*(gamma1[j]+delta));
                Ir[k][j] = Ir[k-1][j-1]/(1+dt*(gamma2[j]+delta));
            }
            intSum3 += dt*(beta1[j]*I1[k-1][j]+beta2[j]*Ir[k-1][j]);
            
        }
       
        symptomaticIncidence[k] = S[k]*intSum3/totalPop;
    }

/*    cout<<"\nR vec size "<<R.size()<<"\n";
    for (auto v: R) { for (auto e: v) cout << e << " "; }
    cout<<"\nI1 vec size "<<I1.size()<<"\n";
    for (auto v: I1) { for (auto e: v) cout << e << " "; }
    cout<<"\nIr vec size "<<Ir.size()<<"\n";
    for (auto v: Ir) { for (auto e: v) cout << e << " "; }
    cout<<"\nS vec size "<<S.size()<<"\n";
    for (auto e: S) { cout << e << " "; }*/
    //cout<<"\nsymptomatic incidence size "<<symptomaticIncidence.size()<<"\n";
    for (auto e: symptomaticIncidence) { cout << e << endl; }
    //cout<<endl;
    return 0;
}
