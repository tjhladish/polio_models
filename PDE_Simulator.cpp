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
   
    assert(argc==2);

//within-host parameters
    const double mu0 = 0.1841; //pathogen growth rate
    const double b0 = 1; //initial pathogen concentration
    const double y0 = 1; //initial antibody concentration
    const double y1 = pow(10,5);//peak antibody concentration
    const double t1 = 28; //duration of infection
    const double mu = (1/t1)*log(y1/y0); //antibody growth rate
    const double c = b0*((mu-mu0)/(y0*(exp((mu-mu0)*t1)-1))); //pathogen clearance rate
    const double r = 1.69; //immunity shape parameter
    const double nu = 1.41; //immunity waning rate
    
//between-host parameters
    const double K1 = 100; //I1 contact rate
    const double K2 = 50; //Ir contact rate
    const double alpha1 = 13; //I1 recovery parameter
    const double alpha2 = 31; //Ir recovery parameter
    const double omega = 0.2; //waning parameter
    const double C = 0.0001; //saturation constant
    const double delta = 0.02; //birth/death rate
    const double totalPop = 101; //total population size
    
//time parameters
    const double dt = 0.1; //time step
    //N*M must be a perfect square
    const unsigned int N = atoi(argv[1]); //delta t * N = final time of interest T (days)
    const unsigned int M = atoi(argv[1]); //delta t * M = final time since infection tau (days)
    //const double T = dt*N;
    const double tau = dt*M;
    
//vectors for linking functions
    vector<double> beta1(M+1,0.0);
    vector<double> beta2(M+1,0.0);
    vector<double> gamma1(M+1,0.0);
    vector<double> gamma2(M+1,0.0);
    vector<double> rho(M+1,0.0);
    
//vectors for within-host solution
    vector<double> pathogen(M);
    vector<double> antibody(M);
cout << "a\n";
    //for(double t = 0; t<tau; t += dt){
    for (unsigned int i = 0; i < M; ++i) {
        const double t = dt*i;
        double b;
        double y;
        if(t<=t1){
            b = b0*exp(mu0*t) - (c*y0*(exp(mu*t)-exp(mu0*t))/(mu-mu0));
            y = y0*exp(mu*t);
        }
        else{
            b = 0;
            y = y1*pow((1+(r-1)*pow(y1,r-1)*nu*(t-t1)),-(1/(r-1)));
        }
        pathogen[i] = b;
        antibody[i] = y;
    }
cout << "b\n";
//Construct linking functions
    for(unsigned int j =0; j<pathogen.size(); ++j){
        const double beta1Link = (K1*pathogen[j])/(pathogen[j]+C);
        const double beta2Link = (K2*pathogen[j])/(pathogen[j]+C);
        const double gamma1Link = (alpha1*antibody[j])/(pathogen[j]+C);
        const double gamma2Link = (alpha2*antibody[j])/(pathogen[j]+C);
        const double rhoLink = omega/(antibody[j]+C);
        
        beta1[j]  = beta1Link;
        beta2[j]  = beta2Link;
        gamma1[j] = gamma1Link;
        gamma2[j] = gamma2Link;
        rho[j]    = rhoLink;
    }
//vectors for between-host solution
    vector<double> S(N,0.0);
    S[0] = 100;
    vector<double> I1(N*M,0.0);
    vector<double> R(N*M,0.0);
    vector<double> Ir(N*M,0.0);
    vector<double> symptomaticIncidence (N,0.0);
    symptomaticIncidence[0] = 1;
//initialize between-host compartments
cout << "c\n";
    for(int i = 0; i<M;i++){
        I1[i] = 1;
        //R[i] = 0;
        //Ir[i] = 0;
    }

//Finite Difference Method
    // use backward Euler difference quotient to approximate time derivatives
    // approximate integrals using right end point rule
    
cout << "d\n";
    for(unsigned int k=0; k<(N-1); k++){//looping through time
        double intSum = 0;
        double intSum1 = 0;
        for(unsigned int j=0; j<M; j++){
            const unsigned int index = M*k+j;
            assert(beta1.size()>j);
            assert(beta2.size()>j);
            assert(gamma1.size()>j);
            assert(gamma2.size()>j);
            assert(I1.size()>index);
            assert(Ir.size()>index);
            intSum  += dt*(beta1[j]*I1[index] + beta2[j]*Ir[index]);
            intSum1 += dt*(gamma1[j]*I1[index]+gamma2[j]*Ir[index]);
            //intSum1,intsum get too large when N,M>20
        }
        
        assert(S.size()>k+1);
        S[k+1] = (S[k] + dt*delta*totalPop)/(1+dt*intSum*(1.0/totalPop)+dt*delta);
        
        double intSum2 =0;
        for(unsigned int j=0; j<M; j++){
            if(j==0){
                assert(R.size()>(M*(k+1))+j);
                R[(M*(k+1))+j] = intSum1;//boundary condition
            }
            assert(R.size()>(M*(k+1))+(j+1));
            assert(R.size()>(M*k)+j);
            assert(rho.size()>j+1);
            assert(rho.size()>j);
            assert(R.size()>(M*(k+1))+j);
            R[(M*(k+1))+(j+1)] = R[(M*k)+j]/(1 + dt*((rho[j+1]/totalPop)*intSum+delta));
            intSum2+= dt*rho[j]*R[(M*(k+1))+j];
        }
        double intSum3=0;
        for(unsigned int j=0; j<M; j++){
            if(j==0){
                //boundary conditions
                assert(I1.size()>(M*(k+1))+j);
                assert(S.size()>k+1);
                assert(Ir.size()>(M*(k+1))+j);
                I1[(M*(k+1))+j] = intSum*S[k+1]/totalPop;
                Ir[(M*(k+1))+j] = (1.0/totalPop)*intSum2*intSum;
            }
           
            assert(I1.size()>(M*(k+1))+(j+1));
            assert(I1.size()>(M*k)+j);
            assert(gamma1.size()>j+1);
            assert(Ir.size()>(M*(k+1))+(j+1));
            assert(Ir.size()>(M*k)+j);
            assert(gamma2.size()>j+1);
            assert(beta1.size()>j);
            assert(I1.size()>(M*k)+j);
            assert(beta2.size()>j);
            assert(Ir.size()>(M*k)+j);
            
            I1[(M*(k+1))+(j+1)] = I1[(M*k)+j]/(1+dt*(gamma1[j+1]+delta));
            Ir[(M*(k+1))+(j+1)] = Ir[(M*k)+j]/(1+dt*(gamma2[j+1]+delta));
            intSum3 += dt*(beta1[j]*I1[(M*k)+j]+beta2[j]*Ir[(M*k)+j]);
            
        }
       
        assert(symptomaticIncidence.size()>k+1);
        assert(S.size()>k+1);
        symptomaticIncidence[k+1] = S[k+1]*intSum3/totalPop;
    }
cout << "e\n";
    /*cout<<"\n R vec\n";
    cout<<"R vec size "<<R.size()<<"\n";
    for(unsigned int i =0; i<R.size();i++){
        cout<<R[i]<<" ";
    }
    cout<<"\n I1 vec\n";
    cout<<"I1 vec size "<<I1.size()<<"\n";
    for(unsigned int i =0; i<I1.size();i++){
        cout<<I1[i]<<" ";
    }
    cout<<"\n Ir vec\n";
    cout<<"Ir vec size "<<Ir.size()<<"\n";
    for(unsigned int i =0; i<Ir.size();i++){
        cout<<Ir[i]<<" ";
    }
    cout<<"\n S vec\n";
    cout<<"S vec size "<<S.size()<<"\n";
    for(unsigned int i =0; i<S.size();i++){
        cout<<S[i]<<" ";
    }
    cout<<"\n symptomatic incidence\n";
    cout<<"symptomatic incidence size "<<symptomaticIncidence.size()<<"\n";
    for(unsigned int i =0; i<symptomaticIncidence.size();i++){
        cout<<symptomaticIncidence[i]<<" ";
    }
    cout<<"endl";*/
    return 0;
}
