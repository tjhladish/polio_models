//
//  PDE_Simulator.cpp
//  
//
//  Created by Celeste on 10/1/17.
//
//
#include "AbcSmc.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <array>
//#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_math.h>
//#include "SIR_H.h"

using namespace std;
const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    
//within-host parameters
    double mu0 = 0.1841; //pathogen growth rate
    double b0 = 1; //initial pathogen concentration
    double y0 = 1; //initial antibody concentration
    double y1 = pow(10,5);//peak antibody concentration
    double t1 = 28; //duration of infection
    double mu = (1/t1)*log(y1/y0); //antibody growth rate
    double c = b0*((mu-mu0)/(y0*(exp((mu-mu0)*t1)-1))); //pathogen clearance rate
    double r = args[0];//1.69; //immunity shape parameter
    double nu = args[1];//1.41; //immunity waning rate
    
//between-host parameters
    double K1 = 100; //I1 contact rate
    double K2 = 50; //Ir contact rate
    double alpha1 = 13; //I1 recovery parameter
    double alpha2 = 31; //Ir recovery parameter
    double omega = 0.2; //waning parameter
    double C = 0.0001; //saturation constant
    double delta = 0.02; //birth/death rate
    double totalPop = 101; //total population size
    
//time parameters
    double dt = 0.1; //time step
    //N*M must be a perfect square
    double N = 5; //delta t * N = final time of interest T (days)
    double M = 5; //delta t * M = final time since infection tau (days)
    double T = dt*N;
    double tau = dt*M;
    

    
//vectors for linking functions
    vector<double> beta1;
    vector<double> beta2;
    vector<double> gamma1;
    vector<double> gamma2;
    vector<double> rho;
    
//vectors for within-host solution
    vector<double> pathogen;
    vector<double> antibody;
    for(double i = 0; i<tau;i = i+ dt){
        double b;
        double y;
        if(i<=t1){
            b = b0*exp(mu0*i) - (c*y0*(exp(mu*i)-exp(mu0*i))/(mu-mu0));
            y = y0*exp(mu*i);
        }
        else{
            b = 0;
            y = y1*pow((1+(r-1)*pow(y1,r-1)*nu*(i-t1)),-(1/(r-1)));
        }
        pathogen.push_back(b);
        antibody.push_back(y);
    }

//Construct linking functions
    for(unsigned int j =0; j<pathogen.size();j++){
        double beta1Link = (K1*pathogen[j])/(pathogen[j]+C);
        double beta2Link = (K2*pathogen[j])/(pathogen[j]+C);
        double gamma1Link = (alpha1*antibody[j])/(pathogen[j]+C);
        double gamma2Link = (alpha2*antibody[j])/(pathogen[j]+C);
        double rhoLink = omega/(antibody[j]+C);
        
        beta1.push_back(beta1Link);
        beta2.push_back(beta2Link);
        gamma1.push_back(gamma1Link);
        gamma2.push_back(gamma2Link);
        rho.push_back(rhoLink);
    }
//vectors for between-host solution
    vector<double> S(N);
    vector<double> I1(N*M);
    vector<double> R(N*M);
    vector<double> Ir(N*M);
    vector<double> symptomaticIncidence (N);
//initialize between-host compartments
    S[0] = 100;
    for(int i = 0; i<N;i++){
        I1[i] = 1;
        R[i] = 0;
        Ir[i] = 0;
    }
    symptomaticIncidence[0] = 1;

//Finite Difference Method
    // use backward Euler difference quotient to approximate time derivatives
    // approximate integrals using right end point rule
    
    for(int k=0; k<(N-1); k++){//looping through time
        double intSum = 0;
        double intSum1 = 0;
        for(int s=0; s<M; s++){
            intSum += dt*(beta1[s]*I1[sqrt(I1.size())*k+s] + beta2[s]*Ir[sqrt(Ir.size())*k+s]);
            intSum1 += dt*(gamma1[s]*I1[sqrt(I1.size())*k+s]+gamma2[s]*Ir[sqrt(Ir.size())*k+s]);
            //intSum1,intsum get too large when N,M>20
        }
        
        S[k+1] = (S[k] + dt*delta*totalPop)/(1+dt*intSum*(1/totalPop)+dt*delta);
        
        double intSum2 =0;
        for(int j=0; j<M; j++){
            if(j==0){
                R[(sqrt(R.size())*(k+1))+j] = intSum1;//boundary condition
            }
            R[(sqrt(R.size())*(k+1))+(j+1)] = R[k+j]/(1 + dt*((rho[j+1]/totalPop)*intSum+delta));
            intSum2+= dt*rho[j]*R[(sqrt(R.size())*(k+1))+j];
        }
        double intSum3=0;
        for(int j=0; j<M; j++){
            if(j==0){
                //boundary conditions
                I1[(sqrt(I1.size())*(k+1))+j] = intSum*S[k+1]/totalPop;
                Ir[(sqrt(Ir.size())*(k+1))+j] = (1/totalPop)*intSum2*intSum;
            }
            
            I1[(sqrt(I1.size())*(k+1))+(j+1)] = I1[(sqrt(I1.size())*k)+j]/(1+dt*(gamma1[j+1]+delta));
            Ir[(sqrt(Ir.size())*(k+1))+(j+1)] = Ir[(sqrt(Ir.size())*k)+j]/(1+dt*(gamma2[j+1]+delta));
            intSum3 += dt*(beta1[j]*I1[(sqrt(I1.size())*k)+j]+beta2[j]*Ir[(sqrt(Ir.size())*k)+j]);
            
        }
        
        symptomaticIncidence[k+1] = S[k+1]*intSum3/totalPop;
        }
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
    cout<<"symptomatic incdience size "<<symptomaticIncidence.size()<<"\n";
    for(unsigned int i =0; i<symptomaticIncidence.size();i++){
        cout<<symptomaticIncidence[i]<<" ";
    }
    cout<<"endl";*/
    vector<double> metrics = {365*symptomaticIncidence.back()};
    return metrics;
}

void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --process --simulate -n <number of simulations per database write>\n\n";

}


int main(int argc, char* argv[]) {

    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) { 
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {  
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {  
            buffer_size = atoi(argv[++i]);
        } else {
            usage(); 
            exit(101);
        }
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    } 

    if (simulate_db) {
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
