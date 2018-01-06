//
//  Polio_boxcar_parameters.h
//  
//
//  Created by Celeste on 12/11/17.
//
//

#ifndef Polio_boxcar_parameters_h
#define Polio_boxcar_parameters_h

#include <math.h>
#include <array>

using namespace std;


const int numComp = 4; //number of model compartments
const double numContacts = 100; //per year
const double birthRate = 0.02;

//parameter for Famulare probability of infection model
const double beta_fam = 14;    //mean from Famulare
const double alpha_fam = .44; //mean from Famulare
const double gamma_fam = .55;//mean from Famulare

//parameters for use in Teunis model

const double antGrowth = 0.2; //antibody growth rate
const double waningShape = 1.2;//r in Teunis model
const double waningRate = 1.5; //nu in Teunis model
array<double, numComp> intAnt = {1,10,100,1000}; //initial antibody concentration for S, R_1, R_2, R_3
array<double, numComp> peakAnt; //peak antibody concentration for S, R_1, R_2, R_3


//waning parameters for extended waning immunity model
double fastWaningRate(double peakAnt){//assumes fast waning period occurs in 50 days
    return (peakAnt*pow((1.0+(waningShape-1.0)*pow(peakAnt,(waningShape-1.0))*waningRate*50.0),(-1.0/(waningShape-1.0)))-peakAnt*pow((1.0+(waningShape-1.0)*pow(peakAnt,(waningShape-1.0))*waningRate*0.0),(-1.0/(waningShape-1.0))))/50.0;
}
double slowWaningRate(double peakAnt){//slow waning time window 50 days past infection to 1000 days past infection
    return (peakAnt*pow((1.0+(waningShape-1.0)*pow(peakAnt,(waningShape-1.0))*waningRate*1000.0),(-1.0/(waningShape-1.0)))-peakAnt*pow((1.0+(waningShape-1.0)*pow(peakAnt,(waningShape-1.0))*waningRate*50.0),(-1.0/(waningShape-1.0))))/950.0;
}

array<double, numComp> setPeakAnt(){ //10 fold boosting
    
    for(unsigned int i = 0; i < peakAnt.size(); i++){
        peakAnt[i] = intAnt[i]*11;
    }
    return peakAnt;
}


array<double,numComp> infectPeriod(){ //duration of infection for I_1, I_2, I_3, I_4 in days
    array<double,numComp> t1;
    for(unsigned int i = 0; i < t1.size(); i++){
        t1[i] = (1.0/antGrowth)*log(peakAnt[i]/intAnt[i]);
    }
    return t1;
}

double probInfection(double antibody){
    return 1-pow(1+(1000.0/beta_fam),-alpha_fam*pow(antibody,-gamma_fam));
}

/*array<double,numComp> probInfection(){//dose response model from Famulare with dose fixed at 1000 TCID50
    array<double,numComp> pInfection;
    for(unsigned int i = 0; i < intAnt.size(); i++){
    cout<<"antibody "<<antibody<<"\n";
    cout<<1-pow(1+(1000.0/beta_fam),(-alpha_fam*pow(intAnt[i],2)))<<"\n";
        pInfection[i] = 1-pow(1+(1000.0/beta_fam),-alpha_fam*pow(intAnt[i],-gamma_fam));
    //cout<<"prob infection "<<1-pow(1+(1000.0/beta_fam),-alpha_fam*pow(antibody,-gamma_fam))<<"\n";
    }
    return pInfection;
}*/



#endif /* Polio_boxcar_parameters_h */
