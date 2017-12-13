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

//parameters for use in Teunis model

const double mu = 0.2; //antibody growth rate
array<double, numComp> intAnt = {1,10,100,1000}; //initial antibody concentration for S, R_1, R_2, R_3
array<double, numComp> peakAnt; //peak antibody concentration for S, R_1, R_2, R_3

array<double, numComp> setPeakAnt(){ //10 fold boosting
    
    for(unsigned int i = 0; i < peakAnt.size(); i++){
        peakAnt[i] = intAnt[i]*11;
    }
    return peakAnt;
}


array<double,numComp> infectPeriod(){ //duration of infection for I_1, I_2, I_3, I_4 in days
    array<double,numComp> t1;
    for(unsigned int i = 0; i < t1.size(); i++){
        t1[i] = (1.0/mu)*log(peakAnt[i]/intAnt[i]);
    }
    return t1;
}

double probInfection(double antibody){//dose response model from Famulare with dose fixed at 1000 TCID50
    double beta = 14;    //mean from Famulare
    double alpha = .44; //mean from Famulare
    double gamma = .55;//mean from Famulare
    
    return 1-pow(1+(1000.0/beta),-alpha*pow(antibody,-gamma));
}

#endif /* Polio_boxcar_parameters_h */
