//
//  Polio_boxcar_simulator.cpp
//  
//
//  Created by Celeste on 12/10/17.
//
//

#include <stdio.h>
#include "Polio_boxcar_model.h"
#include "DIFFEQ_SIM.h"
#include <math.h>
#include "Polio_boxcar_parameters.h"

void usage(){
    cout<<"\n\t Usage: Input model initialization <S> <I1> <I2> <I3> <I4> <R1> <R2> <R3> <R4>\n";
    exit(-1);
}

int main(int argc, char** argv){
    if(argc!= 10) usage();
    
    //Determine peak antibody level for each infected compartment
    setPeakAnt();

    //Determine infectious period for each infected compartment
    infectPeriod();
    
    
    //Constuct model and parameters
    Polio_boxcar_model model(0.02,numContacts*probInfection(1),numContacts*probInfection(10),numContacts*probInfection(100),numContacts*probInfection(1000),1/infectPeriod()[0],1/infectPeriod()[1],1/infectPeriod()[2],1/infectPeriod()[3],.3,.3,.3);
    
    
    //Initialize model using user inputs
    model.initialize(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]));
    
    //Run model
    model.run_simulation();
    
    //Output simulation results
    model.printX();

    
    
    return 0;
}
