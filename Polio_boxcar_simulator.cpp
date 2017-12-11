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

void usage(){
    cout<<"\n\t Usage: Input model initialization <S> <I1> <I2> <I3> <I4> <R1> <R2> <R3> <R4>\n";
    exit(-1);
}

int main(int argc, char** argv){
    if(argc!= 10) usage();
    Polio_boxcar_model model(0.02,1,1,1,1,.5,.5,.5,.5,.003,.003,.003,.003);
    model.initialize(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),atoi(argv[9]));
    model.run_simulation();
    model.printX();

    
    
    return 0;
}
