#include "EventDriven_MassAction_Sim.h"




int main() {
    int numSims = 1000;
    int N        = 100;
    double GAMMA = 13;//1/(double) 6;  //***castes 6 to a double for double division
    double BETA  = 200;//1/(double) 4;
    double KAPPA = .8434;
    double RHO = .02;
    double BIRTH = .02;
    double DEATH = .02;

    
    double vecSumm=0.0;
    
    for(int i=0; i<numSims; i++ ) {
        EventDriven_MassAction_Sim sim(N,GAMMA,BETA,KAPPA,RHO,BIRTH,DEATH);
        sim.rand_infect(1);
        sim.run_simulation();
        //cout << sim.epidemic_size() << endl;
       vecSumm+=sim.FinalTime();
    }
    cout<<(vecSumm/numSims)<<"\n";//avg time to extinction


    return 0;
}
