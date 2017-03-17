#include "EventDriven_MassAction_Sim.h"

<<<<<<< HEAD


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
=======
int main() { 

    int N        = 1000000;
    double GAMMA = 1/(double) 6;
    double BETA  = 1/(double) 4;

    for(int i=0; i<5; i++ ) {
        EventDriven_MassAction_Sim sim(N,GAMMA,BETA);
        sim.rand_infect(100);
        sim.run_simulation();
        cout << sim.epidemic_size() << endl;
    }
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45

    return 0;
}
