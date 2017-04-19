#include "EventDriven_MassAction_Sim.h"


int main() {
    int numSims = 50;
    int N        = 100;
    //double GAMMA = 13;
    
    //double KAPPA = .8434;
    //double RHO = .02;
    double BIRTH = .008;//Pakistan death rate per year
    double DEATH = .008;
    double BETA  = 15;
    double BETAENV = 2;//contact environment once per day
for(int j=0;j<10;j++){

    double vecSumm=0.0;
    double firstInfecSum=0.0;
    vector<double> R0;
    
    for(int i=0; i<numSims; i++) {
        EventDriven_MassAction_Sim sim(N,BETA,BIRTH,DEATH, BETAENV);
        sim.randomizePopulation(1);//randomize with # of infecteds 
        //sim.printPeople();
        sim.runSimulation();
        //sim.printPeople();
        vecSumm+=sim.FinalTime();
        firstInfecSum+=sim.avgFirstInfec();
        R0.push_back(sim.R0());
    }
        double meanR0 = 0.0;
        for(int k=0;k<R0.size();k++){
            meanR0 += R0[k];
            cout<<"in R0 loop\n";
        }
        cout<<"diff in avg R0 "<<((meanR0/numSims)-4)<<"\n";
        if(abs((meanR0/numSims)-8)<.001){
            BETA = BETA;
            BETAENV=BETAENV;
            cout<<"keeping same\n";
        }
        else{
            BETA -= .8*((meanR0/numSims)-8);
            BETAENV-= .2*((meanR0/numSims)-8);
        }
   // cout<<"avg sim time: "<<vecSumm/numSims<<"\n";
   // cout<<"avg age at first infection: "<<firstInfecSum/numSims<<"\n";
    //cout<<"avg num infected: "<<R0Sum/numSims<<"\n";
   std::ofstream myfile1;
    myfile1.open ("/Users/Celeste/Desktop/N=100_R0EST_10sims_test_7.csv");
    myfile1<<vecSumm/numSims<<"\n";
    myfile1<<firstInfecSum/numSims<<"\n";
    for(int i=0;i<R0.size();i++){
        myfile1<<R0[i]<<"\n";
    }
    myfile1.close();
        cout<<"BETA "<<BETA<<"\n";
        cout<<"BETAENV "<<BETAENV<<"\n";
    }

    return 0;
}
