#include "EventDriven_Sim_exp_recovery.hpp"


int main() {
    const int numSims = 1;
    const int N        = 10000;
    const double GAMMA = 13;
    
    const double KAPPA = .4179;
    const double RHO = .2;
    const double BIRTH = 0.02*N;
    const double DEATH = 0.02;
    const double BETA  =135*2;

    std::vector<double> TTE;
    double vecSumm=0.0;
    vector<double> R0;
    double R0sum = 0.0;
    //int seed = 0;
    std::array<double,10000> AVG_I1vec={0};
    std::array<double,10000> AVG_Irvec={0};
    std::array<double,10000> AVG_Svec={0};
    std::array<double,10000> AVG_Rvec={0};
    std::array<double,10000> AVG_Pvec={0};
    int ii=0;
    int births=0;
    int deaths=0;
    int I1inf=0;
    int Irinf=0;
    int rec=0;
    int wane=0;
    
    for(int i=0; i<numSims; i++) {
        EventDriven_MassAction_Sim sim(N,BETA,BIRTH,DEATH,GAMMA,RHO,KAPPA);
        sim.randomizePopulation(2);//randomize with # of infecteds
        sim.runSimulation();
        TTE.push_back(sim.FinalTime());
        vecSumm+=sim.FinalTime();
        ii++;
        R0sum+= sim.NumInfected();
        /*for(int k = 0;k<sim.printVectorI1().size();k++){
                    AVG_I1vec[k]= (AVG_I1vec[k]*i+sim.printVectorI1()[k])/(double)(i+1);
        }
        for(int k = 0;k<sim.printVectorIr().size();k++){
                    AVG_Irvec[k]= (AVG_Irvec[k]*i+sim.printVectorIr()[k])/(double)(i+1);
        }
        for(int k = 0;k<sim.printVectorS().size();k++){
                    AVG_Svec[k]= (AVG_Svec[k]*i+sim.printVectorS()[k])/(double)(i+1);
        }
        for(int k = 0;k<sim.printVectorR().size();k++){
                    AVG_Rvec[k]= (AVG_Rvec[k]*i+sim.printVectorR()[k])/(double)(i+1);
        }
        for(int k = 0;k<sim.printVectorP().size();k++){
                    AVG_Pvec[k]= (AVG_Pvec[k]*i+sim.printVectorP()[k])/(double)(i+1);
        }*/
        for(int k = 0;k<sim.printVectorI1().size();k++){
           // cout<<"k"<<k<<"\n";
           // cout<<"AVG I1 vec[k] before "<<AVG_I1vec[k]<<"\n";
           // cout<<"sim.printVectorI1 "<<sim.printVectorI1()[k]<<"\n";
                if(AVG_I1vec[k]==0){
                    AVG_I1vec[k]=sim.printVectorI1()[k];
                  //  cout<<"in zero loop\n";
                  //  cout<<"AVG I1 vec after "<<AVG_I1vec[k]<<"\n";
                }
                else{
                    if(sim.printVectorI1()[k]!=0){
                        AVG_I1vec[k]= (AVG_I1vec[k]*i+sim.printVectorI1()[k])/(double)(i+1);
                      //  cout<<"in nonzero loop\n";
                      //  cout<<"AVG I1 vec after "<<AVG_I1vec[k]<<"\n";
                    }
                }
        }
        for(int k = 0;k<sim.printVectorIr().size();k++){
            if(AVG_Irvec[k]==0){
                AVG_Irvec[k]=sim.printVectorIr()[k];
            }
            else{
                if(sim.printVectorIr()[k]!=0){
                    AVG_Irvec[k]= (AVG_Irvec[k]*i+sim.printVectorIr()[k])/(double)(i+1);
                }
            }
        }
        for(int k = 0;k<sim.printVectorS().size();k++){
            if(AVG_Svec[k]==0){
                AVG_Svec[k]=sim.printVectorS()[k];
            }
            else{
                if(sim.printVectorS()[k]!=0){
                    AVG_Svec[k]= (AVG_Svec[k]*i+sim.printVectorS()[k])/(double)(i+1);
                }
            }
        }
        for(int k = 0;k<sim.printVectorR().size();k++){
            if(AVG_Rvec[k]==0){
                AVG_Rvec[k]=sim.printVectorR()[k];
            }
            else{
                if(sim.printVectorR()[k]!=0){
                    AVG_Rvec[k]= (AVG_Rvec[k]*i+sim.printVectorR()[k])/(double)(i+1);
                }
            }
        }
        for(int k = 0;k<sim.printVectorP().size();k++){
            if(AVG_Pvec[k]==0){
                AVG_Pvec[k]=sim.printVectorP()[k];
            }
            else{
                if(sim.printVectorP()[k]!=0){
                    AVG_Pvec[k]= (AVG_Pvec[k]*i+sim.printVectorP()[k])/(double)(i+1);
                }
            }
        }
        births+=sim.NumBirths();
        deaths+=sim.NumDeaths();
        I1inf+=sim.NumI1Inf();
        Irinf+=sim.NumIRInf();
        rec+=sim.NumRec();
        wane+=sim.NumWane();
    }
    
    std::ofstream myfile7;
    myfile7.open ("/Users/Celeste/Desktop/C++IBM/N=10000,beta=135,fast_I1vec_valid_time=5_nonhomopop_births.csv");
    for(int i=0;i<AVG_I1vec.size();++i){
        myfile7<<AVG_I1vec[i]<<"\n";
    }
    myfile7.close();
    std::ofstream myfile8;
    myfile8.open ("/Users/Celeste/Desktop/C++IBM/N=10000,beta=135,fast_Irvec_valid_time=5_nonhomopop_births.csv");
    for(int i=0;i<AVG_Irvec.size();++i){
        myfile8<<AVG_Irvec[i]<<"\n";
    }
    myfile8.close();
    
    std::ofstream myfile9;
    myfile9.open ("/Users/Celeste/Desktop/C++IBM/N=10000,beta=135,fast_Svec_valid_time=5_nonhomopop_births.csv");
    for(int i=0;i<AVG_Svec.size();++i){
        myfile9<<AVG_Svec[i]<<"\n";
    }
    myfile9.close();
    std::ofstream myfile10;
    myfile10.open ("/Users/Celeste/Desktop/C++IBM/N=10000,beta=135,fast_Rvec_valid_time=5_nonhomopop_births.csv");
    for(int i=0;i<AVG_Rvec.size();++i){
        myfile10<<AVG_Rvec[i]<<"\n";
    }
    myfile10.close();
    std::ofstream myfile11;
    myfile11.open ("/Users/Celeste/Desktop/C++IBM/N=10000,beta=135,fast_Pvec_valid_time=5_nonhomopop_births.csv");
    for(int i=0;i<AVG_Pvec.size();++i){
        myfile11<<AVG_Pvec[i]<<"\n";
    }
    myfile11.close();


    cout<<"TTE "<<vecSumm/numSims<<"\n";
    
    cout<<"num births "<<births/(double)numSims<<"\n";
    cout<<"num deaths "<<deaths/(double)numSims<<"\n";
    cout<<"num I1 inf "<<I1inf/(double)numSims<<"\n";
    cout<<"num Ir inf "<<Irinf/(double)numSims<<"\n";
    cout<<"num rec "<<rec/(double)numSims<<"\n";
    cout<<"num wane "<<wane/(double)numSims<<"\n";


    return 0;
}
