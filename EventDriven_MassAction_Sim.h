#ifndef EDMASSIM_H
#define EDMASSIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include <random>
#include <assert.h>
#include <array>
#include <fstream>
#include <string>




using namespace std;



class Event {
public:
    double time;
    string type;
    
   // Event(const Event& o) {time = o.time; type=o.type;}//*** this is a copy constructor
    Event(double t=0.0, string e="r"):time(t),type(e){} // *** this is a constructor
   // Event& operator=(const Event& o) {
   //     if(this != &o){
   //         time = o.time;
   //         type=o.type;
   //     }
   //     return *this;
   // } //*** this is copy assignment operator
};

class compTime {
public:
    //two in case we use pointers at some point
    bool operator() (const Event* lhs, const Event* rhs) const {
        
        return (lhs->time>rhs->time); //(->) is a reference operator for pointers
    }
    
    bool operator() (const Event& lhs, const Event& rhs) const {
        return (lhs.time>rhs.time);
    }
};

class Person{
    friend class Event;
    private:
    //attributes
    double m_age;
    double m_timeSinceInfection; //infection means any contact with virus
    double m_titerLevel;
    char m_infectionStatus; //'S','I','V'
    Event event;
    double time =0.0;
    string type= "r";
    priority_queue<Event, vector<Event>, compTime > EventQ;
    
    public:
    //default constructor
    Person(double age=0, double timeSinceInfection=0.0, double titerLevel=1.0, char infectionStatus ='S'):m_age(age), m_timeSinceInfection(timeSinceInfection),m_titerLevel(titerLevel),m_infectionStatus(infectionStatus){
        
    }
    
    void updateEventQ(double time,string type){
        EventQ.push(Event(time,type));
    }
    Event getEventQ(){
        assert(!EventQ.empty());
        return EventQ.top();
    }
    
    void popEventQ(){
        EventQ.pop();
    }
    
    int isEmptyQ(){
        if(EventQ.empty()){
            return 1;
        }
        else{
            return 0;
        }
    }
    
    void setAge(double age){
        m_age = age;
    }
    void setTimeSinceInfection(double timeSinceInfection){
        m_timeSinceInfection +=timeSinceInfection;
    }
    void setTiterLevel(double titerLevel){
        m_titerLevel = titerLevel;//max titer level is 2048
    }
    void setInfectionStatus(char infectionStatus){
        m_infectionStatus = infectionStatus;
    }
    int getAge(){
        return m_age;
    }
    double getTimeSinceInfection(){
        return m_timeSinceInfection;
    }
    double getTiterLevel(){
        return m_titerLevel;
    }
    char getInfectionStatus(){
        return m_infectionStatus;
    }
    
    double probInfGivenDose(double dose){//dose will vary based on if it is WPV or OPV
        if(m_infectionStatus=='S'){
            return 1-pow((1.0+(dose/14.0)),-.44*pow(m_titerLevel,0.55));//put in mean and shape parameters
        }
        else{
            return 0.0;
        }
    }
    void waning(){
        m_titerLevel= std::max(1.0, 3.0*pow(m_timeSinceInfection,-.75));///what is baseline line immmunity one month post infection??--replace 3.0 with this
        return;
    }
    
    void print(int i){
        std::cout<<"Attributes for person: "<<(i+1)<<" : "<<m_age<<" , "<<m_timeSinceInfection<<" , "<<m_titerLevel<<" ,"<<m_infectionStatus<<"\n";
    }
};




class EventDriven_MassAction_Sim {
    public:
                                    // constructor
        EventDriven_MassAction_Sim(const int n, const double gamma, const double beta, const double kappa, const double rho, const double birth, const double death): N(n), rng((random_device())()), GAMMA(gamma), BETA(beta), KAPPA(kappa), RHO(rho), BIRTH(birth), DEATH(death), unif_real(0.0,1.0),unif_int(1,n-1){
            //static constexpr int arraySize=5;
        }

        const int N;// population size
        const double GAMMA;
        const double BETA;
        const double KAPPA;
        const double RHO;
        const double BIRTH;
        const double DEATH;
        const double infDose=5; //number of doses from WPV infection
        const double vaccDose=3; //number of doses from OPV vacc
        static const int arraySize=5;//**can't seem to make 5 a variable***
    

        exponential_distribution<double> exp_gamma;
        exponential_distribution<double> exp_beta;
        exponential_distribution<double> exp_rho;
        uniform_real_distribution<double> unif_real;
        uniform_int_distribution<int> unif_int;
        mt19937 rng;
    

        //create fixed array of number of individuals in population
        array<Person,arraySize> individual;
        array<double,1>previousTime{0};

        double Now=0.0;                 // Current "time" in simulation
        vector<double>finalTime;
    
        void runSimulation(){
            while(nextEvent()>=0){
                if(nextEvent()==0){
                    finalTime.push_back(Now);
                    break;
                }
                continue;
            }
        }
    
        void randomizePopulation(int k,int j){
        //temporary initial conditions: these can be changed at a later date
            for(int i=0;i<arraySize;i++){
                double rr=unif_real(rng);
                if(rr<0.3){
                    individual[i].setAge(1);
                }
                else if(rr<.6){
                    individual[i].setAge(2);
                }
                else{
                    individual[i].setAge(3);
                }
                individual[i].setTimeSinceInfection(unif_real(rng));
                individual[i].setTiterLevel(100.0);
            }
            for(int i=0;i<k;++i){
                infect(individual[i]);//does it matter which individuals in array are initially infected??
            }
            //vaccinate(individual[1+j]);
        }
        
        void printIndividual(){
            for(int k=0;k<arraySize;k++){
                individual[k].print(k);
            }
        }
    


        double FinalTime(){
            return finalTime[0];
        }

        int totalSusceptibles(){
            int sSum=0;
            for(int i=0;i<arraySize;i++){
                if(individual[i].getInfectionStatus()=='S'){
                    sSum+=1;
                }
            }
            return sSum;
        }
        int totalInfecteds(){
            int iSum=0;
            for(int i=0;i<arraySize;i++){
                if(individual[i].getInfectionStatus()=='I'){
                    iSum+=1;
                }
            }
            return iSum;
        }
        double totalInfectionRate(){
            double infectionRate = 0.0;
            for(int i=0;i<arraySize;i++){
                infectionRate+=individual[i].probInfGivenDose(infDose);
            }
            return infectionRate;
        }


        void infect(Person& individual) {
            assert(individual.getInfectionStatus()=='S');
            individual.setInfectionStatus('I');
            individual.setTimeSinceInfection(0.0);
            individual.setTiterLevel(std::min(2048.0,11.0*individual.getTiterLevel()));//boost 10 fold
            //time to recovery
            exponential_distribution<double> exp_gamma(GAMMA);//**temporary recovery rate
            double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (Tc < Tr) {     // does contact occur before recovery?
                individual.updateEventQ(Tc,"inf_c");
                Tc += exp_beta(rng);
            }
            individual.updateEventQ(Tr,"inf_r");
            //time to waning
            exponential_distribution<double> exp_rho(RHO);//**temp waning rate
            double Tw = exp_rho(rng) + Now;
            individual.updateEventQ(Tw, "inf_wane");
            return;
        }
        void vaccinate(Person& individual){
            assert(individual.getInfectionStatus()=='S');
            individual.setInfectionStatus('V');
            individual.setTimeSinceInfection(0.0);
            individual.setTiterLevel(std::min(2048.0,11.0*individual.getTiterLevel()));//boost 10 fold
            //time to recovery
            exponential_distribution<double> exp_gamma(GAMMA);//**temporary recovery rate
            double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (Tc < Tr) {     // does contact occur before recovery?
                individual.updateEventQ(Tc,"inf_c");
                Tc += exp_beta(rng);
            }
            individual.updateEventQ(Tr,"vacc_r");
            exponential_distribution<double> exp_rho(RHO);//**temp waning rate
            double Tw=exp_rho(rng) + Now;
            individual.updateEventQ(Tw, "vacc_wane");
            return;

        }
    
        int nextEvent() {
            int endLoop=0;
            for(int i=0;i<arraySize;i++){//check to see if all event queues are empty
                if(individual[i].isEmptyQ()==0){
                    endLoop++;
                    break;
                }
            }
            if(endLoop==0){
                return 0;
            }
            double minTime=1000000000;
            int minIndex=0;
            for(int i=0;i<arraySize;i++){//find event to occur by finding smallest time in all priority queues
                if (individual[i].isEmptyQ()==0){
                    Event event = individual[i].getEventQ();
                    if(event.time<minTime){
                        minTime=event.time;
                        minIndex=i;
                    }
                }
            }
            Event event = individual[minIndex].getEventQ();
            Now = event.time;
            for(int i=0;i<arraySize;i++){//update individual's demography before event (do we want to make aging an event?)
                individual[i].setTimeSinceInfection(Now-previousTime[0]);
                individual[i].setAge(individual[i].getAge()+(Now-previousTime[0]));
            }
            if(event.type=="inf_c"){//includes contact with infected and vaccinated individual (OPV)
                double r1 = unif_real(rng);
                double r2 = unif_int(rng);
                double sum = 0.0;
                double infectionSum=totalInfectionRate();
                if(r2<=totalSusceptibles()){//if there are enough susceptibles in the pop then a contact will lead to an infection
                    for(int i=0;i<arraySize;i++){
                        if(r1<(sum+(individual[i].probInfGivenDose(infDose)/infectionSum))){
                            infect(individual[i]);
                            break;
                        }
                        sum+=individual[i].probInfGivenDose(infDose)/infectionSum;
                       
                    }
                }
                
            }
          /*  else if(event.type=="vacc_c"){//contact with vaccinated individual
                double r1 = unif_real(rng);
                double r2 = unif_int(rng);
                if(r2<=totalSusceptibles()){
                    for(int i=0;i<arraySize;i++){
                        if(r1<individual[i].probInfGivenDose(vaccDose)){
                            infect(individual[i]);
                            break;
                        }
                        
                    }
                }
            }*/
            else if (event.type == "inf_r") {//recovery from vacc and WPV may be different
                individual[minIndex].setInfectionStatus('R');
            }
            else if(event.type=="vacc_r"){
                individual[minIndex].setInfectionStatus('R');
            }
            else if(event.type=="inf_wane"){
                individual[minIndex].waning();
                individual[minIndex].setInfectionStatus('S');
            }
            else if(event.type=="vacc_wane"){
                individual[minIndex].waning();
                individual[minIndex].setInfectionStatus('S');
            }
            previousTime[0]=event.time;
            individual[minIndex].popEventQ();
                
            
            return 1;
        }
    
};
#endif
