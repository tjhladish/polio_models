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
        EventDriven_MassAction_Sim(const int n, const double gamma, const double beta, const double kappa, const double rho, const double birth, const double death): rng((random_device())()), GAMMA(gamma), BETA(beta), KAPPA(kappa), RHO(rho), BIRTH(birth), DEATH(death), unif_real(0.0,1.0),unif_int(1,n-1){
            previousTime = {0};
            people = vector<Person*>(n);
            for (Person* &p: people) p = new Person;
            Now=0.0;
        }

        ~EventDriven_MassAction_Sim() {
            for (Person* &p: people) delete p;
        }

        const double GAMMA;
        const double BETA;
        const double KAPPA;
        const double RHO;
        const double BIRTH;
        const double DEATH;
        const double infDose=5; //number of doses from WPV infection
        const double vaccDose=3; //number of doses from OPV vacc
    

        exponential_distribution<double> exp_gamma;
        exponential_distribution<double> exp_beta;
        exponential_distribution<double> exp_rho;
        uniform_real_distribution<double> unif_real;
        uniform_int_distribution<int> unif_int;
        mt19937 rng;
    

        //create fixed array of number of individuals in population
        vector<Person*> people;
        array<double,1>previousTime;

        double Now;                 // Current "time" in simulation
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
            for(Person* p: people) {
                double rr=unif_real(rng);
                if(rr<0.3){
                    p->setAge(1);
                }
                else if(rr<.6){
                    p->setAge(2);
                }
                else{
                    p->setAge(3);
                }
                p->setTimeSinceInfection(unif_real(rng));
                p->setTiterLevel(100.0);
            }
            for(int i=0;i<k;++i){
                infect(people[i]);//does it matter which individuals in array are initially infected??
            }
            //vaccinate(individual[1+j]);
        }
        
        void printPeople(){
            for(unsigned int i = 0; i<people.size(); ++i) {
                people[i]->print(i);
            }
        }
    

        double FinalTime(){
            return finalTime[0];
        }

        int totalSusceptibles(){
            int sSum=0;
            for(Person* p: people) {
                if(p->getInfectionStatus()=='S'){
                    sSum+=1;
                }
            }
            return sSum;
        }
        int totalInfecteds(){
            int iSum=0;
            for(Person* p: people) {
                if(p->getInfectionStatus()=='I'){
                    iSum+=1;
                }
            }
            return iSum;
        }
        double totalInfectionRate(){
            double infectionRate = 0.0;
            for(Person* p: people) {
                infectionRate+=p->probInfGivenDose(infDose);
            }
            return infectionRate;
        }


        void infect(Person* p) {
            assert(p->getInfectionStatus()=='S');
            p->setInfectionStatus('I');
            p->setTimeSinceInfection(0.0);
            p->setTiterLevel(std::min(2048.0,11.0*p->getTiterLevel()));//boost 10 fold
            //time to recovery
            exponential_distribution<double> exp_gamma(GAMMA);//**temporary recovery rate
            double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (Tc < Tr) {     // does contact occur before recovery?
                p->updateEventQ(Tc,"inf_c");
                Tc += exp_beta(rng);
            }
            p->updateEventQ(Tr,"inf_r");
            //time to waning
            exponential_distribution<double> exp_rho(RHO);//**temp waning rate
            double Tw = exp_rho(rng) + Now;
            p->updateEventQ(Tw, "inf_wane");
            return;
        }
        void vaccinate(Person* p){
            assert(p->getInfectionStatus()=='S');
            p->setInfectionStatus('V');
            p->setTimeSinceInfection(0.0);
            p->setTiterLevel(std::min(2048.0,11.0*p->getTiterLevel()));//boost 10 fold
            //time to recovery
            exponential_distribution<double> exp_gamma(GAMMA);//**temporary recovery rate
            double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (Tc < Tr) {     // does contact occur before recovery?
                p->updateEventQ(Tc,"inf_c");
                Tc += exp_beta(rng);
            }
            p->updateEventQ(Tr,"vacc_r");
            exponential_distribution<double> exp_rho(RHO);//**temp waning rate
            double Tw=exp_rho(rng) + Now;
            p->updateEventQ(Tw, "vacc_wane");
            return;

        }
    
        int nextEvent() {
            int endLoop=0;
            for(Person* p: people) {//check to see if all event queues are empty
                if(p->isEmptyQ()==0){
                    endLoop++;
                    break;
                }
            }
            if(endLoop==0){
                return 0;
            }
            double minTime=1000000000;
            int minIndex=0;
            for(unsigned int i = 0; i < people.size(); ++i) {//find event to occur by finding smallest time in all priority queues{
                Person* p = people[i];
                if (p->isEmptyQ()==0){
                    Event event = p->getEventQ();
                    if(event.time<minTime){
                        minTime=event.time;
                        minIndex=i;
                    }
                }
            }
            Event event = people[minIndex]->getEventQ();
            Now = event.time;
            for(Person* p: people) {//update individual's demography before event (do we want to make aging an event?)
                p->setTimeSinceInfection(Now-previousTime[0]);
                p->setAge(p->getAge()+(Now-previousTime[0]));
            }
            if(event.type=="inf_c"){//includes contact with infected and vaccinated individual (OPV)
                double r1 = unif_real(rng);
                double r2 = unif_int(rng);
                double sum = 0.0;
                double infectionSum=totalInfectionRate();
                if(r2<=totalSusceptibles()){//if there are enough susceptibles in the pop then a contact will lead to an infection
                    for(Person* p: people) {
                        if(r1<(sum+(p->probInfGivenDose(infDose)/infectionSum))){
                            infect(p);
                            break;
                        }
                        sum+=p->probInfGivenDose(infDose)/infectionSum;
                       
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
                people[minIndex]->setInfectionStatus('R');
            }
            else if(event.type=="vacc_r"){
                people[minIndex]->setInfectionStatus('R');
            }
            else if(event.type=="inf_wane"){
                people[minIndex]->waning();
                people[minIndex]->setInfectionStatus('S');
            }
            else if(event.type=="vacc_wane"){
                people[minIndex]->waning();
                people[minIndex]->setInfectionStatus('S');
            }
            previousTime[0]=event.time;
            people[minIndex]->popEventQ();
            
            return 1;
        }
    
};
#endif
