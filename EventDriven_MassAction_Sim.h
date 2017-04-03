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
#include <time.h>
#include <math.h>
#include "EventDriven_parameters.hpp"




using namespace std;



class Person{
    friend class Event;
    private:
    //attributes
    double m_age;
    double m_timeSinceInfection; //infection means any contact with virus
    double m_titerLevel;
    char m_infectionStatus; //'S','I','V'
    
    public:
    //default constructor
    Person(double age=0, double timeSinceInfection=0.0, double titerLevel=1.0, char infectionStatus ='S'):m_age(age), m_timeSinceInfection(timeSinceInfection),m_titerLevel(titerLevel),m_infectionStatus(infectionStatus){
        
    }
    
    void setAge(double age){
        m_age = age;
    }
    void updateAge(double age){
        m_age +=age;
    }
    double futureAge(double t){
        return m_age+t;
    }
    double deathTime(){
        return maxAge-m_age;
    }
    void reset(){
        m_age = 0.0;
        m_titerLevel = 1.0;
        m_timeSinceInfection=0.0;
        m_infectionStatus = 'S';
    }
    
    void setTimeSinceInfection(double timeSinceInfection){
        m_timeSinceInfection = timeSinceInfection;
    }
    
    void updateTimeSinceInfection(double timeSinceInfection){
        cout<<"in update time since infection loop\n";
        cout<<"double time since infection "<< timeSinceInfection<<"\n";
        cout<<"updated m_timeSinceInfection "<<m_timeSinceInfection+timeSinceInfection<<"\n";
        m_timeSinceInfection +=timeSinceInfection;
    }
    void setTiterLevel(double titerLevel){
        m_titerLevel = std::min(maxTiter,titerLevel);
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
            return 1-pow((1.0+(dose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));//put in mean and shape parameters
        }
        else{
            return 0.0;
        }
    }
    void waning(){
        if(m_timeSinceInfection>=(1/(double)12)){
            m_titerLevel= std::max(1.0, Nab1*pow((m_timeSinceInfection*(double)12),-waningLambda));
        }
        return;
    }
    double peakShedding(){//assumes NAb=1?
        if(m_age<(7/(double)24)){
            return pow(10.0,Smax);
        }
        else{
            return pow(10.0,((Smax-Smin)*exp((7-m_age)/tau)+Smin));
        }
    }
    
    double shedding(double t){//***t needs to be time individual was infected
        if(m_infectionStatus=='I'){
            return .5*erfc((log(t)-(log(muWPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaWPV));
        }
        else if(m_infectionStatus=='V'){
            return .5*erfc((log(t)-(log(muOPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaOPV));
        }
        else{
            return 0;
        }
    }
    
   /* double sheddingOPV(double t){
        return .5*erfc((log(t)-(log(muOPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaOPV));
    }*/
    
    double stoolViralLoad(double t){
        return std::max(pow(10.0,2.6),(peakShedding()/k*log2(m_titerLevel))*(exp(eta-(pow(nu,2)/(double)2)-(pow(log(t)-eta,2)/(double)2*pow(nu+xsi*log(t),2)))/t));
        //**units are in log_10 TCID50/g
    }
    
    void print(int i){
        std::cout<<"Attributes for person: "<<(i+1)<<" : "<<m_age<<" , "<<m_timeSinceInfection<<" , "<<m_titerLevel<<" ,"<<m_infectionStatus<<"\n";
    }
};

class Event {
    friend class Person;
public:
    double time;
    string type;
    Person* people;

    
    // Event(const Event& o) {time = o.time; type=o.type;}//*** this is a copy constructor
    Event(double t=0.0, string e="r",Person* p = 0):time(t),type(e),people(p){} // *** this is a constructor
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


class EventDriven_MassAction_Sim {
    public:
                                    // constructor
        EventDriven_MassAction_Sim(const int n, const double gamma, const double beta, const double kappa, const double rho, const double birth, const double death): rng((random_device())()), GAMMA(gamma), BETA(beta), KAPPA(kappa), RHO(rho), BIRTH(birth), DEATH(death), unif_real(0.0,1.0),unif_int(1,n-1){
            previousTime = {0};
            people = vector<Person*>(n);
            for (Person* &p: people) p = new Person;
            Now=0.0;
            Environment=100000.0;
            numInfected=0;
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
    
        exponential_distribution<double> exp_gamma;
        exponential_distribution<double> exp_beta;
        exponential_distribution<double> exp_rho;
        exponential_distribution<double> exp_death;
        exponential_distribution<double> exp_betaEnvironment;
        exponential_distribution<double> exp_checkEnvironment;
        uniform_real_distribution<double> unif_real;
        uniform_int_distribution<int> unif_int;
        mt19937 rng;
    
        //containers to keep track of various pieces of information
        vector<Person*> people;
        priority_queue<Event, vector<Event>, compTime > EventQ;
        array<double,1>previousTime={0};
        array<double,1>finalTime;
        vector<double> timeOfParalyticCase;

    
    
        double Now;                 // Current "time" in simulation
        double Environment;
        int numInfected;
    

    
        void runSimulation(){
            while(nextEvent()>=0){
                if(nextEvent()==0){
                    finalTime[0]=Now;
                    break;
                }
                continue;
            }
        }
    
        void randomizePopulation(int k,int j){
        //temporary initial conditions: these can be changed at a later date
            discrete_distribution<int> age {0,1,1,1};
            for(Person* p: people) {
                p->setAge(age(rng));
                p->setTimeSinceInfection(unif_real(rng));
                p->setTiterLevel(1000.0);
                death(p); //set when each individual will die
            }
            for(int i=0;i<k;++i){
                infect(people[i]);//does it matter which individuals in array are initially infected??
            }
            //vaccinate(individual[1+j]);
            //time to testing the environment (put a person as a placeholder--nothing happens to this person)
            exponential_distribution<double> checkEnvironment(GAMMA);//**temp check rate--needs to change
            double Te = checkEnvironment(rng);
            EventQ.push(Event(Te,"check_env",people[1]));
        }
        
        void printPeople(){
            for(unsigned int i = 0; i<people.size(); ++i) {
                people[i]->print(i);
            }
        }
    

        double FinalTime(){
            return TTE;
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
      /*  int totalInfecteds(){
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
        }*/


        void infect(Person* p) {
            assert(p->getInfectionStatus()=='S');
            cout<<"infection\n";
        //    cout<<"age : "<<p->getAge()<<"\n";
        //    cout<<"titer level "<<p->getTiterLevel()<<"\n";
        //    cout<<"time since infection "<<p->getTimeSinceInfection()<<"\n";
            numInfected++;
            p->setInfectionStatus('I');
            p->setTimeSinceInfection(0.0);
     //       cout<<"time since infection after set "<<p->getTimeSinceInfection()<<"\n";
     //       cout<<"titer level before set: "<<p->getTiterLevel()<<"\n";
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
    //        cout<<"titer level after set: "<<p->getTiterLevel()<<"\n";
            //is this a paralytic case?
            double r1 = unif_real(rng);
            if(r1<PIR){
                timeOfParalyticCase.push_back(Now);
            }
            //time to recovery
            exponential_distribution<double> exp_gamma(GAMMA);//**temporary recovery rate
            double Tr = exp_gamma(rng) + Now;
            // time to next human-human contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (Tc < Tr) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,"inf_c",p));
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,"inf_r",p));
            return;
        }
        void vaccinate(Person* p){
            assert(p->getInfectionStatus()=='S');
            p->setInfectionStatus('V');
            p->setTimeSinceInfection(0.0);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            //time to recovery
            exponential_distribution<double> exp_gamma(GAMMA);//**temporary recovery rate
            double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (Tc < Tr) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,"inf_c",p));
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,"vacc_r",p));
            return;

        }
    
        void infectByEnvironment(Person* p){
            assert(p->getInfectionStatus()=='S');
            numInfected++;
            p->setInfectionStatus('I');
            p->setTimeSinceInfection(0.0);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            //is this a paralytic case?
            double r1 = unif_real(rng);
            if(r1<PIR){
                timeOfParalyticCase.push_back(Now);
                cout<<"paralytic case\n";
            }
            //time to recovery
            exponential_distribution<double> exp_gamma(DEATH);//**temporary recovery rate
            double Tr = exp_gamma(rng) + Now;
            // time to next human-human contact
            exponential_distribution<double> exp_beta(DEATH);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (Tc < Tr) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,"inf_c",p));
                Tc += exp_beta(rng);
                cout<<"contact\n";
            }
            EventQ.push(Event(Tc,"inf_r",p));
            cout<<"recovery\n";
            return;
        }
        void death(Person* p){
            exponential_distribution<double> exp_death(DEATH);
            double deathAge = exp_death(rng);
            double Td = deathAge + Now;
            if(p->futureAge(deathAge)>100.0){
                Td = p->deathTime();
            }
            cout<<"death time " <<Td<<"\n";
            EventQ.push(Event(Td,"death",p));
            return;
        }
    
        double environmentalSurveillance(){
            if((Environment/(latrineVolume+Environment))>detectionRate){
                cout<<"detected pathogen in water!!!!!!!!!!!!!\n";
                return 0;
            }
            else{
                return 1;
            }
        }
    
        int nextEvent() {
            if(numInfected==0 or EventQ.empty()) return 0;
            Event event = EventQ.top();
            Now = event.time;
            cout<<"Now "<<Now<<"\n";
            cout<<"event "<<event.type<<"\n";
            for(int i=0;i<people.size();i++){
                cout<<"person "<<(i+1)<<" titer level "<<people[i]->getTiterLevel()<<"\n";
            }
            for(Person* p: people) {
                p->updateAge((Now-previousTime[0]));
                p->updateTimeSinceInfection((Now-previousTime[0]));
                // does this person contact the environment?
                exponential_distribution<double> exp_betaEnvironment(DEATH);//**temporary contact rate
                double rand2 = unif_real(rng);
                double Tce = exp_betaEnvironment(rng)+Now;
                if(rand2 < exp_betaEnvironment(rng)){
                    EventQ.push(Event(Tce,"env_c",p));
                }
           /*     double rand = unif_real(rng);
                if(rand<p->shedding(p->getTimeSinceInfection()) and p->getTiterLevel()!=0){
                    Environment+=128*(365*(Now-previousTime[0]))*p->stoolViralLoad(p->getTimeSinceInfection()*365);//avg 46720 g of feces per year (128 g/day)
                  //  cout<<"Now - previous time: "<<Now-previousTime[0]<<"\n";
                  //  cout<<"Environment: "<<128*(365*(Now-previousTime[0]))*p->stoolViralLoad(p->getTimeSinceInfection()*365)<<"\n";
                }*/
            }
            Person* individual = event.people;
            if(event.type=="inf_c"){//includes contact with infected and vaccinated individual (OPV)
                cout<<"in contact\n";
                double r1 = unif_real(rng);
                double r2 = unif_int(rng);
                if(r2<=totalSusceptibles()){//if there are enough susceptibles in the pop then a contact will lead to an infection
                    for(Person* p: people) {
                        cout<<"time since infection: "<<p->getTimeSinceInfection()<<"\n";
                      //  cout<<"update to time since infection: "<<abs(Now-(p->getTimeSinceInfection()))<<"\n";
                     //   cout<<"time since infection after update: "<<p->getTimeSinceInfection()<<"\n";
                     //   cout<<"titer level before waning: "<<p->getTiterLevel()<<"\n";
                        p->waning();
                    //    cout<<"titer level after waning: "<<p->getTiterLevel()<<"\n";
                        if(r1<(p->probInfGivenDose(infDose))){//infects the first person with low enough titer level**does this matter?
                            infect(p);
                            break;
                        }
                       
                    }
                }
                
            }
            else if(event.type=="vacc_c"){//contact with vaccinated individual
                double r1 = unif_real(rng);
                double r2 = unif_int(rng);
                if(r2<=totalSusceptibles()){
                    for(Person* p: people){
                        p->updateTimeSinceInfection(abs(Now-(p->getTimeSinceInfection())));
                        p->waning();
                        if(r1<p->probInfGivenDose(vaccDose)){
                            infect(p);
                            break;
                        }
                        
                    }
                }
            }
            else if (event.type == "inf_r") {//recovery from vacc and WPV may be different
                individual->setInfectionStatus('S');
                TTE=Now;
                numInfected--;
            }
            else if(event.type=="vacc_r"){
                individual->setInfectionStatus('S');
                TTE=Now;
                numInfected--;
            }
            else if(event.type=="death"){
                cout<<"death\n";
                individual->reset();//keeps population constant
                death(individual);//give each new person a death
            }
            else if(event.type=="env_c"){
                //first update Environment
                latrineVolume+=people.size()*feces*urine;//assume constant for now
                for(Person* p: people){
                    if(p->getTiterLevel()!=0){
                        Environment+=gramsFeces*p->shedding((Now-p->getTimeSinceInfection()));
                    }
                }
                double rand1=unif_real(rng);
                exponential_distribution<double> virusDeath(delta);
                if(rand1<virusDeath(rng) and Environment>=0){
                    Environment-=inactivationRate;
                }
                double rand2 = unif_real(rng);
                if(rand2<(Environment/(latrineVolume+Environment)) and individual->getInfectionStatus()=='S'){
                    infectByEnvironment(individual);
                }
                
            }
            else if(event.type=="check_env"){
                environmentalSurveillance();
            }
            previousTime[0]=event.time;
            EventQ.pop();
            
            return 1;
        }
    
};
#endif
