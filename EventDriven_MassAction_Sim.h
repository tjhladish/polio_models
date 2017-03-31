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




using namespace std;



class Person{
    friend class Event;
    private:
    //attributes
    double m_age;
    double m_timeSinceInfection; //infection means any contact with virus
    double m_titerLevel;
    char m_infectionStatus; //'S','I','V'
    //Event event;
    //double time =0.0;
    //string type= "r";
    //priority_queue<Event, vector<Event>, compTime > EventQ;
    
    public:
    //default constructor
    Person(double age=0, double timeSinceInfection=0.0, double titerLevel=0.0, char infectionStatus ='S'):m_age(age), m_timeSinceInfection(timeSinceInfection),m_titerLevel(titerLevel),m_infectionStatus(infectionStatus){
        
    }
    
    void setAge(double age){
        m_age = age;
    }
    void updateAge(double age){
        m_age +=age;
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
        if(m_timeSinceInfection!=0){
        m_titerLevel= std::max(1.0, 3.0*pow(m_timeSinceInfection,-.75));
        }///what is baseline line immmunity one month post infection??--replace 3.0 with this
        else{
            m_titerLevel=0;
        }
        return;
    }
    double peakShedding(){
        if(m_age<(7/(double)24)){
            return 6.7;//I put in the mean***need to change
        }
        else{
            return ((6.7-4.3)*exp((7-m_age)/(10/(double)24))+4.3);//I put in the mean***need to change
        }
    }
    
    double shedding(double t){//***t needs to be time individual was infected
        return .5*erfc((log(t)-(log(30.3)-log(1.16)*log(m_titerLevel)))/sqrt(2.0)*log(1.86));
    }//I put in the mean of these parameters***need to change
    
    double stoolViralLoad(double t){
        return std::max(2.6,((1.0-.056*log2(m_titerLevel))*log(peakShedding()))*(exp(1.65-(pow(.17,2)/(double)2)-(pow(log(t)-1.65,2)/(double)2*pow(.17+.32*log(t),2)))/t));
        //I put in the mean of these parameters***need to change**units are in log_10
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
        const double PIR = 0.001; //type 3 paralysis incidence rate
        double TTE =0;
        double delta = 1/(double)86;
    

        exponential_distribution<double> exp_gamma;
        exponential_distribution<double> exp_beta;
        exponential_distribution<double> exp_rho;
        exponential_distribution<double> exp_death;
        exponential_distribution<double> exp_betaEnvironment;
        uniform_real_distribution<double> unif_real;
        uniform_int_distribution<int> unif_int;
        mt19937 rng;
    
    

        //containers to keep track of various pieces of information
        vector<Person*> people;
        array<double,1>previousTime;
        priority_queue<Event, vector<Event>, compTime > EventQ;
        array<double,1>finalTime;
        vector<double> timeOfParalyticCase;
    
        double Now;                 // Current "time" in simulation
        double Environment;
        double inactivationRate = .0304*365; //assumes anaerobic, nonsteril soil avg temp 23 deg C (Hurst paper)**units are in log_10
    
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
                p->setTiterLevel(100.0);
              //  death(p); //set when each individual will die
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
            cout<<"in infect\n";
            p->setInfectionStatus('I');
            p->setTimeSinceInfection(0.0);
            p->setTiterLevel(std::min(2048.0,11.0*p->getTiterLevel()));//boost 10 fold
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
            //time to waning
           // exponential_distribution<double> exp_rho(RHO);//**temp waning rate
           // double Tw = exp_rho(rng) + Now;
           // EventQ.push(Event(Tw,"inf_wane",p));
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
                EventQ.push(Event(Tc,"inf_c",p));
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,"vacc_r",p));
           // exponential_distribution<double> exp_rho(RHO);//**temp waning rate
           // double Tw=exp_rho(rng) + Now;
           // EventQ.push(Event(Tw,"inf_wane",p));
            return;

        }
    
        void infectByEnvironment(Person* p){
            assert(p->getInfectionStatus()=='S');
            cout<<"in infectByEnvironment\n";
            p->setInfectionStatus('I');
            p->setTimeSinceInfection(0.0);
            p->setTiterLevel(std::min(2048.0,11.0*p->getTiterLevel()));//boost 10 fold
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
            //time to waning
            // exponential_distribution<double> exp_rho(RHO);//**temp waning rate
            // double Tw = exp_rho(rng) + Now;
            // EventQ.push(Event(Tw,"inf_wane",p));
            return;
        }
      /*  void death(Person* p){
            //time to death
            exponential_distribution<double> exp_death(DEATH);
            double Td = exp_death(rng) + Now;
            Td > 100.0 ? Td = 100.0: Td = Td; //max lifespan is 100
            EventQ.push(Event(Td,"death",p));
            cout<<"death time: "<<Td<<"\n";
            return;
        }*/
    
        int nextEvent() {
            if(EventQ.empty()) return 0;
            Event event = EventQ.top();
            Now = event.time;
            for(Person* p: people) {//update individual's demography before event (do we want to make aging an event?)
                double r4 = unif_real(rng);
                if(r4 < DEATH or (p->getAge()>=100)){///temporary death***need to change
                   *p = Person();//keeps population constant
                }
                p->setTimeSinceInfection(Now-previousTime[0]);
                p->updateAge((Now-previousTime[0]));
                p->waning();//should this come before or after event occurs?
                // does this person contact the environment?
                exponential_distribution<double> exp_betaEnvironment(DEATH);//**temporary contact rate
                double rand2 = unif_real(rng);
                double Tce = exp_betaEnvironment(rng)+Now;
                if(rand2 < exp_betaEnvironment(rng)){
                    EventQ.push(Event(Tce,"env_c",p));
                    break;
                }
                double rand = unif_real(rng);
                if(rand<p->shedding(p->getTimeSinceInfection()) and p->getTiterLevel()!=0){
                    Environment+=128*(365*(Now-previousTime[0]))*p->stoolViralLoad(p->getTimeSinceInfection()*365);//avg 46720 g of feces per year (128 g/day)
                    cout<<"Now - previous time: "<<Now-previousTime[0]<<"\n";
                    cout<<"Environment: "<<128*(365*(Now-previousTime[0]))*p->stoolViralLoad(p->getTimeSinceInfection()*365)<<"\n";
                }
            }
            //first update Environment
            double rand1=unif_real(rng);
            if(rand1<delta and Environment>=0){
                Environment-=inactivationRate;
            }
            Person &individual = *event.people;
            if(event.type=="inf_c"){//includes contact with infected and vaccinated individual (OPV)
                double r1 = unif_real(rng);
                double r2 = unif_int(rng);
                if(r2<=totalSusceptibles()){//if there are enough susceptibles in the pop then a contact will lead to an infection
                    for(Person* p: people) {
                        if(r1<(p->probInfGivenDose(infDose))){//infects any susceptible whose antibody titer is low enough
                            infect(p);
                        }
                       
                    }
                }
                
            }
            else if(event.type=="vacc_c"){//contact with vaccinated individual
                double r1 = unif_real(rng);
                double r2 = unif_int(rng);
                if(r2<=totalSusceptibles()){
                    for(Person* p: people){
                        if(r1<p->probInfGivenDose(vaccDose)){
                            infect(p);
                        }
                        
                    }
                }
            }
            else if (event.type == "inf_r") {//recovery from vacc and WPV may be different
                individual.setInfectionStatus('S');
                TTE=Now;
            }
            else if(event.type=="vacc_r"){
                individual.setInfectionStatus('S');
                TTE=Now;
            }
          //  else if(event.type=="death"){
          //      individual=Person();//keeps population constant
               // death(&individual);//give each new person a death
         //   }
            else if(event.type=="env_c"){
                double r1=unif_real(rng);
               for(Person* p: people) {
                   if(r1 < .0000000005 && p->getInfectionStatus()=='S'){
          //          if(r1<(p->probInfGivenDose(infDose))){//infects any susceptible whose antibody titer is low enough
                        infectByEnvironment(p);
                    }
                    
                }
            }
          /*  else if(event.type=="inf_wane"){
                individual.waning();
                individual.setInfectionStatus('S');
            }
            else if(event.type=="vacc_wane"){
                individual.waning();
                individual.setInfectionStatus('S');
            }*/
            previousTime[0]=event.time;
            EventQ.pop();
            
            return 1;
        }
    
};
#endif
