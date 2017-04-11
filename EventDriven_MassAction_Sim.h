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
#include <limits>
#include "EventDriven_parameters.hpp"




using namespace std;



class Person{
    friend class Event;
    private:
    //attributes
    double m_age;
    double m_titerLevel;
    double m_timeAtInfection;//infection means either OPV or WPV
    string m_infectionStatus; //'NS','I1','S','I','V'
    const int m_index;
    
    public:
    //default constructor
    Person(int idx, double age=0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), string infectionStatus =""):m_age(age),m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection),m_infectionStatus(infectionStatus),m_index(idx){
        
    }
    int getIndex () const {
        return m_index;
    }
    
    void setAge(double age){
        m_age = age;
    }
    void updateAge(double age){
        m_age +=age;
    }
    double deathTime(){
        return maxAge-m_age;
    }
    void reset(){
        m_age = 0.0;
        m_titerLevel = 1.0;
      //  m_timeSinceInfection=0.0;
        m_timeAtInfection=numeric_limits<double>::max();
        m_infectionStatus = "";
    }
    
 /*   void setTimeSinceInfection(double timeSinceInfection){
        m_timeSinceInfection = timeSinceInfection;
    }
    
    void updateTimeSinceInfection(double timeSinceInfection){
       // cout<<"in update time since infection loop\n";
       // cout<<"double time since infection "<< timeSinceInfection<<"\n";
       // cout<<"updated m_timeSinceInfection "<<m_timeSinceInfection+timeSinceInfection<<"\n";
        m_timeSinceInfection +=timeSinceInfection;
    }*/
    void setTiterLevel(double titerLevel){
        m_titerLevel = std::min(maxTiter,titerLevel);
    }
    void setTimeAtInfection(double t){
        m_timeAtInfection = t;
    }
    void setInfectionStatus(string infectionStatus){
        m_infectionStatus = infectionStatus;
    }
    double sheddingThreshold(string infectionStatus){
        if(infectionStatus=="I1" or infectionStatus=="I"){
            return WPVrecThresh;
        }
        else{
            return OPVrecThresh;
        }
    }
    int getAge(){
        return m_age;
    }
 /*   double getTimeSinceInfection(){
        return m_timeSinceInfection;
    }*/
    double getTiterLevel(){
        return m_titerLevel;
    }
    double getTimeAtInfection(){
        return m_timeAtInfection;
    }
    string getInfectionStatus(){
        return m_infectionStatus;
    }
    
    double probInfGivenDose(string infstat){//dose will vary based on if it is WPV or OPV
        if(infstat=="I1" or infstat=="I"){
            return 1-pow((1.0+(infDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));//put in mean and shape parameters
        }
        else if(infstat == "V"){//assumes m_titerLevel>0
            return 1-pow((1.0+(vaccDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));//put in mean and shape parameters
        }
        else{
            return 1.0;
        }
    }
    void waning(double t){
        if(t>=(1/(double)12)){//only wanes after one month post infection
            m_titerLevel= std::max(1.0, m_titerLevel*pow((t*12),-waningLambda));
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
    
    double probShedding(double t){//***t needs to be time individual was infected
        if(m_infectionStatus=="I" or m_infectionStatus=="I1"){
            return .5*erfc((log(t)-(log(muWPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaWPV));
        }
        else if(m_infectionStatus=="V"){
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
        //**units are in TCID50/g
    }
    
    void print(int i){
        std::cout<<"Attributes for person: "<<(i+1)<<" : "<<m_age<<" , "<<m_titerLevel<<" ,"<<m_infectionStatus<<" , "<<m_timeAtInfection<<"\n";
    }
};

class Event {
    friend class Person;
public:
    double time;
    string type;
    Person* person;

    
    // Event(const Event& o) {time = o.time; type=o.type;}//*** this is a copy constructor
    Event(double t=0.0, string e="r",Person* p = 0):time(t),type(e),person(p){} // *** this is a constructor
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
        EventDriven_MassAction_Sim(const int n, const double gamma, const double beta, const double kappa, const double rho, const double birth, const double death, const double betaenv): rng((random_device())()), GAMMA(gamma), BETA(beta), KAPPA(kappa), RHO(rho), BIRTH(birth), DEATH(death), BETAENV(betaenv),unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
            previousTime = {0};
            people = vector<Person*>(n);
            for (Person* &p: people) p = new Person(people_counter++);
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
        const double BETAENV;
    
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
        array<double,1>previousTime;
        array<double,1>finalTime;
        vector<double> timeOfParalyticCase;

    
        int people_counter; 
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
              //  p->setTimeSinceInfection(unif_real(rng));
                p->setTiterLevel(100.0);
                p->setTimeAtInfection(0.0);
                death(p); //set when each individual will die
            }
            for(int i=0;i<k;++i){
                infect(people[i]);//does it matter which individuals in array are initially infected??
            }
            //vaccinate(individual[1+j]);
            exponential_distribution<double> checkEnvironment(GAMMA);//**temp check rate--needs to change
            double Te = checkEnvironment(rng);
            EventQ.push(Event(Te,"check_env",nullptr));//nullptr since checking env does not involve a person
        }
        
        void printPeople(){
            for(unsigned int i = 0; i<people.size(); ++i) {
                people[i]->print(i);
            }
        }
    

        double FinalTime(){
            return TTE;
        }

      /*  int totalSusceptibles(){
            int sSum=0;
            for(Person* p: people) {
                if(p->getInfectionStatus()=="S"){
                    sSum+=1;
                }
            }
            return sSum;
        }*/
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
            cout<<"infection\n";
            numInfected++;
            if(p->getInfectionStatus()==""){
                p->setInfectionStatus("I1");
            }
            else{
                p->setInfectionStatus("I");
            }
            p->setTimeAtInfection(Now);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            //is this a paralytic case?
            double r1 = unif_real(rng);
            if(r1<PIR){
                timeOfParalyticCase.push_back(Now);
            }
            //time to recovery(non-active infection)
            //exponential_distribution<double> exp_gammaWPV(GAMMAWPV);
            //double Tr = exp_gamma(rng) + Now;
            // time to next human-human contact
            exponential_distribution<double> exp_beta(BETA);
            double Tc = exp_beta(rng) + Now;
            while (p->probShedding((Tc-p->getTimeAtInfection())*365)>WPVrecThresh) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,"inf_c",p));
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,"inf_r",p));//recovery used for decrementing # of infecteds
            return;
        }
    /*void infect(Person* p) {
        cout<<"infection\n";
        //    cout<<"age : "<<p->getAge()<<"\n";
        //    cout<<"titer level "<<p->getTiterLevel()<<"\n";
        //    cout<<"time since infection "<<p->getTimeSinceInfection()<<"\n";
        numInfected++;
        if(p->getInfectionStatus()==""){
            p->setInfectionStatus("I1");
        }
        else{
            p->setInfectionStatus("I");
        }
        // p->setTimeSinceInfection(0.0);
        p->setTimeAtInfection(Now);
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
    }*/
        void vaccinate(Person* p){
            cout<<"in vaccinate\n";
            p->setInfectionStatus("V");
            numInfected++;
            p->setTimeAtInfection(Now);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            //time to recovery
           // exponential_distribution<double> exp_gamma(GAMMA);//**temporary recovery rate
           // double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (p->probShedding((Tc-p->getTimeAtInfection())*365)>OPVrecThresh) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,"inf_c",p));
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,"vacc_r",p));
            return;

        }
  /*  void vaccinate(Person* p){
        p->setInfectionStatus("V");
        numInfected++;
        //    p->setTimeSinceInfection(0.0);
        p->setTimeAtInfection(Now);
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
        
    }*/
    
        void infectByEnvironment(Person* p){
            numInfected++;
            if(p->getInfectionStatus()==""){
                p->setInfectionStatus("I1");
            }
            else{
                p->setInfectionStatus("I");
            }
            p->setTimeAtInfection(Now);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            //is this a paralytic case?
            double r1 = unif_real(rng);
            if(r1<PIR){
                timeOfParalyticCase.push_back(Now);
                cout<<"paralytic case\n";
            }
            //time to recovery
           // exponential_distribution<double> exp_gamma(DEATH);//**temporary recovery rate
          //  double Tr = exp_gamma(rng) + Now;
            // time to next human-human contact
            exponential_distribution<double> exp_beta(BETA);//contact rate
            double Tc = exp_beta(rng) + Now;
            while (p->probShedding((Tc-p->getTimeAtInfection())*365)>WPVrecThresh) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,"inf_c",p));
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,"inf_r",p));
            cout<<"recovery\n";
            return;
        }
  /*  void infectByEnvironment(Person* p){
        numInfected++;
        if(p->getInfectionStatus()==""){
            p->setInfectionStatus("I1");
        }
        else{
            p->setInfectionStatus("I");
        }
        p->setTimeAtInfection(Now);
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
    }*/
        void death(Person* p){
            exponential_distribution<double> exp_death(DEATH);
            double deathAge = exp_death(rng);
            double Td = deathAge + Now;
            cout<<"Now "<<Now<<"\n";
            
            if((p->getAge()+deathAge)>100.0){
                Td = p->deathTime();
            }
            cout<<"death time " <<Td<<"\n";
            EventQ.push(Event(Td,"death",p));
            return;
        }
    
        void environmentalSurveillance(){
            if((Environment/(latrineVolume+Environment))>detectionRate){//Environment means virus particles
                cout<<"detected pathogen in water\n";
            }
            else{
                cout<<"did not detect pathogen in water\n";
            }
            return;
        }
    
        int nextEvent() {
            if(numInfected==0 or EventQ.empty()) return 0;
            Event event = EventQ.top();
            Now = event.time;
            //cout<<"Now "<<Now<<"\n";
            //cout<<"event "<<event.type<<"\n";
            double timeStep = Now - previousTime[0];
            minusEnvironment += inactivationRate*timeStep; //used only if environment contact event is selected
          /*  for(int i=0;i<people.size();i++){
                cout<<"person "<<(i+1)<<" titer level "<<people[i]->getTiterLevel()<<"\n";
            }*/
            latrineVolume+=(people.size()*feces*urine - evapRate)*timeStep*365;
            for(Person* p: people) {
                p->updateAge(timeStep);
                // does this person contact the environment?
                exponential_distribution<double> exp_betaEnvironment(BETAENV);//**temporary contact rate
                double rand2 = unif_real(rng);
                double TceStep = exp_betaEnvironment(rng);
                double Tce = TceStep+Now;
                if(rand2 < TceStep){
                    EventQ.push(Event(Tce,"env_c",p));
                }
           /*     double rand = unif_real(rng);
                if(rand<p->shedding(p->getTimeSinceInfection()) and p->getTiterLevel()!=0){
                    Environment+=128*(365*(Now-previousTime[0]))*p->stoolViralLoad(p->getTimeSinceInfection()*365);//avg 46720 g of feces per year (128 g/day)
                  //  cout<<"Now - previous time: "<<Now-previousTime[0]<<"\n";
                  //  cout<<"Environment: "<<128*(365*(Now-previousTime[0]))*p->stoolViralLoad(p->getTimeSinceInfection()*365)<<"\n";
                }*/
            }
            Person* individual = event.person;
            if(event.type=="inf_c"){//includes contact with infected and vaccinated individual (OPV)
                cout<<"in contact\n";
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
               // cout<<"contact index "<<contact_idx<<"\n";
               // cout<<"contact person "<<contact->getAge()<<"\n";
                double r1 = unif_real(rng);
                if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                contact->waning(Now - contact->getTimeAtInfection());
                }
                if(r1 < contact->probInfGivenDose(contact->getInfectionStatus())){
                    infect(contact);
                }

            }
            else if(event.type=="vacc_c"){//contact with vaccinated individual
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
                cout<<"contact index "<<contact_idx<<"\n";
                cout<<"contact person "<<people[contact_idx]<<"\n";
                double r1 = unif_real(rng);
                if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                contact->waning(Now - contact->getTimeAtInfection());
                }
                if(r1 < contact->probInfGivenDose(contact->getInfectionStatus())){
                    infect(contact);
                }
            }
            else if (event.type == "inf_r") {//recovery from vacc and WPV may be different
                //individual->setInfectionStatus("S");
                TTE=Now;
                numInfected--;
            }
            else if(event.type=="vacc_r"){
                //individual->setInfectionStatus("S");
                TTE=Now;
                numInfected--;
            }
            else if(event.type=="death"){
                cout<<"death\n";
                if(individual->probShedding(Now)<individual->sheddingThreshold(individual->getInfectionStatus())){
                    numInfected--;
                }
                individual->reset();//keeps population constant
                death(individual);//give each new person a death
            }
            else if(event.type=="env_c"){
                //first update Environment
                double r3 = unif_real(rng);
                for(Person* p: people){
                    if(p->getTiterLevel()!=0 and r3 < p->probShedding((Now - p->getTimeAtInfection())*365)){
                        Environment+=gramsFeces*p->stoolViralLoad((Now - p->getTimeAtInfection())*365);
                    }
                }
                Environment = std::max(0.0, Environment - minusEnvironment);
                double rand2 = unif_real(rng);
                if(rand2<individual->probInfGivenDose(individual->getInfectionStatus())){
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
