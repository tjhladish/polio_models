//
//  EventDriven_Sim_nonexp_recovery_waningfn_ENV.hpp
//  RtoC++practice
//
//  Created by Celeste on 5/3/17.
//  Copyright © 2017 Celeste. All rights reserved.
//

#ifndef EDMASSIM_H
#define EDMASSIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <queue>
#include <random>
#include <assert.h>
#include <array>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <limits>
#include <numeric>
#include "EventDriven_parameters.hpp"


using namespace std;

enum EventType {INFECTIOUS_CONTACT, INFECTION_RECOVERY, DEATH, BEGIN_SHEDDING, CHECK_ENVIRONMENT, ENVIRONMENT_CONTACT, AGING, SHED_INTO_ENVIRONMENT,NUM_OF_EVENT_TYPES};

enum InfectionStatus{IDC,IE,NA};
//infection status means the most recent cause of infection
//NA means never infected
//DC suffix -> infection by direct contact
//E suffix -> infection by contacting environment

enum AgeClass{AGE5,AGE15,AGE100}; //discretize age class
//Age5 = 0-5
//Age15 = 5-15
//Age100 = 15+

class Person{
    friend class Event;
private:
    //Famulare specific attributes
    double m_initialTiterLevel;//need initial titer level after infection for waning fn
    double m_titerLevel;
    double m_timeAtInfection;
    double m_timeToRec=0.0;
    double m_timetoDeath=0.0;
    double m_timetoBirth=0.0;


    //General attributes
    InfectionStatus m_infectionStatus;
    int m_numInfectionsDC; //DC means direct contact
    int m_numInfectionsEnv;
    const int m_index;
    double m_timeToShed;
    int m_age;
    AgeClass m_ageClass;
    double m_previousAgingTime = 0.0;

    //Teunis specific attributes
    double m_durationInfection;
    double m_initialPathogen;
    double m_initialAntibody;
    double m_peakAntibody;
    double m_antibodyLevel;
    double m_pathogenLevel;

public:
    //default constructor
    Person(int idx, double initialTiterLevel = 1.0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), InfectionStatus infStat = NA, int numInf=0, int numInfEnv=0, double timeToShed = 0.0, int age =0.0, AgeClass a = AGE5, double durInf=0, double intPat = 0.0, double intAnt=1.0, double pkAnt=1.0, double antLvl = 1.0, double patLvl = 0.0):m_initialTiterLevel(initialTiterLevel), m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfectionsDC(numInf),m_numInfectionsEnv(numInfEnv),m_index(idx), m_timeToShed(timeToShed), m_age(age), m_ageClass(a), m_durationInfection(durInf), m_initialPathogen(intPat), m_initialAntibody(intAnt), m_peakAntibody(pkAnt), m_antibodyLevel(antLvl), m_pathogenLevel(patLvl){

    }
    //General Functions
    int getIndex () const {
        return m_index;
    }
    int getAge(){
        return m_age;
    }
    void setAge(int a){
        m_age = a;
    }
    int getAgeClass(){
        return m_ageClass;
    }
    void setAgeClass(AgeClass a){
        m_ageClass = a;
    }
    void setPreviousAgingTime(double time){
        m_previousAgingTime = time;
    }
    double getPreviousAgingTime(){
        return m_previousAgingTime;
    }
    int getNumInfectionsDC(){
        return m_numInfectionsDC;
    }
    int getNumInfectionsEnv(){
        return m_numInfectionsEnv;
    }
    void setInfectionStatus(InfectionStatus infectionStatus){
        m_infectionStatus = infectionStatus;
    }
    void setNumInfectionsDC(int numInf){
        m_numInfectionsDC = numInf;
    }
    void setNumInfectionsEnv(int numInfEnv){
        m_numInfectionsEnv = numInfEnv;
    }
    InfectionStatus getInfectionStatus(){
        return m_infectionStatus;
    }
    double getTimeToShed(){
        return m_timeToShed;
    }
    void setTimeToShed(double time){
        m_timeToShed = time;
    }
    void setRecoveryTime(double rtime){
        m_timeToRec = rtime;
    }
    void setDeathTime(double dtime){
        m_timetoDeath = dtime;
    }
    void setBirthTime(double btime){
        m_timetoBirth = btime;
    }
    double getBirthtime(){
        return m_timetoBirth;
    }
    double getRecoveryTime(){
        return m_timeToRec;
    }
    double getDeathTime(){
        return m_timetoDeath;
    }
    void reset(int waningImmunityScenario){
        setNumInfectionsDC(0);
        setNumInfectionsEnv(0);
        setInfectionStatus(NA);
        setAge(0);
        setAgeClass(AGE5);
        setTimeToShed(0.0);
        if(waningImmunityScenario==1){
            setInitialTiterLevel(minTiter);
            setTimeAtInfection(numeric_limits<double>::max());
        }
        else{
            setInitialAntibody(1.0);
            setInitialPathogen(0.0);
            setDurationInfection();//sets peak antibody level
        }
        
    }
    double convertToDays(double t){
        return t*365;
    }
    double convertToMonths(double t){
        return t*12;
    }


    //Famulare specific functions
    void setTiterLevel(double titerLevel){
        m_titerLevel = std::min(maxTiter,titerLevel);
    }
    void setInitialTiterLevel(double titerLevel){
        m_initialTiterLevel = std::min(maxTiter, titerLevel);
        m_titerLevel = m_initialTiterLevel; //set equal for the purposes of further calculation
    }
    void setTimeAtInfection(double t){
        m_timeAtInfection = t;
    }
    double getTiterLevel(){
        return m_titerLevel;
    }
    double getInitialTiterLevel(){
        return m_initialTiterLevel;
    }
    double getTimeAtInfection(){
        return m_timeAtInfection;
    }
    void waningFamulare(double t){
        const double tnew = convertToMonths(t-m_timeAtInfection); //time needs to be in months post infection
        if(tnew>=1){//only wanes after one month post infection
            m_titerLevel= max(minTiter, m_initialTiterLevel*pow((tnew),-waningLambda));
        }
        return;
    }
    double probInfGivenDoseFamulare(InfectionStatus infstat){//dose will vary based on cause of infection, used mean of all shape parameters
        if(infstat == IDC){
            return 1-pow((1.0+(infDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else if(infstat == IE){
            return 1-pow((1.0+(envDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else{
            return 0;
        }
    }
    double probShedding(double t){
        const double tnew = convertToDays(t - m_timeAtInfection);//***t needs to be time since infection (days)
        return .5*erfc((log(tnew)-(log(muWPV)-log(deltaShedding)*log(m_titerLevel)))/(sqrt(2.0)*log(sigmaWPV)));
    }
    double peakShedding(){//age needs to be converted to months
        if(m_age<newBorn){
            return Smax;
        }
        else{
            return ((Smax-Smin)*exp((convertToMonths(newBorn)-(convertToMonths(m_age)))/tau)+Smin);
        }
    }

    double stoolViralLoad(double t){
        const double tnew = convertToDays(t-m_timeAtInfection);
        return max(pow(10.0,2.6),pow(10.0,((1.0-k*log2(m_titerLevel))*log(peakShedding())))*(exp(eta-(pow(nu,2.0)/2.0)-(pow(log(tnew)-eta,2.0)/(2.0*pow(nu+xsi*log(tnew),2.0))))/tnew));
            //**units are in TCID50/g
    }


    //Teunis specific functions
    double getDurationInfection(){
        return m_durationInfection;
    }
    void setDurationInfection(){
        //sets duration of infection, peak antibody level
        m_durationInfection = (1/(antibodyGrowth - pathogenGrowth))*log(1+((antibodyGrowth - pathogenGrowth)*m_initialPathogen)/(double)(pathogenClearance*m_initialAntibody));
        m_peakAntibody = m_initialAntibody*exp(antibodyGrowth*m_durationInfection);
    }
    double getPeakAntibody(){
        return m_peakAntibody;
    }
    void setInitialPathogen(double pat){
        m_initialPathogen = pat;
    }
    double timePeakPathogen(){
        //assume this is the time individual sheds into environment
        return (1/(antibodyGrowth - pathogenGrowth))*log((pathogenGrowth/antibodyGrowth)*(1+((antibodyGrowth-pathogenGrowth)*m_initialPathogen)/(pathogenClearance*m_initialAntibody)));
    }
    double getInitialAntibody(){
        return m_initialAntibody;
    }
    void setInitialAntibody(double ant){
        m_initialAntibody = ant;
    }
    double getInitialPathogen(){
        return m_initialPathogen;
    }
    double getAntibodyLevel(double t){
        //will this function name be confusing if we use it for both waning and boosting?
        const double tnew = convertToDays(t);
        if(t<=m_durationInfection){
            m_antibodyLevel = m_initialAntibody*exp(antibodyGrowth*t);
        }
        else{
            m_antibodyLevel = m_peakAntibody*pow((1+(r - 1)*pow(m_peakAntibody,(r-1))*antibodyDecay*(tnew-m_durationInfection)),(-1/(r-1)));
        }
        return m_antibodyLevel;
    }
    double getPathogenLevel(double t){
        const double tnew = convertToDays(t);
        if(tnew <= m_durationInfection){
            m_pathogenLevel = m_initialPathogen*exp(pathogenGrowth*tnew) - (pathogenClearance*m_initialAntibody*(exp(antibodyGrowth*tnew)-exp(pathogenGrowth*tnew))/(antibodyGrowth-pathogenGrowth));
        }
        else{
            m_pathogenLevel = 0.0;
        }
        return m_pathogenLevel;
    }
    /*void waningTeunis(double t){
        const double tnew = convertToDays(t);
        m_antibodyLevel = m_peakAntibody*pow((1+ (r - 1)*pow(m_peakAntibody,(r-1))*antibodyDecay*(tnew-m_durationInfection)),(-1/(r-1)));
    }*/

};

class Event {
    friend class Person;
public:
    double time;
    EventType type;
    Person* person;

    Event(double t, EventType e, Person* p):time(t),type(e),person(p){}

};

class compTime {
public:
    //two in case we use pointers at some point
    bool operator() (const Event* lhs, const Event* rhs) const {

        return (lhs->time>rhs->time);
    }

    bool operator() (const Event& lhs, const Event& rhs) const {
        return (lhs.time>rhs.time);
    }
};


class EventDriven_MassAction_Sim {
public:
    // constructor
    EventDriven_MassAction_Sim(const int n, const double beta, const double death, const double maxRunTime, const int wanIm, const double virIntWat,int seed=(random_device())()): rng(seed), N(n), BETA(beta), DEATH_RATE(death), maxRunTime(maxRunTime), waningImmunityScenario(wanIm), propVirusinWater(virIntWat),unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
        people = vector<Person*>(n);
        for (Person* &p: people) p = new Person(people_counter++);
        Now=0.0;
        ii=0;//counter for displaying events in queue
        IDCvec={0};
        IEvec={0};
        NonInfvec={0};
        timevec={0};
        antTit50={0};
        antTit100={0};
        antTit2048={0};
        ageDist = {0};
        event_counter = vector<int>(NUM_OF_EVENT_TYPES, 0);
        numBirths=0;
        numDeaths=0;
        numDCInf=0;
        numEInf=0;
        delta=0;
        NonInfsum=n;
        IDCsum=0;
        IEsum=0;
        sheddingPeople.clear();
    }

    ~EventDriven_MassAction_Sim() {
        for (Person* &p: people) delete p;
    }


    mt19937 rng;
    const int N; //total population size
    const double BETA;
    const double DEATH_RATE;
    double virusCon;
    double environment;
    const double maxRunTime;
    const int waningImmunityScenario;
    const double propVirusinWater;
    

    exponential_distribution<double> exp_beta;
    exponential_distribution<double> exp_death;
    uniform_real_distribution<double> unif_real;
    uniform_int_distribution<int> unif_int;

    //containers to keep track of various pieces of information
    vector<Person*> people;
    priority_queue<Event, vector<Event>, compTime> EventQ;
    unordered_set<Person*> sheddingPeople;
    vector<int> ageAtFirstInfect;


    int people_counter;
    const double reportingThreshold = 0.01;//used to determine how often to count compartments
    double Now; //simulation time
    double previousTime; //used for time step calculation

    //used for displaying events in queue
    double counter=0.0;
    int ii;
    vector<int> event_counter;

    //used to keep track of number in each compartment -- to be relaxed in next code iteration
    int NonInfsum; //all individuals not infected
    int IDCsum;
    int IEsum;

    //vectors used for export of compartment counts to csv at end of sim
    vector<double> IDCvec;
    vector<double> IEvec;
    vector<double> NonInfvec;
    vector<double> timevec;
    vector<double> antTit50;
    vector<double> antTit100;
    vector<double> antTit2048;
    vector<int> ageDist;


    //used to keep track of number of events that occur
    int numBirths;
    int numDeaths;
    int numDCInf;
    int numEInf;

    //used to determine how often to add compartment counts to above vectors
    double delta;
    
    double ageSum;
    double meanAge;//used for determining death times


    void runSimulation(){
        IDCvec[0] = IDCsum;
        IEvec[0] = IEsum;
        NonInfvec[0] = NonInfsum;
        timevec[0]=0;

        while(nextEvent() and EventQ.top().time < maxRunTime) {
        }
        if(EventQ.top().time >= maxRunTime){
            int Age_00_04   = 0;
            int Age_05_09   = 0;
            int Age_10_14   = 0;
            int Age_15_19   = 0;
            int Age_20_24   = 0;
            int Age_25_29   = 0;
            int Age_30_34   = 0;
            int Age_35_39   = 0;
            int Age_40_44   = 0;
            int Age_45_49   = 0;
            int Age_50_54   = 0;
            int Age_55_59   = 0;
            int Age_60_64   = 0;
            int Age_65_69   = 0;
            int Age_70_74   = 0;
            int Age_75_79   = 0;
            int Age_80_84   = 0;
            int Age_85_89   = 0;
            int Age_90_94   = 0;
            int Age_95_99   = 0;
            for(Person* p: people){
                if(p->getAge()<=4){
                    Age_00_04++;
                }
                else if(p->getAge()<=9){
                    Age_05_09++;
                }
                else if(p->getAge()<=14){
                    Age_10_14++;
                }
                else if(p->getAge()<=19){
                    Age_15_19++;
                }
                else if(p->getAge()<=24){
                    Age_20_24++;
                }
                else if(p->getAge()<=29){
                    Age_25_29++;
                }
                else if(p->getAge()<=34){
                    Age_30_34++;
                }
                else if(p->getAge()<=39){
                    Age_35_39++;
                }
                else if(p->getAge()<=44){
                    Age_40_44++;
                }
                else if(p->getAge()<=49){
                    Age_45_49++;
                }
                else if(p->getAge()<=54){
                    Age_50_54++;
                }
                else if(p->getAge()<=59){
                    Age_55_59++;
                }
                else if(p->getAge()<=64){
                    Age_60_64++;
                }
                else if(p->getAge()<=69){
                    Age_65_69++;
                }
                else if(p->getAge()<=74){
                    Age_70_74++;
                }
                else if(p->getAge()<=79){
                    Age_75_79++;
                }
                else if(p->getAge()<=84){
                    Age_80_84++;
                }
                else if(p->getAge()<=89){
                    Age_85_89++;
                }
                else if(p->getAge()<=94){
                    Age_90_94++;
                }
                else{
                    Age_95_99++;
                }
            }
            ageDist.resize(20);
            ageDist[0] = Age_00_04;
            ageDist[1] = Age_05_09;
            ageDist[2] = Age_10_14;
            ageDist[3] = Age_15_19;
            ageDist[4] = Age_20_24;
            ageDist[5] = Age_25_29;
            ageDist[6] = Age_30_34;
            ageDist[7] = Age_35_39;
            ageDist[8] = Age_40_44;
            ageDist[9] = Age_45_49;
            ageDist[10] = Age_50_54;
            ageDist[11] = Age_55_59;
            ageDist[12] = Age_60_64;
            ageDist[13] = Age_65_69;
            ageDist[14] = Age_70_74;
            ageDist[15] = Age_75_79;
            ageDist[16] = Age_80_84;
            ageDist[17] = Age_85_89;
            ageDist[18] = Age_90_94;
            ageDist[19] = Age_95_99;
        }
    }

    void randomizePopulation(int infected){
        int TiterLevel50=0;
        int TiterLevel100=0;
        int TiterLevel2048=0;
        int Age_00_04   = 0;
        int Age_05_09   = 0;
        int Age_10_14   = 0;
        int Age_15_19   = 0;
        int Age_20_24   = 0;
        int Age_25_29   = 0;
        int Age_30_34   = 0;
        int Age_35_39   = 0;
        int Age_40_44   = 0;
        int Age_45_49   = 0;
        int Age_50_54   = 0;
        int Age_55_59   = 0;
        int Age_60_64   = 0;
        int Age_65_69   = 0;
        int Age_70_74   = 0;
        int Age_75_79   = 0;
        int Age_80_84   = 0;
        int Age_85_89   = 0;
        int Age_90_94   = 0;
        int Age_95_99   = 0;

        for(Person* p: people) {

            //age distribution based on Pakistan data (pbs.gov.pk - last updated 4/18/18)
            //5 year age buckets
            //age distribution within age buckets taken from fitting age distribution to exponential function
            discrete_distribution<int> age {18.33,15.16,12.54,10.38,8.58,7.10,5.88,4.86,4.02,3.33,2.75,2.28,1.88,1.56,1.29,1.07,0.88,0.73,0.60,0.50};
            int personAge = age(rng);
            p->setAge(chooseAge(personAge,"age"));
            if(p->getAge()<=5){
                p->setAgeClass(AGE5);
            }
            else if(p->getAge()<=15){
                p->setAgeClass(AGE15);
            }
            else{
                p->setAgeClass(AGE100);
            }

            p->setInfectionStatus(NA);
            p->setInitialTiterLevel(maxTiter);
            p->setTimeAtInfection(numeric_limits<double>::max());
            p->setBirthTime(Now);//used to distinguish between old and new events
            deathTime(p);
            
            //double agingYear = unif_real(rng);
            //p->setAgeEventTime(Now);
            p->setPreviousAgingTime(Now);
            EventQ.emplace(Now+1.0,AGING,p);
            event_counter[AGING]++;

            //set environment contact--occurs daily
            //random so that at initialization everyone isn't contacting at exact same time
            double envCon = unif_real(rng)/(double)365;
            EventQ.emplace(envCon,ENVIRONMENT_CONTACT,p);
            event_counter[ENVIRONMENT_CONTACT]++;
            
            //beginning age distribution
            if(p->getAge()<=4){
                Age_00_04++;
            }
            else if(p->getAge()<=9){
                Age_05_09++;
            }
            else if(p->getAge()<=14){
                Age_10_14++;
            }
            else if(p->getAge()<=19){
                Age_15_19++;
            }
            else if(p->getAge()<=24){
                Age_20_24++;
            }
            else if(p->getAge()<=29){
                Age_25_29++;
            }
            else if(p->getAge()<=34){
                Age_30_34++;
            }
            else if(p->getAge()<=39){
                Age_35_39++;
            }
            else if(p->getAge()<=44){
                Age_40_44++;
            }
            else if(p->getAge()<=49){
                Age_45_49++;
            }
            else if(p->getAge()<=54){
                Age_50_54++;
            }
            else if(p->getAge()<=59){
                Age_55_59++;
            }
            else if(p->getAge()<=64){
                Age_60_64++;
            }
            else if(p->getAge()<=69){
                Age_65_69++;
            }
            else if(p->getAge()<=74){
                Age_70_74++;
            }
            else if(p->getAge()<=79){
                Age_75_79++;
            }
            else if(p->getAge()<=84){
                Age_80_84++;
            }
            else if(p->getAge()<=89){
                Age_85_89++;
            }
            else if(p->getAge()<=94){
                Age_90_94++;
            }
            else{
                Age_95_99++;
            }

            //antibody titer level distribution
            if(p->getTiterLevel()<=50){
                TiterLevel50++;
            }
            else if(p->getTiterLevel()<=100){
                TiterLevel100++;
            }
            else{
                TiterLevel2048++;
            }
        }
        
        for(int i = 0; i < infected; i++){
            if(waningImmunityScenario==1){//maybe it would be better to make this a better descriptor than 1 for Fam and 2 for Teunis
                NonInfsum--;
                IDCsum++;
                infectByDirectContactFamulare(people[i]);
            }
            else{
                NonInfsum--;
                IDCsum++;
                people[i]->setInitialPathogen(10);//input variable to change
                infectByDirectContactTeunis(people[i]);//infect first "infected" number of people in vector
            }
        }

        //Environmental Surveillance occurs monthly
        double chkEnv = unif_real(rng)/(double)12;
        EventQ.emplace(chkEnv,CHECK_ENVIRONMENT,nullptr);
        event_counter[CHECK_ENVIRONMENT]++;
        
        //meanAge = ageSum/people.size();

        antTit50[0]=TiterLevel50;
        antTit100[0]=TiterLevel100;
        antTit2048[0]=TiterLevel2048;
        cout<<"num 95+ "<<Age_95_99<<"\n";
    }

    //use these functions to retrieve end of simulation data
    const double maxVectorLength = maxRunTime/reportingThreshold;
    vector<double> printVectorIDC(){
        if(IDCvec.size() < maxVectorLength){
            IDCvec.resize(maxVectorLength,IDCvec.back());
        }
        return IDCvec;
    }
    vector<double> printVectorNonInf(){
        if(NonInfvec.size() < maxVectorLength){
            NonInfvec.resize(maxVectorLength,NonInfvec.back());
        }
        return NonInfvec;
    }
    vector<double> printVectorIE(){
        if(IEvec.size() < maxVectorLength){
            IEvec.resize(maxVectorLength,IEvec.back());
        }
        return IEvec;
    }
    vector<double> printTimeVector(){
        if(timevec.size() < maxVectorLength){
            timevec.resize(maxVectorLength,timevec.back());
        }
        return timevec;
    }
    vector<int> printAgeDist(){
        return ageDist;
    }
    vector<double> printAntTit50(){
        return antTit50;
    }
    vector<double> printAntTit100(){
        return antTit100;
    }
    vector<double> printAntTit2048(){
        return antTit2048;
    }
    vector<int> printAvgAgeAtFirstInfect(){
        return ageAtFirstInfect;
    }
    
    int NumBirths(){
        return numBirths;
    }
    int NumDeaths(){
        return numDeaths;
    }
    int NumDCInf(){
        return numDCInf;
    }
    int NumEInf(){
        return numEInf;
    }
    int numNonInf(){
        return NonInfsum;
    }
    int numInfDC(){
        return IDCsum;
    }
    int numInfE(){
        return IEsum;
    }
    
    int chooseAge(int personAge,string event){
        int age;
        vector<double> ageVec = {0,0,0,0,0};
        double ageVecSum = 0.0;
        int index = 0;
        //use person age to find the age group and put those weights into ageVec
        for(unsigned int i = 0; i < numAgeBuckets; i++){
            index = personAge*numAgeBuckets + i;
            if(event=="death"){
                ageVec[i] = age_Death[index];
            }
            else if(event=="age"){
                ageVec[i] = age_Aging[index];
            }
            ageVecSum +=ageVec[i];
        }
        //choose age withing age group
        double rand = unif_real(rng);
        double ageCumSum = 0.0;
        int ageGroupIndex = 0;
        for(unsigned int i = 0; i < ageVec.size(); i++){
            ageCumSum += ageVec[i];
            if(rand < (ageCumSum/ageVecSum)){
                ageGroupIndex = i;
            }
        }
        age = personAge*numAgeBuckets + ageGroupIndex;
        ageVec.clear();
        return age;
    }
    
    void infectByDirectContactFamulare(Person* p) {

        //update number of infections
        p->setNumInfectionsDC(p->getNumInfectionsDC()+1);
        numDCInf++;
        //find average age at first infection
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv())==1 && Now >= maxRunTime/2){
            ageAtFirstInfect.push_back(p->getAge());
        }
        p->setInfectionStatus(IDC);
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold

        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);//assume shedding begins a day later--can't be instantaneous
        event_counter[BEGIN_SHEDDING]++;
        p->setTimeToShed(sheddingTime);

        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;

        while (p->probShedding(Tc)>WPVrecThresh) {
            if(Tc > sheddingTime){//only make infectious contacts after shedding begins
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                event_counter[INFECTIOUS_CONTACT]++;
            }
            Tc += exp_beta(rng);
        }
        //time to recovery
        double Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        p->setRecoveryTime(Tr);
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;
        return;
    }
    void infectByDirectContactTeunis(Person* p){

        //update number of infections
        p->setNumInfectionsDC(p->getNumInfectionsDC()+1);
        numDCInf++;
        //find average age at first infection
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv())==1 && Now >= maxRunTime/2){
            ageAtFirstInfect.push_back(p->getAge());
        }
        p->setInfectionStatus(IDC);
        p->setDurationInfection();

        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);//assume shedding begins a day later
        event_counter[BEGIN_SHEDDING]++;
        p->setTimeToShed(sheddingTime);//need to prevent simultaneous infections

        //time to recovery
        double Tr = p->getDurationInfection() + Now;
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;

        //time to next human-human contact
        exponential_distribution<double> exp_beta(BETA);
        double Tc = exp_beta(rng) + Now;

        while (Tr>Tc) {
            EventQ.emplace(Tc, INFECTIOUS_CONTACT,p);
            event_counter[INFECTIOUS_CONTACT]++;
            Tc+=exp_beta(rng);
        }

        //time to shed into environment
        //EventQ.emplace(p->timePeakPathogen(),SHED_INTO_ENVIRONMENT,p);
        //event_counter[SHED_INTO_ENVIRONMENT]++;

    }

    void infectByEnvironment(Person* p) {

        //update number of infections
        p->setNumInfectionsEnv(p->getNumInfectionsEnv()+1);
        numEInf++;
        //find average age at first infection
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv())==1 && Now >= maxRunTime/2){
            ageAtFirstInfect.push_back(p->getAge());
        }
        p->setInfectionStatus(IE);
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold

        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);//assume shedding begins a day later--can't be instantaneous
        event_counter[BEGIN_SHEDDING]++;
        p->setTimeToShed(sheddingTime);

        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;

        while (p->probShedding(Tc)>WPVrecThresh) {
            if(Tc > sheddingTime){//only make infectious contacts after shedding begins
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                event_counter[INFECTIOUS_CONTACT]++;
            }
            Tc += exp_beta(rng);
        }
        //time to recovery
        double Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;
        return;
    }

    void environmentalSurveillance(){
        /*if(virusCon>detectionRate){
            cout<<"detected pathogen in water\n";
        }
        else{
            cout<<"did not detect pathogen in water\n";
        }*/

        //after each ES event generate time to next one--monthly
        EventQ.emplace(Now+(1/(double)12),CHECK_ENVIRONMENT,nullptr);
        event_counter[CHECK_ENVIRONMENT]++;
        return;
    }

    void deathTime(Person* p){
        
        vector<double> deathCDF {0.178,0.324,0.446,0.546,0.630,0.698,0.755,0.802,0.841,0.874,0.900,0.922,0.941,0.956,0.968,0.979,0.987,0.994,1.000};
        double rand = unif_real(rng);
        int deathAgeGroup;//do i want to put a default value just in case?
        for(unsigned int i = 0; i < deathCDF.size(); i++){
            if(rand < deathCDF[i]){
                deathAgeGroup = i;
                break;
            }
        }
        int deathAge = chooseAge(deathAgeGroup,"death");
        int deathTime = deathAge - p->getAge();
        assert(deathTime + p->getAge() <=99);
        double Td;
        if(deathTime < 0){
            Td = Now+.01;
        }
        else{
            Td = deathTime + Now;
        }
        EventQ.emplace(Td,DEATH,p);
        event_counter[DEATH]++;
        return;
    }
    
    void infectiousContact(Person* p){
        switch (waningImmunityScenario) {
            case 1://Famulare waning
            {
                //choose person to contact
                int contact_idx = unif_int(rng);
                if(contact_idx >= p->getIndex()) contact_idx++;
                Person* contact = people[contact_idx];
                
                //no simultaneous infections
                if(contact->getTimeToShed() < Now and sheddingPeople.count(contact)==0){
                    
                    //wane immunity if applicable
                    if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                        contact->waningFamulare(Now);
                    }
                    
                    //check susceptibility
                    double r2 = unif_real(rng);
                    if(r2<contact->probInfGivenDoseFamulare(IDC)){
                        NonInfsum--;
                        IDCsum++;
                        infectByDirectContactFamulare(contact);
                    }
                }
                break;
                
            }
            case 2://Teunis waning
            {
                
                //choose person to contact
                int contact_idx = unif_int(rng);
                if(contact_idx >= p->getIndex()) contact_idx++;
                Person* contact = people[contact_idx];
                
                //no simultaneous infections
                if(contact->getTimeToShed() < Now and sheddingPeople.count(contact)==0){
                    
                    //check susceptibility
                    double dose = p->getPathogenLevel(Now);//assumes all of individual's pathogen is transmitted
                    if((pathogenGrowth*dose - pathogenClearance*contact->getAntibodyLevel(Now)) > 0){
                        NonInfsum--;
                        IDCsum++;
                        infectByDirectContactTeunis(contact);
                    }
                }
                break;
            }
            default:
            {
                cout<<"incorrect input\n";
                exit(-1);
                break;
            }
        }
    }
    
    void infectionRecovery(Person* p){
        //no else statement to prevent dead person from fulfilling recovery time
        if(p->getInfectionStatus()==IDC){
            IDCsum--;
            NonInfsum++;
            sheddingPeople.erase(p);
            p->setInfectionStatus(NA);
        }
        else if(p->getInfectionStatus()==IE){
            IEsum--;
            NonInfsum++;
            sheddingPeople.erase(p);
            p->setInfectionStatus(NA);
        }
    }
    
    void environmentContact(Person* p){
        //define a dose of virus from environment
        virusCon = propVirusinWater*environment;
        
        //assumes all virus particles shed end up in water source -- will relax this assumption in another code iteration
        envDose = 2*virusCon; //2 is num L water drank per day
        
        //no simultaneous infections
        if(p->getTimeToShed() < Now and sheddingPeople.count(p)==0){
            
            //check susceptibility
            double r2 = unif_real(rng);
            if(r2<p->probInfGivenDoseFamulare(IE)){
                NonInfsum--;
                IEsum++;
                infectByEnvironment(p);
            }
        }
        
        //update environment with additions
        if(sheddingPeople.count(p) > 0){
            environment+=gramsFeces*p->stoolViralLoad(Now);
        }
        
        //set next time to contact environment
        EventQ.emplace(Now + (1/(double)365),ENVIRONMENT_CONTACT,p);
        event_counter[ENVIRONMENT_CONTACT]++;

    }
    
    void deathEvent(Person* p){
        numDeaths++;
        numBirths++;
        switch (p->getInfectionStatus()) {
            case NA:
            {
                //don't need to update compartment counts
                break;
            }
            case IDC:
            {
                IDCsum--;
                NonInfsum++;
                break;
            }
            case IE:
            {
                IEsum--;
                NonInfsum++;
                break;
            }
        }
        sheddingPeople.erase(p);
        p->reset(waningImmunityScenario);
        p->setBirthTime(Now);//this will be used to determine if events are old or new
        deathTime(p);
        p->setPreviousAgingTime(Now);
        agingTime(p);
    }
    
    void updateEnvironment(double timeStep){
        if(timeStep > 0 and environment > 0){
            environment-=exp(inactivationRate*timeStep*365);
            if(environment < 0){
                environment = 0;
            }
        }
    }
    
    void agingTime(Person* p){
        EventQ.emplace(Now + 1.0, AGING,p);
        event_counter[AGING]++;
    }


    int nextEvent() {
        Event event = EventQ.top();
        EventQ.pop();
        Now = event.time;

        //take compartment counts for output
        if(Now-delta>reportingThreshold){
            IDCvec.push_back(IDCsum);
            IEvec.push_back(IEsum);
            NonInfvec.push_back(NonInfsum);
            timevec.push_back(Now);
            int TiterLevel50 = 0;
            int TiterLevel100 = 0;
            int TiterLevel2048 = 0;
            for(Person* p: people){
                if(p->getTiterLevel()<=50){
                    TiterLevel50++;
                }
                else if(p->getTiterLevel()<=100){
                    TiterLevel100++;
                }
                else{
                    TiterLevel2048++;
                }
            }
            antTit50.push_back(TiterLevel50);
            antTit100.push_back(TiterLevel100);
            antTit2048.push_back(TiterLevel2048);
            delta=Now;
        }
        //display events in queue
        /*if(Now>counter){
          cout<<"Loop "<< ii<<"\n";
          cout<<"Now "<<Now<<"\n";
          cout<<" queue size "<<EventQ.size()<<"\n";
          cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT] << endl;
          cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY] << endl;
          cout<< "\tBEGIN_SHEDDING: " << event_counter[BEGIN_SHEDDING] << endl;
          cout<<"\tENVIRONMENT_CONTACT " << event_counter[ENVIRONMENT_CONTACT]<< endl;
          cout<<"\tCHECK_ENVIRONMENT " << event_counter[CHECK_ENVIRONMENT]<< endl;
          cout << "\tDEATH: " << event_counter[DEATH] << endl;
          cout << "\tAGING: "<< event_counter[AGING] << endl;
          counter+=.1;
          ii++;
          }*/
        //update virus particles in environment
        double timeStep = Now - previousTime;
        
        updateEnvironment(timeStep);

        Person* individual = event.person;
        if(event.type==INFECTIOUS_CONTACT and sheddingPeople.count(individual)>0){
            infectiousContact(individual);

        }
        else if (event.type==INFECTION_RECOVERY and sheddingPeople.count(individual)>0) {
            infectionRecovery(individual);
        }
        else if(event.type == BEGIN_SHEDDING and individual->getInfectionStatus()!=NA){
            sheddingPeople.insert(individual);
        }
        else if(event.type == AGING){
        if(individual->getPreviousAgingTime()<=Now-1.0){
            agingTime(individual);//sets next time to age
            individual->setAge(individual->getAge() + 1);
            individual->setPreviousAgingTime(Now);
        
            //sets age class - may be useful later
            if(individual->getAge()<=5){
                individual->setAgeClass(AGE5);
            }
            else if(individual->getAge()<=15){
                individual->setAgeClass(AGE15);
            }
            else{
                individual->setAgeClass(AGE100);
            }
        }
        }
        else if(event.type == ENVIRONMENT_CONTACT){
            environmentContact(individual);
        }
        else if(event.type==CHECK_ENVIRONMENT){
            environmentalSurveillance();
        }
        else if(event.type == DEATH){
            deathEvent(individual);
        }
        event_counter[event.type]--;
        previousTime = Now;
        return 1;
    }
};
#endif

