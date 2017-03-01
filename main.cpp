//
//  main.cpp
//  Practice
//
//  Created by Celeste on 2/16/17.
//  Copyright © 2017 Celeste. All rights reserved.
//

#include <iostream>
#include <array>
#include <random>
#include <algorithm>
#include <math.h>


//create a class for each person in the simulation
class Person{
private:
    //attributes:
    int m_age;
    double m_timeSinceInfection;//infection means any contact with virus
    double m_titerLevel; //does this need to be a double?
    char m_infectionStatus; //'S','I','V','R'
public:
    //default constructor
    Person(int age=0, double timeSinceInfection=0.0,double titerLevel=1.0,char infectionStatus='S'):m_age(age),m_timeSinceInfection(timeSinceInfection),m_titerLevel(titerLevel),m_infectionStatus(infectionStatus){
        //defaults to age 0, 0 time since infection, 1 titer level, and naive susceptible
        
    }
    //set values for attributes
    void setAge(int age){
        m_age = age;
    }
    void setTimeSinceInfection(double timeSinceInfection){
        m_timeSinceInfection += timeSinceInfection;
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
    double waning(){
        if(m_infectionStatus=='R'){
            return std::max(1.0, 3*pow(m_titerLevel,-.75));
        }
        else{
            return 0.0;
        }
    }
    
    void print(int i){
        std::cout<<"Attributes for person: "<<(i+1)<<" : "<<m_age<<" , "<<m_timeSinceInfection<<" , "<<m_titerLevel<<" ,"<<m_infectionStatus<<"\n";
    }
    
};



int main(){
    //define number of individuals
    const int numPeople = 5;
    //define number of times to run sim
    const int numSims = 5;
    
    //create fixed array of number of class Person
    std::array<Person,numPeople> people;
    
    //generate a random real number
    std::random_device rd;//this is the seed
    std::mt19937 gen(rd());//this generates the rn with the above seed
    std::uniform_real_distribution<> dis(0.0, 1.0);//generates unif rn
    std::uniform_real_distribution<> dis1(1.0, 2048.0);//generates unif rn for titer level
    
    //randomly initialize population
    for (Person& p: people) { 
        //first sets age
        double rr=dis(gen);
        if(rr<.8){
            p.setAge(1);
        } else if (rr < .9) {
            p.setAge(2);
        } else {
            p.setAge(3);
        }
        //second sets time since last infection
        p.setTimeSinceInfection(dis(gen));
        
        //third sets titer level
        p.setTiterLevel(dis1(gen));

        //fourth sets infection status
    //    double r2 = dis(gen);
    //    if(r2<.1){
    //        p.setInfectionStatus('I');
    //    }
    //    else if(r2<.7){
    //        p.setInfectionStatus('R');
    //    }
    }
    //set one infected individual in the population
    people[1].setInfectionStatus('I'); // why not person 0? selecting the second person seems odd.
  
    //parameters
    //const int beta = 200;
    const int recoverInf = 13; //recover from WPV
    const int recoverVacc= 10; //recover from OPV
    // const double wane = .2;
    // const double age = .0025;
    // const double lambda = .75;
    // const double phi = .3; //vaccination rate
    double infDose = 5; //number of doses from WPV infection
    double vaccDose = 3; //number of doses from OPV/vacc
    
    //events
    for(int k=0;k<numSims;++k){
        //resets sums after each simulation
        int sSum =0;
        int iSum=0;
        int vSum=0;
        int rSum=0;
        double time =0.0;
        double ageSum=0;
        double waneSum=0;
        double infectionSum=0;
        double vaccinationSum=0;
        for (Person& p: people) { 
            //first get total num in each class
            switch(p.getInfectionStatus()){
                case 'S':
                    sSum+=(1/p.getTiterLevel()); //lower titer level=more susc
                    break;
                case 'I':
                    iSum+=1;
                    break;
                case 'V':
                    vSum+=1;
                    break;
                case 'R':
                    rSum+=1;
                    break;
            }
            //get total aging rate
            ageSum += 1.0/p.getAge();
            //get total waning rate
            waneSum += p.waning();
            //get total infection rate
            infectionSum += p.probInfGivenDose(infDose);
            //get total vaccination rate
            if(p.getAge()==1){
                vaccinationSum+=p.probInfGivenDose(vaccDose);
            }
        }

        // get total WPV recovery rate
        double recoverInfSum = recoverInf*iSum;
        //get total OPV recovery rate
        double recoverVaccSum = recoverVacc*vSum;

        //sum all rates and choose the event that will occur
        double eventSum = ageSum + infectionSum + vaccinationSum + recoverInfSum + recoverVaccSum + waneSum;
        double r4 = dis(gen);

        //chooses WPV infection to occur
        if(r4 < (infectionSum/eventSum)){
            //now choose the individual to which the infection will occur
            std::cout<<"infection\n";
            double r5=dis(gen);
            double sum=0.0;
            for (Person& p: people) { 
                if(r5<(sum+(p.probInfGivenDose(infDose)/infectionSum))){
                    p.setInfectionStatus('I');
                    p.setTimeSinceInfection(0.0);
                    p.setTiterLevel(2048.0);//boosts to max titer level
                    break;
                }
                sum+=(p.probInfGivenDose(infDose)/infectionSum);
            }
        }
        //chooses recovery from WPV to occur
        else if(r4<((recoverInfSum+infectionSum)/eventSum)){
            //now choose the individual that will recover
            std::cout<<"recover\n";
            for (Person& p: people) { 
                if(p.getTimeSinceInfection()>0.11 and p.getInfectionStatus()=='I'){
                    p.setInfectionStatus('R');
                    break;
                }
            }
        }
        //chooses waning to occur
        else if(r4<((waneSum+recoverInfSum+infectionSum)/eventSum)){
            //now choose the individual whose immunity wanes
            std::cout<<"wane\n";
            double r6=dis(gen);
            double sum=0.0;
            for (Person& p: people) { 
                if(r6<(sum+(p.waning()/waneSum))){
                    p.setTiterLevel(p.waning());
                    p.setInfectionStatus('S');
                    break;
                }
                sum+=(p.waning()/waneSum);
            }
        }
        //chooses vaccination to occur
        else if(r4<((vaccinationSum+waneSum+recoverInfSum+infectionSum)/eventSum)){
            std::cout<<"vacc\n";
            //now choose the individual to which the vaccination occurs
            double r8 = dis(gen);
            double sum =0.0;
            for (Person& p: people) { 
                if(r8<(sum+(p.probInfGivenDose(vaccDose)/vaccinationSum))){
                    p.setInfectionStatus('V');
                    p.setTimeSinceInfection(0.0);
                    p.setTiterLevel(2048.0); //set to max titer level
                    break;
                }
                sum+=(p.probInfGivenDose(vaccDose)/vaccinationSum);
            }
        }
        //chooses recover from vaccination to occur
        else if(r4<((recoverVaccSum+vaccinationSum+waneSum+recoverInfSum+infectionSum)/eventSum)){
            std::cout<<"recover from vacc\n";
            //now choose individual that recovers from OPV
            //double r9 = dis(gen);
            for (Person& p: people) { //.11 is 40 days
                if(p.getTimeSinceInfection()>.11 and p.getInfectionStatus()=='V'){
                    p.setInfectionStatus('R');
                    break;
                }

            }
        }
        //chooses aging to occur
        else{
            //now choose which individual will age
            std::cout<<"age\n";
            double r7 = dis(gen);
            double sum = 0.0;
            for (Person& p: people) {
                if(r7<(sum+(p.getAge()/ageSum))){
                    p.setAge((p.getAge()+1));
                    break;
                }
                sum+=(p.getAge()/ageSum);
            }
        }
        //generate exp rn
        std::exponential_distribution<>rng((eventSum));
        double t1 = rng(gen);

        //update all individual's time since infection
        for (Person& p: people) {
            p.setTimeSinceInfection(t1);
        }
        //update simulation time
        time+=t1;
        std::cout<<time<<"\n";
        for (unsigned int i = 0; i < people.size(); ++i) {
            people[i].print(i);
        }
    }

    return 0;
}
