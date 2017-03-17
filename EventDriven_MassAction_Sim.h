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
#include <algorithm>



using namespace std;

class Event {
    public:
        double time;
        char type;

        Event(const Event& o) {  time = o.time; type=o.type; }//*** this is a copy constructor
        Event(double t, char e) { time=t; type=e; } // *** this is a constructor
        Event& operator=(const Event& o) {
            if(this != &o){
            time = o.time;
            type=o.type;
            }
            return *this;
        } //*** this is copy assignment operator
};

class compTime {
    public:
        bool operator() (const Event* lhs, const Event* rhs) const {

            return (lhs->time>rhs->time); //(->) is an assignment operator for pointers
        }

        bool operator() (const Event& lhs, const Event& rhs) const {
            return (lhs.time>rhs.time);
        }
};

class EventDriven_MassAction_Sim {
    public:
                                    // constructor
        EventDriven_MassAction_Sim( int n, double gamma, double beta, double kappa, double rho, double birth, double death): N(n), rng((random_device())()), GAMMA(gamma), BETA(beta), KAPPA(kappa), RHO(rho), BIRTH(birth), DEATH(death), unif_int(1,n-1) {
            reset();
        }

        int N;                      // population size
        double GAMMA;
        double BETA;
        double KAPPA;
        double RHO;
        double BIRTH;
        double DEATH;
    
       // exponential_distribution<double> exp_gamma;
       // exponential_distribution<double> exp_beta;
        /*exponential_distribution<double> exp_kappa;
        exponential_distribution<double> exp_rho;
        exponential_distribution<double> exp_birth
        exponential_distribution<double> exp_death;*/

        exponential_distribution<double> exp_gamma;
        exponential_distribution<double> exp_beta;
        uniform_int_distribution<int> unif_int;

                                    // event queue
        priority_queue<Event, vector<Event>, compTime > EventQ;
        //***** priority_queue<Type of object stored in queue, container used to implement queue, comparison function used to determine the smallest element>

    
        vector<int> Compartments;   // S, I1, R, P, Ir compartments, with counts for each
        //**** i want to change vector to array since there will be a fixed number of compartments
    
        double Now;                 // Current "time" in simulation
        vector<double>finalTime;


        mt19937 rng;

        void run_simulation() {
            while (next_event()>=0) {//**changed while loop structure to find avg time to end of sim
                if(next_event()==0){
                    finalTime.push_back(Now);
                    break;
                }
            }
        }

        int epidemic_size() {
            return Compartments[2]; // Recovered class
        }

        double FinalTime(){
            return finalTime[0];
        }


        void reset(){
            Now = 0.0;
            Compartments.clear();
            Compartments.resize(5,0);//***** resize to 5 compartments with 0 individuals in them
            Compartments[0] = N; //**** puts entire population in S
        }
       
        void rand_infect(int k) {   // randomly infect k people
            for (int i = 0; i < k; i++) {
                infectI1();
            }
            return;
        }


        void infectI1() {
            assert(Compartments[0] > 0);//*** aborts if there are no susceptibles
            Compartments[0]--;      // decrement susceptibles
            Compartments[1]++;      // increment I1
            // time to recovery
            exponential_distribution<double> exp_gamma(GAMMA*Compartments[1]);
            double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(BETA*(Compartments[0]/N)*(Compartments[1]+KAPPA*Compartments[4]));
            double Tc = exp_beta(rng) + Now;
            while ( Tc < Tr ) {     // does contact occur before recovery?
                add_event(Tc, 'r'); 
                Tc += exp_beta(rng);
            }
            add_event(Tr, 'c' );
            return;
        }
    
        void infectIr() {
            assert(Compartments[3] > 0);//*** aborts if there are no partially susceptibles
            Compartments[3]--;      // decrement susceptibles
            Compartments[4]++;      // increment I1
            // time to recovery
            exponential_distribution<double> exp_gamma((GAMMA/KAPPA)*Compartments[4]);
            double Tr = exp_gamma(rng) + Now;
            // time to next contact
            exponential_distribution<double> exp_beta(KAPPA*BETA*(Compartments[3]/N)*(Compartments[1]+KAPPA*Compartments[4]));
            double Tc = exp_beta(rng) + Now;
            while ( Tc < Tr ) {     // does contact occur before recovery?
                add_event(Tc, 'p'); // potential transmission event (p means P contact event)
                Tc += exp_beta(rng);
            }
            add_event(Tr, 'i' ); //***(i means Ir recovers)
            return;
        }


        bool is_susceptible(int x) {
            if (Compartments[0] >= x) return true;
            else return false;
        }

        bool is_partsusc(int y){
            if(Compartments[3] >= y) return true;
            else return false;
        }


        int next_event() {
            if ( EventQ.empty() ) return 0;
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
            if (event.type == 'r') {// I1 recovery event
                Compartments[1]--;      // decrement Infected class
                Compartments[2]++;      // increment Recovered class
            }
            if(event.type == 'i'){//***Ir recovery event
                Compartments[4]--;
                Compartments[2]++;
            }
            if(event.type=='p'){
                int pcontact = unif_int(rng);
                if(is_partsusc(pcontact)){
                    infectIr();
                }
            }
            else {                    // event type must be 'c'
                // N-2 because person can't self-infect, and because randint includes endpoints
                int contact = unif_int(rng); // add 1 b/c there's no person 0
                if ( is_susceptible(contact) ){
                infectI1();
                }
            }
            return 1;
        }

        void add_event( double time, char type) {
            EventQ.push( Event(time,type) );
            return;
        }

};
#endif
