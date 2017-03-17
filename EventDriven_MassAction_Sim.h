#ifndef EDMASSIM_H
#define EDMASSIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include <random>
#include <assert.h>
<<<<<<< HEAD
#include <array>
#include <fstream>
#include <string>
#include <algorithm>

=======
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45

using namespace std;

class Event {
    public:
        double time;
        char type;
<<<<<<< HEAD
        Event(const Event& o) {  time = o.time; type=o.type; }//*** this is a copy constructor
        Event(double t, char e) { time=t; type=e; } // *** this is a constructor
        Event& operator=(const Event& o) {
            if(this != &o){
            time = o.time;
            type=o.type;
            }
            return *this;
        } //*** this is copy assignment operator
=======
        Event(const Event& o) {  time = o.time; type=o.type; }
        Event(double t, char e) { time=t; type=e; }
        Event& operator=(const Event& o) { time = o.time; type=o.type; }
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
};

class compTime {
    public:
        bool operator() (const Event* lhs, const Event* rhs) const {
<<<<<<< HEAD
            return (lhs->time>rhs->time); //(->) is an assignment operator for pointers
        }

        bool operator() (const Event& lhs, const Event& rhs) const {
            return (lhs.time>rhs.time);//returns true if lhs is ordered before rhs
=======
            return (lhs->time>rhs->time);
        }

        bool operator() (const Event& lhs, const Event& rhs) const {
            return (lhs.time>rhs.time);
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
        }
};

class EventDriven_MassAction_Sim {
    public:
                                    // constructor
<<<<<<< HEAD
        EventDriven_MassAction_Sim( int n, double gamma, double beta, double kappa, double rho, double birth, double death): N(n), GAMMA(gamma),BETA(beta),KAPPA(kappa),RHO(rho),BIRTH(birth),DEATH(death),rng((random_device())()), unif_int(1,n-1) {
=======
        EventDriven_MassAction_Sim( int n, double gamma, double beta): N(n), rng((random_device())()), exp_gamma(gamma), exp_beta(beta), unif_int(1,n-1) {
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
            reset();
        }

        int N;                      // population size
<<<<<<< HEAD
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
=======
        exponential_distribution<double> exp_gamma;
        exponential_distribution<double> exp_beta;
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
        uniform_int_distribution<int> unif_int;

                                    // event queue
        priority_queue<Event, vector<Event>, compTime > EventQ;
<<<<<<< HEAD
        //***** priority_queue<Type of object stored in queue, container used to implement queue, comparison function used to determine the smallest element>

    
        vector<int> Compartments;   // S, I1, R, P, Ir compartments, with counts for each
        //**** i want to change vector to array since there will be a fixed number of compartments
    
        double Now;                 // Current "time" in simulation
        vector<double>finalTime;
=======
        vector<int> Compartments;   // S, I, R compartments, with counts for each
        double Now;                 // Current "time" in simulation
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45

        mt19937 rng;

        void run_simulation() {
<<<<<<< HEAD
            while (next_event()>=0) {//**changed while loop structure to find avg time to end of sim
                if(next_event()==0){
                    finalTime.push_back(Now);
                    break;
                }
=======
            while (next_event()) {
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
                continue;
            }
        }

        int epidemic_size() {
            return Compartments[2]; // Recovered class
        }
<<<<<<< HEAD
        double FinalTime(){
            return finalTime[0];
        }


        void reset() {
            Now = 0.0;

            Compartments.clear();
            Compartments.resize(5,0);//***** resize to 5 compartments with 0 individuals in them
            Compartments[0] = N; //**** puts entire population in S
=======

        int reset() {
            Now = 0.0;

            Compartments.clear();
            Compartments.resize(3,0);
            Compartments[0] = N;
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
            
            EventQ = priority_queue<Event, vector<Event>, compTime > ();
        }

        void rand_infect(int k) {   // randomly infect k people
            for (int i = 0; i < k; i++) {
<<<<<<< HEAD
                infectI1();
=======
                infect();
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
            }
            return;
        }

<<<<<<< HEAD
        void infectI1() {
            assert(Compartments[0] > 0);//*** aborts if there are no susceptibles
            Compartments[0]--;      // decrement susceptibles
            Compartments[1]++;      // increment I1
                                    // time to recovery
            exponential_distribution<double> exp_gamma(GAMMA*Compartments[1]);
            double Tr = exp_gamma(rng) + Now;
                                    // time to next contact
            exponential_distribution<double> exp_beta(BETA*(Compartments[0]/N)*(Compartments[1]+KAPPA*Compartments[4]));
=======
        void infect() {
            assert(Compartments[0] > 0);
            Compartments[0]--;      // decrement susceptibles
            Compartments[1]++;      // increment infecteds
                                    // time to recovery
            double Tr = exp_gamma(rng) + Now;
                                    // time to next contact
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
            double Tc = exp_beta(rng) + Now;
            while ( Tc < Tr ) {     // does contact occur before recovery?
                add_event(Tc, 'c'); // potential transmission event
                Tc += exp_beta(rng);
            }
            add_event(Tr, 'r' );
            return;
        }
<<<<<<< HEAD
    
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
=======
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45

        bool is_susceptible(int x) {
            if (Compartments[0] >= x) return true;
            else return false;
        }
<<<<<<< HEAD
        bool is_partsusc(int y){
            if(Compartments[3] >= y) return true;
            else return false;
        }

        int next_event() {
            if ( EventQ.empty() ){
                return 0;
            }
=======

        int next_event() {
            if ( EventQ.empty() ) return 0;
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
<<<<<<< HEAD
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
=======
            if (event.type == 'r') {    // recovery event
                Compartments[1]--;      // decrement Infected class
                Compartments[2]++;      // increment Recovered class
            } else {                    // event type must be 'c'
                // N-2 because person can't self-infect, and because randint includes endpoints
                int contact = unif_int(rng); // add 1 b/c there's no person 0
                if ( is_susceptible(contact) ) infect();
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
            }
            return 1;
        }

        void add_event( double time, char type) {
            EventQ.push( Event(time,type) );
            return;
        }

<<<<<<< HEAD

=======
>>>>>>> ec27e3c4b803ef378741285cd9ebebc7c8ee5a45
};
#endif
