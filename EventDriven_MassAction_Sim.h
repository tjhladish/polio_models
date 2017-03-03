#ifndef EDMASSIM_H
#define EDMASSIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include <random>
#include <assert.h>

using namespace std;

class Event {
    public:
        double time;
        char type;
        Event(const Event& o) {  time = o.time; type=o.type; }
        Event(double t, char e) { time=t; type=e; }
        Event& operator=(const Event& o) { time = o.time; type=o.type; }
};

class compTime {
    public:
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
        EventDriven_MassAction_Sim( int n, double gamma, double beta): N(n), rng((random_device())()), exp_gamma(gamma), exp_beta(beta), unif_int(1,n-1) {
            reset();
        }

        int N;                      // population size
        exponential_distribution<double> exp_gamma;
        exponential_distribution<double> exp_beta;
        uniform_int_distribution<int> unif_int;

                                    // event queue
        priority_queue<Event, vector<Event>, compTime > EventQ;
        vector<int> Compartments;   // S, I, R compartments, with counts for each
        double Now;                 // Current "time" in simulation

        mt19937 rng;

        void run_simulation() {
            while (next_event()) {
                continue;
            }
        }

        int epidemic_size() {
            return Compartments[2]; // Recovered class
        }

        int reset() {
            Now = 0.0;

            Compartments.clear();
            Compartments.resize(3,0);
            Compartments[0] = N;
            
            EventQ = priority_queue<Event, vector<Event>, compTime > ();
        }

        void rand_infect(int k) {   // randomly infect k people
            for (int i = 0; i < k; i++) {
                infect();
            }
            return;
        }

        void infect() {
            assert(Compartments[0] > 0);
            Compartments[0]--;      // decrement susceptibles
            Compartments[1]++;      // increment infecteds
                                    // time to recovery
            double Tr = exp_gamma(rng) + Now;
                                    // time to next contact
            double Tc = exp_beta(rng) + Now;
            while ( Tc < Tr ) {     // does contact occur before recovery?
                add_event(Tc, 'c'); // potential transmission event
                Tc += exp_beta(rng);
            }
            add_event(Tr, 'r' );
            return;
        }

        bool is_susceptible(int x) {
            if (Compartments[0] >= x) return true;
            else return false;
        }

        int next_event() {
            if ( EventQ.empty() ) return 0;
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
            if (event.type == 'r') {    // recovery event
                Compartments[1]--;      // decrement Infected class
                Compartments[2]++;      // increment Recovered class
            } else {                    // event type must be 'c'
                // N-2 because person can't self-infect, and because randint includes endpoints
                int contact = unif_int(rng); // add 1 b/c there's no person 0
                if ( is_susceptible(contact) ) infect();
            }
            return 1;
        }

        void add_event( double time, char type) {
            EventQ.push( Event(time,type) );
            return;
        }

};
#endif
