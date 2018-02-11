#ifndef DIFFEQ_SIM_H
#define DIFFEQ_SIM_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;


class DiffEq_Sim {
    private:
        double t;      //initial time
        double h;      //time step
        double tmax;   //max time
        double hmin;
        double erTol = 1e-10;

    public:
        unsigned int nbins;
        double* x;

        DiffEq_Sim() {
            t    = 0.0;      //initial time
            h    = 0.1;      //time step
            tmax = 2000;
            hmin = 0.2;
            nbins = 0;
        };

        DiffEq_Sim(unsigned int _nbins):nbins(_nbins) {
            t    = 0.0;      //initial time
            h    = 0.1;      //time step
            tmax = 100;
            hmin = 0.2;
        };

        ~DiffEq_Sim() {};

        void printX() { for(unsigned int i=0; i < nbins-1; i++) { cout << x[i] << " ";} cout << x[nbins-1] << endl; }

        vector<double> get_state() {
            vector<double> C;
            C.assign(x, x + nbins);
            return C;
        }

        inline double get_time() { return t; }


        virtual void initialize() {}
        virtual void derivative(const double x[], double dxdt[]){}

        static int function(double t, double const x[], double dxdt[], void *params) {
            DiffEq_Sim* model = static_cast <DiffEq_Sim*> (params);
            model->derivative(x, dxdt);
            return GSL_SUCCESS;
        }


        int run_simulation() {
            gsl_odeiv_evolve*  e = gsl_odeiv_evolve_alloc(nbins);
            gsl_odeiv_control* c = gsl_odeiv_control_y_new(erTol, 0);
            gsl_odeiv_step*    s = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, nbins);
            gsl_odeiv_system sys = {function, NULL, nbins, this };
            while (t < tmax) {  //convergence check here
                int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, tmax, &h, x);
                if (status != GSL_SUCCESS) { return status; }
            }
            return 0;
        }

        int step_simulation( double stepsize ) {
            gsl_odeiv_evolve*  e = gsl_odeiv_evolve_alloc(nbins);
            gsl_odeiv_control* c = gsl_odeiv_control_y_new(erTol, 0);
            gsl_odeiv_step*    s = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, nbins);
            gsl_odeiv_system sys = {function, NULL, nbins, this };

            double tstop = t+stepsize;
            while (t < tstop) {
                int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, tstop, &h, x);
                if (status != GSL_SUCCESS) { return status; }
            }
            return 0;
        }

        /*
           double* advance_simulation(double I_lim) {
           while (t < 2000) {  //convergence check here
           int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
           if (y[1] < I_lim) { return y; }
           }
           printY();
           return y;
           }
         */
};

#endif
