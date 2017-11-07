#include <cmath>
#include <iostream>
#include <string>

#include <gsl/gsl_integration.h>

// compiled on roc using:
// g++ -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic -I/home/tjhladish/work/AbcSmc/gsl_local/include/ integration_test2.cpp -o integ_test -lm -L/home/tjhladish/work/AbcSmc/gsl_local/lib/ -lgsl -lgslcblas -lpthread -ldl

using namespace std;

double integ_pi1(const double tau_max);
double integ_pi2(const double tau_max);
double integ_pi3(const double tau_max);
double integ_gamma1pi1(const double tau_max);
double integ_gamma2pi3(const double tau_max);
double R0(double tau, void *params);
double pi1(double tau, void *params);
double pi2(double tau, void *params);
double pi3(double tau,void *params);

struct f_params {
    double y0=1;
    double b0=1;
    double mu=0.318;
    double mu0=.16;
    double K=0.001;
    double delta=0.02;
    double contactRate1 =135.0/365;
    double contactRate2 = 100.0/365;
    double recover1=13.0/365;
    double recover2=31.0/365;
    double waneRate=0.02/365;
    double y1=10;
    double r=1.49;
    double nu=1.69;
    double t1=28;
    double c = b0*(mu-mu0)/(y0*(exp((mu-mu0)*t1)-1));
    double popSize = 10;
    double T = 28;//infectivity time
    double S = 120; //susceptiblility time

};

double beta1(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return p.contactRate1*(p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau))/(p.mu-p.mu0)))/((p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau))/(p.mu-p.mu0)))+p.K);
}
double gamma1(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return p.recover1*p.y0*exp(p.mu*tau)/(p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau)/(p.mu-p.mu0)))+p.K);
}
double gamma2(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return p.recover2*p.y0*exp(p.mu*tau)/(p.b0*exp(p.mu0*tau)-(p.c*p.y0*(exp(p.mu*tau)-exp(p.mu0*tau)/(p.mu-p.mu0)))+p.K);
}
double gamma1pi1(double tau,void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return gamma1(tau,&p)*pi1(tau,&p);
}
double gamma2pi3(double tau,void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return gamma2(tau,&p)*pi3(tau,&p);
}
double pi1(double tau, void *params) {
    f_params &p= *reinterpret_cast<f_params *>(params);
    //solved the integral below symbollically using mathematica
    return exp(-(p.K*p.recover1*tau+tau*p.delta+(p.recover1*(tau*p.mu0-log(-p.c*exp(tau*p.mu)*p.y0+exp(tau*p.mu0)*(p.c*p.y0+p.b0*(p.mu-p.mu0)))+log(p.b0*(p.mu-p.mu0))))/p.c));
}
double eq_I1(double tau, void *params){
    f_params &p = *reinterpret_cast<f_params *>(params);
    return p.mu*p.popSize*(1.0-1.0/R0(tau,&p))*pi1(tau,&p);
}
double pi2(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    //approximate the function using a Taylor series and then integrate that
    //integrated taylor series is summed below
    //double nTerm=0;
    double taylorSum=0;
    //double tol = 1e-2;
    //int n =0;//series counter

    int factorial;
    for(int n=0;n<5;n++){
        if(n==0){
            factorial=1;
        }
        else{
            factorial*=n;
        }

        taylorSum += pow(p.waneRate*p.mu,n)*pow(p.y1,n*p.r-1)*pow(p.nu,n)*pow(tau,n+1)/factorial;//this assumes concentration of antibodies is never zero
    }

    return exp(-taylorSum*(R0(tau,&p))+p.delta);
}
double eq_R(double tau,void *params){
    f_params &p = *reinterpret_cast<f_params *>(params);
    return (p.mu*p.popSize*(1.0-1.0/R0(tau,&p))*integ_gamma1pi1(p.T)+p.popSize*integ_gamma2pi3(p.T)*((1.0-1.0/R0(tau,&p))*(1.0-p.mu*(integ_pi1(p.T)+integ_gamma1pi1(p.T)*integ_pi2(p.S))))/integ_gamma2pi3(p.T)*integ_pi2(p.S)+integ_pi3(p.T))*pi2(tau,&p);
}
double pi3(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    //solved the integral below symbollically using mathematica
    return exp(-(p.K*p.recover2*tau+tau*p.delta+(p.recover2*(tau*p.mu0-log(-p.c*exp(tau*p.mu)*p.y0+exp(tau*p.mu0)*(p.c*p.y0+p.b0*(p.mu-p.mu0)))+log(p.b0*(p.mu-p.mu0))))/p.c));
}
double eq_Ir(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return p.popSize*(1.0-1.0/R0(tau,&p))*(1.0-p.mu*(integ_pi1(p.T)+integ_gamma1pi1(p.T)*integ_pi2(p.S)))/(integ_gamma2pi3(p.T)*integ_pi2(p.S)+integ_pi3(p.T))*pi3(tau,&p);
}

double R0(double tau, void *params){
    f_params &p= *reinterpret_cast<f_params *>(params);
    return beta1(tau,&p)*pi1(tau,&p);
}



/*double integ_eqI1(const double tau) {
    f_params params;


    gsl_function F;
    F.function = &eq_I1;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=1e-4;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}*/


/*double integ_eqR(const double tau) {
    f_params params;


    gsl_function F;
    F.function = &eq_R;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=1e-4;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}*/


double integ_eqIr(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &eq_Ir;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_gamma1pi1(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &gamma1pi1;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_gamma2pi3(const double tau) {
    f_params params;


    gsl_function F;
    F.function = &gamma2pi3;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_R0(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &R0;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-1;
    const double epsrel=1e-1;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_pi1(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &pi1;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-4;
    const double epsrel=0;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_pi2(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &pi2;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-1;
    const double epsrel=1e-1;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


double integ_pi3(const double tau) {
    f_params params;

    gsl_function F;
    F.function = &pi3;
    F.params = reinterpret_cast<void *>(&params);

    double result, error;
    size_t neval;

    const double xlow=0;
    const double xhigh=tau;
    const double epsabs=1e-1;
    const double epsrel=1e-1;

    gsl_integration_qng (&F,
            xlow,
            xhigh,
            epsabs,
            epsrel,
            &result,
            &error,
            &neval);
    return result;
}


int main(void) {
    const double step = 1.0;
    for (double tau_max = 0; tau_max < 28; tau_max += step){
        double result = integ_pi2(tau_max);
        //double result = integ_eqIr(tau);
        cout << "tau, Result " << tau_max << ", " << result << endl;

    }

}
