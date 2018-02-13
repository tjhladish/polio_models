#include "Cholera_Sim.h"

int main() { 
    const double b       = 0.02*1/12;//birth rate
    const double beta    = 10.0;//transmission rate
    const double C       = 0.0002;//prob of symptomatic infection
    const double epsilon = 30.0/583;//waning rate from symptomatic infection
    const double gamma   = 30.0/14;//recovery rate from symptomatic infection
    const double mu      = 0.02*1/12;//death rate
    const double rho     = 30.0/102.2;//recovery from asymptomatic infection

    //const double N       = 790590.0/4;
    const double S       = 33509.20478;
    const double I       = 0.4518533614;
    const double Y       = 164119.6172;
    const double R       = 18.22614054;

    Cholera_Sim sim(b, beta, C, epsilon, gamma, mu, rho);
    sim.initialize(S, I, Y, R);

    const int max_time   = 168;

    for (int i = 0; i < max_time; ++i) {
        sim.printX();
        sim.step_simulation(1);
    }
    sim.printX();

    //sim.run_simulation();

    return 0;
}
