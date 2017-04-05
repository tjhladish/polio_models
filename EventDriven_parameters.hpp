//
//  EventDriven_parameters.hpp
//  EventDrivenModel
//
//  Created by Celeste on 4/1/17.
//  Copyright © 2017 Celeste. All rights reserved.
//

#ifndef EventDriven_parameters_hpp
#define EventDriven_parameters_hpp

#include <stdio.h>
#include <math.h>


const double infDose=5; //number of doses from WPV infection
const double vaccDose=3; //number of doses from OPV vacc
const double PIR = 0.001; //type 3 paralysis incidence rate
const double maxTiter = 2048.0;
const double maxAge=100.0;
double TTE =0;

//Environment parameters
double delta = 180/(double)365;//mean virus death rate in sand saturated with septic liquor (2003 - WHO Virus survival report)
//assume latrine has dimensions 10x1x1 in units m^3 and evaporation occurs at a rate of 12 mm/day
double evapRate = .12; //**units in L/day


//Environment contact parameters
double inactivationRate = (pow(10,.0304))*365; //assumes anaerobic, nonsterile sandy-loam soil avg temp 23 deg C (Hurst paper)
double latrineVolume = 0.0;
const double feces = .128; //Liters/day
const double urine = 1.4; //Liters/day
const double gramsFeces = 128;

//shedding parameters (these are all means from Famulare paper) - units in days
const double muOPV = 30.3;
const double sigmaOPV = 1.86;
const double muWPV = 43.0;
const double sigmaWPV = 1.69;
const double deltaShedding = 1.16;

//peak shedding concentration (these are all means from Famulare paper)
const double Smax = 6.7; //units TCID50/gram
const double Smin = 4.3; //units TCID50/gram
const double tau = 10/(double)24; //units years

//stool viral load parameters (these are all means from Famulare paper)
const double k = 0.056; //unitless
const double eta = 1.65; //unitless
const double nu = .17; //unitless
const double xsi = 0.32; //unitless

//probability of infection given dose parameters (these are all means from Famulare paper)
const double betaDose = 14.0; //TCID50
const double alphaDose = .44; //shape parameter
const double gammaDose = .55; //immunity-dependent shape parameter

//waning parameters (these are all means from Famulare paper)
const double Nab1 = 1000.0; //baseline immunity one month pose immunization?
const double waningLambda = .75;//unitless

//environmental detection parameters (4/3/2017 - these are just guesses)
double detectionRate = 1/(double)5; //can detect 1 virus particle in 5 liters of water?


#endif /* EventDriven_parameters_hpp */
