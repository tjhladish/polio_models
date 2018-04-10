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



const double PIR                  = 0.001; //type 3 paralysis incidence rate
const double maxTiter             = 2048.0;
const double minTiter             = 1.0;
const double maxAge               = 70.0;// (max age in simulation so that mean of dist not zero (in death function)
const double meanAge              = 18;//Nigeria mean age
double TTE                        = 0;
const double vaccRate             = 0.01;
//thresholds used to determine when an individual recovers -- needed to end contact events
const double WPVrecThresh         = 0.17;//used mean shedding duration (for WPV) in prob shedding at t function .17
const double OPVrecThresh         = 0.24;//used mean shedding duration (for OPV)in prob shedding at t function .24
const double propToVacc           = 0.25;

//Environment contact parameters
const double inactivationRate     = 2*pow(10,0.223161);//TCID50 per day assumes surface water at 22 deg C (converted from PFU to TCID5)(Hurst,1989)
//const double feces = .128; //Liters/day
//const double urine = 1.4; //Liters/day
const double gramsFeces           = 128.0;//avg grams of feces produced per day
const double chkEnvRate           = 12.0;//once per month
//const double fiftyPerInf=pow(10,7);//units TCID50/L used saturating incidence function -- concentration of PV at which infection rate is 50% of max
//const double propVirusinWater     = 1e-9;

//probability of shedding parameters (these are all means from Famulare paper) - units in days
const double muOPV                = 30.3;
const double sigmaOPV             = 1.86;
const double muWPV                = 43.0;
const double sigmaWPV             = 1.69;
const double deltaShedding        = 1.16;

//peak shedding concentration (these are all means from Famulare paper)
const double Smax                 = 6.7; //units TCID50/gram
const double Smin                 = 4.3; //units TCID50/gram
const double tau                  = 10.0; //units months
const double newBorn              = 7.0/12.0;//age in months

//stool viral load parameters (these are all means from Famulare paper)
const double k                    = 0.056; //unitless
const double eta                  = 1.65; //unitless
const double nu                   = 0.17; //unitless
const double xsi                  = 0.32; //unitless

//probability of infection given dose parameters (these are all means from Famulare paper)
const double betaDose             = 14.0; //TCID50
const double alphaDose            = 0.44; // unitless shape parameter
const double gammaDose            = 0.55; //unitless immunity-dependent shape parameter
const double infDose              = 1000.0; //units TCID50 per infection (upper bound on number of virus particles need to initiate infection)
const double vaccDose             = 1e5; //units TCID50 per vaccine
double envDose;//changes based on concentration of virus in environment (water source)

//waning parameters (these are means from Famulare paper)
const double waningLambda         = 0.75;//unitless

//environmental detection parameters (4/3/2017 - these are just guesses)
const double detectionRate        = 1.0/5.0; //can detect 1 virus particle in 5 liters of water?

//death rate parameters (indexmundi.com, Pakistan as of 2016)
const double infantMortalityRate  = 53.9/1000.0; //0-1 year olds
const double childMortalityRate   = 87.0/1000.0; // <5 year olds
const double adultMortalityRate   = 6.4/1000.0; //>5

//Teunis boosting/waning model
const double pathogenClearance    = 5.0;//needs to be fit
const double pathogenGrowth       = 4.0;//needs to be fit
const double antibodyGrowth       = 5.0;//needs to be fit

//shape parameter for waning fn
const double r                    = 2.19; //needs to be fit
const double antibodyDecay        = 1.49; //needs to be fit


#endif /* EventDriven_parameters_hpp */
