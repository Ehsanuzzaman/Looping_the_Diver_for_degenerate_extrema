#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <stdbool.h>
#include <stddef.h>
#include "diver.h"
#include "simpson_integration.h"


const int         nPar                 = 8;                            // Dimensionality of the parameter space
const double      lowerbounds[]        = {-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.};      // Lower boundaries of parameter space
const double      upperbounds[]        = { 1., 1.,1.,1.,1., 1.,1.,1.};      // Upper boundaries of parameter space
const char        path[]               = "example_c/output/example";   // Path to save samples, resume files, etc
const int         nDerived             = 0;                            // Number of derived quantities to output
const int         nDiscrete            = 0;                            // Number of parameters that are to be treated as discrete
const int         discrete[]           = {};                           // Indices of discrete parameters, Fortran style, i.e. starting at 1!!
const bool        partitionDiscrete    = false;                        // Split the population evenly amongst discrete parameters and evolve separately
const int         maxciv               = 1;                            // Maximum number of civilisations
const int         maxgen               = 100;                          // Maximum number of generations per civilisation
const int         NP                   = 160;                         // Population size (individuals per generation)
const int         nF                   = 1;                            // Size of the array indicating scale factors
const double      F[]                  = {0.6};                        // Scale factor(s).  Note that this must be entered as an array.
const double      Cr                   = 0.9;                          // Crossover factor
const double      lambda               = 0.8;                          // Mixing factor between best and rand/current
const bool        current              = false;                        // Use current vector for mutation
const bool        expon                = false;                        // Use exponential crossover
const int         bndry                = 3;                            // Boundary constraint: 1=brick wall, 2=random re-initialization, 3=reflection
const bool        jDE                  = true;                         // Use self-adaptive choices for rand/1/bin parameters as per Brest et al 2006
const bool        lambdajDE            = true;                         // Use self-adaptive rand-to-best/1/bin parameters; based on Brest et al 2006
const double      convthresh           = 1.e-6;                        // Threshold for gen-level convergence: smoothed fractional improvement in the mean population value
const int         convsteps            = 10;                           // Number of steps to smooth over when checking convergence
const bool        removeDuplicates     = true;                         // Weed out duplicate vectors within a single generation
const bool        doBayesian           = false;                        // Calculate approximate log evidence and posterior weightings
const double      maxNodePop           = 1.9;                          // Population at which node is partitioned in binary space partitioning for posterior
const double      Ztolerance           = 1.e-3;                        // Input tolerance in log-evidence
const int         savecount            = 100;                          // Save progress every savecount generations
const bool        resume               = false;                        // Restart from a previous run
const bool        outputSamples        = false;                        // Write output .raw and .sam (if nDerived != 0) files
const int         init_pop_strategy    = 0;                            // Initialisation strategy: 0=one shot, 1=n-shot, 2=n-shot with error if no valid vectors found.
const bool        discard_unfit_points = false;                        // Recalculate any trial vector whose fitness is above max_acceptable_value
const int         max_init_attempts    = 10000;                        // Maximum number of times to try to find a valid vector for each slot in the initial population.
const double      max_acceptable_val   = 1e6;                          // Maximum fitness to accept for the initial generation if init_population_strategy > 0, or any generation if discard_unfit_points = true.
//const int         seed                 = 1111;                      // base seed for random number generation; non-positive or absent means seed from the system clock
const int         verbose              = 1;                            // Output verbosity: 0=only error messages, 1=basic info, 2=civ-level info, 3+=population info


double pi = 3.1416;


//////////////////////////////////////////////////Your fuction should start from here\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
/////////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\






 double mN = 800.;
 double mZ2 = (1.6*2000.);
 double me = 0.0005 ;
 double mtau = 1.77686;
 double mHplus = 160.;
 double mHzero = 125.10;
 double mHpp = 400.;

/////////////////////////////////////////////////////////////
 double s_dfhzero = 1.6e-11;
 double r_dfhzero = 2e-4;


 double s_dfw2 = 0.16;
 double r_dfw2 = 6.25e-14;



 double s_dfhplus = 9.8e-12;
 double r_dfhplus = 25.;
 


 double s_dfhpp = 1.6e-12;
 double r_dfhpp = 1.97e-5;
 
 

 double s_dfz2 = 2.44e-14;
 double r_dfz2 = 3.1e-7;



////////////////////////////////////////////////////////////



double dfhplus_integral( double x, double s, double r){
    return x*(1-x)/(1-x+(r*x)-(s*x*(1-x)));
}

double dfhplus(double hRr, double hLi, double hRi, double hLr){

    return  -mN *((hRr*hLi)-(hRi*hLr))* trp(dfhplus_integral,s_dfhplus,r_dfhplus, 0.0,1.0,1000) /(16*pi*pi*(mHplus*mHplus));

}


/////////////////////////////////////////////////////////






double dfhpp_integral( double x, double s, double r){
    return (2*x*(1-x)/(1-x+(r*x)-(s*x*(1-x))))+(x*x/(1-x+(r*x)-(s*x*(1-x))));
}

double dfhpp(double hRr, double hLi, double hRi, double hLr){

    return -mtau *((hRr*hLi)-(hRi*hLr))* trp(dfhpp_integral,s_dfhpp,r_dfhpp, 0.0,1.0,1000) /(16*pi*pi*(mHpp*mHpp));

}






/////////////////////////////////////////////////////////






double dfhzero_integral( double x, double s, double r){
    return x*x/(1-x+(r*x)-(s*x*(1-x)));
}

double dfhzero(double hRr, double hLi, double hRi, double hLr){

    return  -mtau *((hRr*hLi)-(hRi*hLr))* trp(dfhzero_integral,s_dfhzero,r_dfhzero, 0.0,1.0,1000) /(16*pi*pi*(mHzero*mHzero));

}





/////////////////////////////////////////////////////////










double dfw2_integral( double x, double s, double r){
    return (1-x)*((3*(1-x))-(s*x*x))/(1-x+(r*x)-(s*x*(1-x)));
}

double dfw2 (double wRr, double wLi, double wRi, double wLr){

    return 0.5* mN *((wRr*wLi)-(wRi*wLr))* trp(dfw2_integral,s_dfw2,r_dfw2, 0.0,1.0,1000) /(16*pi*pi*(2000*2000));

}


/////////////////////////////////////////////////////////








double dfz2_integral( double x, double s, double r){
    return x*((4*(1-x))+(x*(r-s)))/(1-x+(r*x)-(s*x*(1-x)));
}

double dfz2 (double wRr, double wLi, double wRi, double wLr){

    return  mtau *((wRr*wLi)-(wRi*wLr))* trp(dfz2_integral,s_dfz2,r_dfz2, 0.0,1.0,1000) /(16*pi*pi*(mZ2*mZ2));

}





//////////////////////////////////////////////////////








///////////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
///////////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


//Function to be minimized. 

double input_funtion(double params[], const int param_dim, int *fcall, bool *quit, const bool validvector, void** context)
{ double result;
    double dftot;
	
  
 //Function should be written below.  
  
  dftot = dfhpp( params[0],  params[1],  params[2],  params[3])+ dfhplus( params[0],  params[1],  params[2],  params[3])+dfhzero( params[0],  params[1],  params[2],  params[3])+dfw2 ( params[4],  params[5],  params[6],  params[7])+dfz2 (params[4],  params[5],  params[6],  params[7]);
  //dftot = dfhpp( hRr,  hLi,  hRi,  hLr)+ dfhplus( hRr,  hLi,  hRi,  hLr)+dfhzero( hRr,  hLi,  hRi,  hLr)+dfw2 ( wRr,  wLi,  wRi,  wLr)+dfz2 ( wRr,  wLi,  wRi,  wLr);

  result = pow((dftot-(5e-15)),2);
  
 //End of the function. 

  if (!validvector) result = DBL_MAX;
  *fcall += 1;
  *quit = false;
  return result;
}


int main(int argc, char** argv)
{

int seed;

sscanf(argv[1],"%d",&seed);

  void* context = &input_funtion; //Not actually used in this example.
  cdiver(input_funtion, nPar, lowerbounds, upperbounds, path, nDerived, nDiscrete, discrete, partitionDiscrete,
         maxciv, maxgen, NP, nF, F, Cr, lambda, current, expon, bndry, jDE, lambdajDE, convthresh,
         convsteps, removeDuplicates, doBayesian, NULL, maxNodePop, Ztolerance, savecount, resume,
         outputSamples, init_pop_strategy, discard_unfit_points, max_init_attempts, max_acceptable_val, seed, context, verbose);
         //Note that prior, maxNodePop and Ztolerance are just ignored if doBayesian = false



}




