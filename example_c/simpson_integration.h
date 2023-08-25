#ifndef SIMPSON_H
#define SIMPSON_H


//trp(function, lower limit of variable, upper limit of variable, number of divisions(for accuracy purpose) )

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <float.h>
#include <stdbool.h>
#include <stddef.h>


double trp(double func( double x, double s, double r),double s, double r, double initial_x, double final_x, double divisions){
        
        double h;
        h = (final_x-initial_x)/divisions;
        double sum = 0.0;
        for(int i=0; i< (divisions/2);i++){
        sum = sum + (h/3)*(func((initial_x+(2*i*h)),s,r)+4*func((initial_x+((2*i+1)*h)),s,r)+func((initial_x+((2*i+2)*h)),s,r));
        }
        
        return sum;
}








#endif

