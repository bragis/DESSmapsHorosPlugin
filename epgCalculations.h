//
//  epgCalculations.h
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/26/16.
//
//

#ifndef epgCalculations_h
#define epgCalculations_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void computeEchoesEPG(double T1, double T2, double TR, double TE, double alpha_deg, double G, double Tg, double D, int Nstates, double *ReS1, double *ImS1, double *ReS2, double *ImS2);

double SC(double ReSA,double ImSA,double ReSB,double ImSB);

#endif /* epgCalculations_h */
