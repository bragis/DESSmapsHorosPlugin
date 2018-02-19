//
//  estimateADC.h
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/26/16.
//
//

#ifndef estimateADC_h
#define estimateADC_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "epgCalculations.h"

struct diffusionArrayPoint {
    double DD;
    double RR;
};

void generateDiffusionRatioArray(float TR_H, float TE_H, float alpha_H, float G_H, float Tg_H, float TR_L, float TE_L, float alpha_L, float G_L, float Tg_L, double T1start, double T1end, int NT1points, double T2start, double T2end, int NT2points, double Dstart, double Dend, int NDpoints, struct diffusionArrayPoint* diffusionRatioArray);

void computeADCmap(float* SH1, float* SH2, float* SL1, float* SL2, int imageSize, int estADCtype, float* estimatedADC, struct diffusionArrayPoint* diffusionRatioArray, int NDpoints);

#endif /* estimateADC_h */
