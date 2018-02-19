//
//  estimateT2.c
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/26/16.
//
//

#include "estimateT2.h"



// t2FitWelschNew_notCorr = -2*(TR-TE)./log(echo2L./(echo1L * (sind(alpha_deg_L/2)).^2 * (1 + exp(-TR/assumedT1))/(1 - cosd(alpha_deg_L)*exp(-TR/assumedT1))));


void computeT2map(float TR, float TE, float alphaL_deg, float assumedT1, float* SL1, float* SL2, int imageSize, int estT2type, float* estimatedT2) {
    
    float alphaL_rad = alphaL_deg * M_PI/180;
    int pp;     // Running index to loop through the pixels in the image.
    
//    int ppp = 118*256+50;
//    float testT2est = -2 * (TR-TE)/log(SL2[ppp]/(SL1[ppp]*pow(sin(alphaL_rad/2),2) * (1 + exp(-TR/assumedT1))/(1-cos(alphaL_rad)*exp(-TR/assumedT1))));
//    printf("SL1[50,118]: %g\n", SL1[ppp]);
//    printf("SL2[50,118]: %g\n", SL2[ppp]);
//    printf("T2est50,118]: %g\n", testT2est);
//    printf("sin(alphaL_rad/2): %g\n", sin(alphaL_rad/2));
//    printf("pow(sin(alphaL_rad/2),2): %g\n", pow(sin(alphaL_rad/2),2));
//    printf("exp(-TR/assumedT1): %g\n", exp(-TR/assumedT1));
//    printf("TR: %g\n", TR);
//    printf("TE: %g\n", TE);
//    printf("(1 + exp(-TR/assumedT1))/(1-cos(alphaL_rad)*exp(-TR/assumedT1)): %g\n", (1 + exp(-TR/assumedT1))/(1-cos(alphaL_rad)*exp(-TR/assumedT1)));
    
    
    for (pp=0; pp<imageSize; pp++){
//        if (SL1[pp] > noiseThreshold){
            estimatedT2[pp] = -2 * (TR-TE)/log(SL2[pp]/(SL1[pp]*pow(sin(alphaL_rad/2),2) * (1 + exp(-TR/assumedT1))/(1-cos(alphaL_rad)*exp(-TR/assumedT1))));
            estimatedT2[pp] *= 1000;   // Need to scale by 1000 because we are exporting to Dicoms.
//        } else {
//            estimatedT2[pp] = 0;
//        }

    }
    
}