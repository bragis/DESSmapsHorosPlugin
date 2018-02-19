//
//  maskNoise.c
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/28/16.
//
//

#include "noiseMasking.h"

void maskNoise(float* noisyImage, float* referenceImage, float maskingPercentage, int xsize, int ysize) {

    // Find the maximum value of the reference image (this might for example be the S1L DESS image).
    int pp;
    float noiseThreshold;
    float maxValueReference = 0;
    for (pp = 0; pp < xsize*ysize; pp++){
        if (referenceImage[pp] > maxValueReference){
            maxValueReference = referenceImage[pp];
        }
    }
    
    // Set the noise threshold to a percentage of the maximum
    noiseThreshold = maskingPercentage*maxValueReference;
    
    // Zero out everything where the reference is below the threshold (these pixels are assumed to be noise).
    for (pp = 0; pp < xsize*ysize; pp++){
        if (referenceImage[pp] <= noiseThreshold){
            noisyImage[pp] = 0;
        }
    }
    
}