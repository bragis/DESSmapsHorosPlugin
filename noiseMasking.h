//
//  maskNoise.h
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/28/16.
//
//

#ifndef noiseMasking_h
#define noiseMasking_h

#include <stdio.h>

void maskNoise(float* noisyImage, float* referenceImage, float maskingPercentage, int xsize, int ysize);

#endif /* noiseMasking_h */
