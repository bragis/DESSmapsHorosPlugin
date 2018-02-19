//
//  estimateT2.h
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/26/16.
//
//

#ifndef estimateT2_h
#define estimateT2_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void computeT2map(float TR, float TE, float alphaL, float assumedT1, float* SL1, float* SL2, int imageSize, int estT2type, float* estimatedT2);


#endif /* estimateT2_h */



