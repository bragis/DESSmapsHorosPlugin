//
//  estimateADC.c
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/26/16.
//
//

#include "estimateADC.h"


//void computeADCmap(float* SH1, float* SH2, float* SL1, float* SL2, int imageSize, int estADCtype, float noiseThreshold, float* estimatedADC) {
void generateDiffusionRatioArray(float TR_H, float TE_H, float alpha_H, float G_H, float Tg_H, float TR_L, float TE_L, float alpha_L, float G_L, float Tg_L, double T1start, double T1end, int NT1points, double T2start, double T2end, int NT2points, double Dstart, double Dend, int NDpoints, struct diffusionArrayPoint* diffusionRatioArray){
    printf("Have entered generateDiffusionRatioArray function.\n");
    
//    double T1start = 1.15;
//    double T1end = 1.25;
//    int NT1points = 2;
//    double T2start = 0.01;
//    double T2end = 0.08;
//    int NT2points = 70;
//    double Dstart = 0.01;
//    double Dend = 0.08;
//    int NDpoints = 70;
    printf("TR_H: %g\n", TR_H);
    printf("TE_H: %g\n", TE_H);
    printf("alpha_H: %g\n", alpha_H);
    printf("G_H: %g\n", G_H);
    printf("Tg_L: %g\n", Tg_L);
    printf("TR_L: %g\n", TR_L);
    printf("TE_L: %g\n", TE_L);
    printf("alpha_L: %g\n", alpha_L);
    printf("G_L: %g\n", G_L);
    printf("T1start: %g\n", T1start);
    printf("T1end: %g\n", T1end);
    printf("NT1points: %i\n", NT1points);
    printf("T2start: %g\n", T2start);
    printf("T2end: %g\n", T2end);
    printf("NT2points: %i\n", NT2points);
    printf("Dstart: %g\n", Dstart);
    printf("Dend: %g\n", Dend);
    printf("NDpoints: %i\n", NDpoints);
    
    int Nstates = 6;
    
    int ii,jj,kk;
    double T1_ijk, T2_ijk, D_ijk;
    double ReS1H,ImS1H,ReS2H,ImS2H,ReS1L,ImS1L,ReS2L,ImS2L;
    ReS1H = 0;
    ImS1H = 0;
    ReS2H = 0;
    ImS2H = 0;
    ReS1L = 0;
    ImS1L = 0;
    ReS2L = 0;
    ImS2L = 0;
    
//    double ReS1HSum,ImS1HSum,ReS2HSum,ImS2HSum,ReS1LSum,ImS1LSum,ReS2LSum,ImS2LSum;
//    ReS1HSum = 0;
//    ImS1HSum = 0;
//    ReS2HSum = 0;
//    ImS2HSum = 0;
//    ReS1LSum = 0;
//    ImS1LSum = 0;
//    ReS2LSum = 0;
//    ImS2LSum = 0;
    double diffusionRatioSum = 0;
    int index = 0;
    double T1step = ((double) (T1end-T1start))/( (double) (NT1points-1) );
    double T2step = ((double) (T2end-T2start))/( (double) (NT2points-1) );
    double Dstep = ((double) (Dend-Dstart))/( (double) (NDpoints-1) );
    //struct diffusionArrayPoint* diffusionRatioArray;
    
//    printf("Start allocating memory for the diffusion ratio array.\n");
//    diffusionRatioArray = malloc(NDpoints*sizeof(struct diffusionArrayPoint));
//    if (diffusionRatioArray == NULL) {
//        printf("Could not allocate memory for diffusion ratio array!!!\n");
//    }
//    printf("Done allocating memory for the diffusion ratio array.\n");
    
    
    //    for (i=0; i<NT1points; i++){
    //        T1_ijk = T1start + i*T1step;
    printf("Start constructing diffusion ratio array.\n");
    for (ii=0; ii<NDpoints; ii++){
        D_ijk = Dstart + ii*Dstep;
        printf("ii %i\n", ii);
        
        for (jj=0; jj<NT2points; jj++){
            T2_ijk = T2start + jj*T2step;
            
            for (kk=0; kk<NT1points; kk++){
                T1_ijk = T1start + kk*T1step;
                
                computeEchoesEPG(T1_ijk,T2_ijk,TR_H,TE_H,alpha_H,G_H,Tg_H,D_ijk,Nstates,&ReS1H,&ImS1H,&ReS2H,&ImS2H);
                computeEchoesEPG(T1_ijk,T2_ijk,TR_L,TE_L,alpha_L,G_L,Tg_L,D_ijk,Nstates,&ReS1L,&ImS1L,&ReS2L,&ImS2L);
                
                diffusionRatioSum += SC(ReS2H,ImS2H,ReS1H,ImS1H)*SC(ReS1L,ImS1L,ReS2L,ImS2L);
//                ReS1HSum += ReS1H;
//                ImS1HSum += ImS1H;
//                ReS2HSum += ReS2H;
//                ImS2HSum += ImS2H;
//                ReS1LSum += ReS1L;
//                ImS1LSum += ImS1L;
//                ReS2LSum += ReS2L;
//                ImS2LSum += ImS2L;
                
                //index = (i*NT2points*NDDpoints + j*NDDpoints + k);
                
                // Store values of T1, T2 and DD in the space array
                //T1T2DDspace[index].T1 = T1_ijk;
                //T1T2DDspace[index].T2 = T2_ijk;
                //T1T2DDspace[index].DD = DD_ijk;
                
                // Store echo ratios in the space
                // //T1T2DDspace[index].SC1 = SC(ReS2m,ImS2m,ReS1p,ImS1p);
                // //T1T2DDspace[index].SC2 = SC(ReS2p,ImS2p,ReS1p,ImS1p);
                // //T1T2DDspace[index].SC3 = SC(ReS2m,ImS2m,ReS1m,ImS1m);
                //T1T2DDspace[index].SC1 = SC(ReS1m,ImS1m,ReS1p,ImS1p);
                //T1T2DDspace[index].SC2 = SC(ReS2m,ImS2m,ReS2p,ImS2p);
                //T1T2DDspace[index].SC3 = SC(ReS1p,ImS1p,ReS2p,ImS2p);
                
                index++;
            }
        }
        diffusionRatioArray[ii].DD = D_ijk;
        //diffusionRatioArray[ii].RR = SC(ReS2H,ImS2H,ReS1H,ImS1H)*SC(ReS1L,ImS1L,ReS2L,ImS2L);
        // Compute the diffusion ratio. The idea is to average over T1 and T2 for each D value. However, we can just as well use sums instead of averages, it cancels out.
//        diffusionRatioArray[ii].RR = SC(ReS2HSum,ImS2HSum,ReS1HSum,ImS1HSum)*SC(ReS1LSum,ImS1LSum,ReS2LSum,ImS2LSum);
//        ReS1HSum = 0;
//        ImS1HSum = 0;
//        ReS2HSum = 0;
//        ImS2HSum = 0;
//        ReS1LSum = 0;
//        ImS1LSum = 0;
//        ReS2LSum = 0;
//        ImS2LSum = 0;
        diffusionRatioArray[ii].RR = diffusionRatioSum/((double) (NT1points*NT2points));
        diffusionRatioSum = 0;
        
        printf("diffusionRatioArray[ii]: %g, %g\n", diffusionRatioArray[ii].DD, diffusionRatioArray[ii].RR);
    }
    printf("Done constructing diffusion ratio array.\n");
}



//void computeADCmap(float TR_H, float TE_H, float alpha_H, float G_H, float Tg_H, float TR_L, float TE_L, float alpha_L, float G_L, float Tg_L, float* SH1, float* SH2, float* SL1, float* SL2, int imageSize, int estADCtype, float noiseThreshold, float* estimatedADC) {
void computeADCmap(float* SH1, float* SH2, float* SL1, float* SL2, int imageSize, int estADCtype, float* estimatedADC, struct diffusionArrayPoint* diffusionRatioArray, int NDpoints) {

    printf("Have entered computeADCmap function.\n");
    
//    double T1start = 1.15;
//    double T1end = 1.25;
//    double NT1points = 2;
//    double T2start = 0.01;
//    double T2end = 0.08;
//    double NT2points = 70;
//    double Dstart = 0.01;
//    double Dend = 0.08;
//    double NDpoints = 70;
//    
//    int Nstates = 6;
//    
//    int ii,jj,kk;
//    double T1_ijk, T2_ijk, D_ijk;
//    double ReS1H,ImS1H,ReS2H,ImS2H,ReS1L,ImS1L,ReS2L,ImS2L;
//    ReS1H = 0;
//    ImS1H = 0;
//    ReS2H = 0;
//    ImS2H = 0;
//    ReS1L = 0;
//    ImS1L = 0;
//    ReS2L = 0;
//    ImS2L = 0;
//    
//    double ReS1HSum,ImS1HSum,ReS2HSum,ImS2HSum,ReS1LSum,ImS1LSum,ReS2LSum,ImS2LSum;
//    ReS1HSum = 0;
//    ImS1HSum = 0;
//    ReS2HSum = 0;
//    ImS2HSum = 0;
//    ReS1LSum = 0;
//    ImS1LSum = 0;
//    ReS2LSum = 0;
//    ImS2LSum = 0;
//    int index = 0;
//    double T1step = ((double) (T1end-T1start))/( (double) (NT1points-1) );
//    double T2step = ((double) (T2end-T2start))/( (double) (NT2points-1) );
//    double Dstep = ((double) (Dend-Dstart))/( (double) (NDpoints-1) );
//    struct diffusionArrayPoint* diffusionRatioArray;
//    
//    printf("Start allocating memory for the diffusion ratio array.\n");
//    diffusionRatioArray = malloc(NDpoints*sizeof(struct diffusionArrayPoint));
//    if (diffusionRatioArray == NULL) {
//        printf("Could not allocate memory for diffusion ratio array!!!\n");
//    }
//    printf("Done allocating memory for the diffusion ratio array.\n");
//    
//    
////    for (i=0; i<NT1points; i++){
////        T1_ijk = T1start + i*T1step;
//    printf("Start constructing diffusion ratio array.\n");
//    for (ii=0; ii<NDpoints; ii++){
//        D_ijk = Dstart + ii*Dstep;
//        printf("ii %i\n", ii);
//    
//        for (jj=0; jj<NT2points; jj++){
//            T2_ijk = T2start + jj*T2step;
//            
//            for (kk=0; kk<NT1points; kk++){
//                T1_ijk = T1start + kk*T1step;
//                
//                computeEchoesEPG(T1_ijk,T2_ijk,TR_H,TE_H,alpha_H,G_H,Tg_H,D_ijk,Nstates,&ReS1H,&ImS1H,&ReS2H,&ImS2H);
//                computeEchoesEPG(T1_ijk,T2_ijk,TR_L,TE_L,alpha_L,G_L,Tg_L,D_ijk,Nstates,&ReS1L,&ImS1L,&ReS2L,&ImS2L);
//                
//                ReS1HSum += ReS1H;
//                ImS1HSum += ImS1H;
//                ReS2HSum += ReS2H;
//                ImS2HSum += ImS2H;
//                ReS1LSum += ReS1L;
//                ImS1LSum += ImS1L;
//                ReS2LSum += ReS2L;
//                ImS2LSum += ImS2L;
//
//                //index = (i*NT2points*NDDpoints + j*NDDpoints + k);
//                
//                // Store values of T1, T2 and DD in the space array
//                //T1T2DDspace[index].T1 = T1_ijk;
//                //T1T2DDspace[index].T2 = T2_ijk;
//                //T1T2DDspace[index].DD = DD_ijk;
//                
//                // Store echo ratios in the space
//                // //T1T2DDspace[index].SC1 = SC(ReS2m,ImS2m,ReS1p,ImS1p);
//                // //T1T2DDspace[index].SC2 = SC(ReS2p,ImS2p,ReS1p,ImS1p);
//                // //T1T2DDspace[index].SC3 = SC(ReS2m,ImS2m,ReS1m,ImS1m);
//                //T1T2DDspace[index].SC1 = SC(ReS1m,ImS1m,ReS1p,ImS1p);
//                //T1T2DDspace[index].SC2 = SC(ReS2m,ImS2m,ReS2p,ImS2p);
//                //T1T2DDspace[index].SC3 = SC(ReS1p,ImS1p,ReS2p,ImS2p);
//                
//                index++;
//            }
//        }
//        diffusionRatioArray[ii].DD = D_ijk;
//        diffusionRatioArray[ii].RR = SC(ReS2H,ImS2H,ReS1H,ImS1H)*SC(ReS1L,ImS1L,ReS2L,ImS2L);
//    }
//    printf("Done constructing diffusion ratio array.\n");
    
    // We've now generated an array of [D,diffusionRatio(D)]. Now we do the fit. This is done by looping through all pixels, and for every pixel, compute the estimated diffusion ratio, compare with the one in the array, and find the entry that gives the best match. This should give our diffusivity D.
    int ii, pp, ii_min;
    double E_ii, E_min, measuredDratio_pp;

    printf("Start fitting pixels.\n");
    //pp = 135*256+95;
    for (pp=0; pp<imageSize; pp++) {
        E_min = 1e6;
        ii_min = 0;
        //printf("pp %i\n", pp);
        // Compute the measured diffusion ratio
        measuredDratio_pp = (((double) SH2[pp])*((double) SL1[pp]))/(((double) SH1[pp])*((double) SL2[pp]));
        //printf("measuredDratio_pp[95,135] = %g\n", measuredDratio_pp);
    
        // Loop through the ratio array, compare to the measured results.
        for (ii=0; ii<NDpoints; ii++) {
            E_ii = pow(measuredDratio_pp - diffusionRatioArray[ii].RR,2);
            if (E_ii < E_min){
                E_min = E_ii;
                ii_min = ii;
            }
        }
        //printf("ii_min[95,135] = %g\n", measuredDratio_pp);
        
        estimatedADC[pp] = diffusionRatioArray[ii_min].DD*1e12;  // Need to scale by 1e9 so that I can display as Dicom.
    }
    printf("Done fitting pixels.\n");
}