//
//  epgCalculations.c
//  DESSmaps
//
//  Created by Bragi Sveinsson on 4/26/16.
//
//

#include "epgCalculations.h"


////////////////////////////////////////////////////////////////////////////////////
// SC:		Compute the SC = |SA|/|SB|, where SA and SB are two echo measurements.
//		This quantity is used as a metric to determine how good our DESS fit is.
//
// Inputs:	The real and imaginary components of the signals SA and SB
//
// Returns:	|SA|/|SB|
////////////////////////////////////////////////////////////////////////////////////
double SC(double ReSA,double ImSA,double ReSB,double ImSB){
    return sqrt((pow(ReSA,2)+pow(ImSA,2))/(pow(ReSB,2)+pow(ImSB,2)));
}

int min(int a, int b){
    if (a>b) return b;
    else return a;
}

int max(int a, int b){
    if (a>b) return a;
    else return b;
}


//////////////////////////////////////////////////////////////////////////
// matrixMultiply:	Computes the matrix multiplication C=A*B. A,B and C
//			are stored as C++ arrays.
//
// inputs:		A:	The matrix A in the equation above
//			Arows:	Number of rows in A
//			Acols:  Number of columns in A
//			B:	The matrix B in the equation above
//			Brows:	Number of rows in B
//			Bcols:	Number of columns in B
//			C:	We store the result of the multiplication in this array.
//////////////////////////////////////////////////////////////////////////
void matrixMultiply(double* A, int Arows, int Acols, double* B, int Brows, int Bcols, double* C){
    
    double sum;
    int i,j,k;
    
    if (Acols == Brows){
        for (i=0; i<Arows; i++){
            for (j=0; j<Bcols; j++){
                sum = 0;
                for (k=0; k<Acols; k++){
                    sum += A[i*Acols+k]*B[k*Bcols+j];
                }
                C[i*Bcols+j] = sum;
            }
        }
    } else{
        printf("Matrix dimensions must match for matrix multiplication.\n");
    }
}


//////////////////////////////////////////////////////////////////////////
// matrixSum:		Computes the matrix sum C=A+B. A,B and C
//			are stored as C++ arrays.
//
// inputs:		A:	The matrix A in the equation above
//			Arows:	Number of rows in A
//			Acols:  Number of columns in A
//			B:	The matrix B in the equation above
//			Brows:	Number of rows in B
//			Bcols:	Number of columns in B
//			C:	We store the result of the summation in this array.
//////////////////////////////////////////////////////////////////////////
void matrixSum(double* A, int Arows, int Acols, double* B, int Brows, int Bcols, double* C){
    
    int i,j;
    
    if ((Acols == Bcols)&&(Arows==Brows)){
        for (i=0; i<Arows; i++){
            for (j=0; j<Acols; j++){
                C[i*Acols+j] = A[i*Acols+j] + B[i*Bcols+j];
            }
        }
    } else{
        printf("Matrix dimensions must match for matrix summation.\n");
    }
}


////////////////////////////////////////////////////////////
// Solve matrix equation A*x = b using LU decomposition.
// This is a slightly altered version of the LU decomposition
// method contained in jama_LU.h, developed by the National
// Institute of Standards and Technology (NIST)
////////////////////////////////////////////////////////////
void matrixSolve(double* A, double* x, double* b, int n, int m){
    
    int piv[n];
    double LUcolj[n];
    double LU_[n*m];
    int kmax, p;
    double s;
    
    
    
    int i;
    for (i=0; i<n*m; i++){
        LU_[i] = A[i];
    }
    
    int i11;
    for (i11=0; i11<n; i11++){
        piv[i11] = i11;
    }
    int pivsign = 1;
    
    int j;
    for (j=0; j<m;j++){
        
        
        int i10;
        for (i10=0; i10<n; i10++){
            LUcolj[i10] = LU_[i10*m+j];
        }
        
        // Apply previous transformations
        int i8;
        for (i8=0; i8<n; i8++){
            
            kmax = i8;
            if (j<i8) kmax = j;
            s = 0.0f;
            
            int k;
            for (k=0; k<kmax; k++){
                s += LU_[i8*m+k]*LUcolj[k];
            }
            
            LU_[i8*m+j] = LUcolj[i8] -= s;
        }
        
        // Find pivot and exchange if necessary
        
        p = j;
        
        int i6;
        for (i6=j+1; i6<n; i6++){
            if (fabs(LUcolj[i6]) > fabs(LUcolj[p])){
                p=i6;
            }
        }
        
        if (p!=j){
            int k=0;
            for (k=0; k<n; k++){
                double t = LU_[p*m+k];
                LU_[p*m+k] = LU_[j*m+k];
                LU_[j*m+k] = t;
            }
            k = piv[p];
            piv[p] = piv[j];
            piv[j] = k;
            pivsign = -pivsign;
        }
        
        // Compute multipliers
        if ((j<m) && (LU_[j*m+j] != 0.0)){
            int i;
            for (i=j+1; i<n; i++){
                LU_[i*m+j] /= LU_[j*m+j];
                //		printf("LU_[%i,%i] = %g\n", i,j,LU_[i*m+j]);
            }
        }
    }
    
    
    //x = permute_copy(b,piv);
    int i3;
    for (i3=0; i3<n; i3++){
        x[i3] = b[piv[i3]];
    }
    
    // Solve L*y = b(piv)
    int k;
    for (k=0; k<m; k++){
        int i;
        for (i=k+1; i<m; i++){
            x[i] -= x[k]*LU_[i*m+k];
        }
    }
    
    // Solve U*x = y;
    int k0;
    for (k0=m-1; k0>=0; k0--){
        x[k0] /= LU_[k0*m+k0];
        int i;
        for (i=0; i<k0; i++){
            x[i] -= x[k0]*LU_[i*m+k0];
        }
    }
}


////////////////////////////////////////////////////////////
// The function matrixSolve does not work for skinny matrices
// (only works for square matrices?). To solve A*x = b, where
// A is skinny, we find the least squares solution
//   x = A'b,
// where A' is the Moore-Penrose pseudo-inverse:
//   A' = (AT A)^(-1) AT
// (AT is A transposed)
// We do this by computing
//   ATA = AT*A and ATb = AT*b
// and then using matrixSolve on
//   ATAx = ATb.
////////////////////////////////////////////////////////////
void matrixLSsolve(double* A, double* x, double* b, int n, int m){
    double ATA[m*m];
    double ATb[m];
    double AT[m*n];
    int i, j;
    
    // Create AT, the transpose of A
    for (i=0; i<m; i++){
        for (j=0; j<n; j++){
            AT[i*n+j] = A[j*m+i];
            //printf("%g ", AT[i*n+j]);
        }
        //printf("\n");
    }
    //printf("\n");
    
    /*printf("A=\n");
     for (i=0; i<n; i++){
     for (j=0; j<m; j++){
     printf("%g ", A[i*m+j]);
     }
     printf("\n");
     }
     printf("\n");
     
     printf("b=\n");
     for (i=0; i<n; i++){
     printf("%g\n", b[i]);
     }*/
    
    // Create ATA = AT*A
    matrixMultiply(AT,m,n,A,n,m,ATA);
    
    /*for (i=0; i<m; i++){
     for (j=0; j<m; j++){
     printf("%g ", ATA[i*m+j]);
     }
     printf("\n");
     }
     printf("\n");*/
    
    // Create ATb = AT*b
    matrixMultiply(AT,m,n,b,n,1,ATb);
    
    /*for (i=0; i<m; i++){
     printf("%g ", ATb[i]);
     printf("\n");
     }
     printf("\n");*/
    
    // Solve ATA*x = ATb, i.e. find the least
    // squares solution to A*x = b
    matrixSolve(ATA,x,ATb,m,m);
    
    /*for (i=0; i<m; i++){
     printf("%g ", x[i]);
     printf("\n");
     }
     printf("\n");*/
}



/////////////////////////////////////////////////////////////////////////////////
// epg_arr_rf:	Computes the effects of the rf pulse on the state vector, i.e. computes
//		matrices A and B such that the effect of the pulse is described by
//		A*F+B. Assumes use of C++ arrays for matrices and vectors.
//
// inputs:	alpha:		flip angle (in degrees)
//		phi:		rf pulse orientation (normally zero)
//		Nstates:	The number of steady states assumed in the EPG computations
//		A:		Variable to store the matrix A in the operation A*F+B
//		B:		Variable to store the matrix B in the operation A*F+B
//////////////////////////////////////////////////////////////////////////////////
void epg_arr_rf(double alpha, double phi, int Nstates, double* A, double* B){
    
    double ReRR00 = pow((double) cos(alpha/2), (double) 2.0);
    //double ImRR00 = 0.0;
    double ReRR01 = cos(2*phi)*pow((double) sin(alpha/2),(double) 2.0);
    double ImRR01 = sin(2*phi)*pow((double) sin(alpha/2),(double) 2.0);
    double ReRR02 = sin(phi)*sin(alpha);
    double ImRR02 = -cos(phi)*sin(alpha);
    double ReRR10 = ReRR01;
    double ImRR10 = -ImRR01;
    double ReRR11 = ReRR00;
    //double ImRR11 = ImRR00;
    double ReRR12 = ReRR02;
    double ImRR12 = -ImRR02;
    double ReRR20 = -0.5*ReRR02;
    double ImRR20 = 0.5*ImRR02;
    double ReRR21 = -0.5*ReRR02;
    double ImRR21 = -0.5*ImRR02;
    double ReRR22 = cos(alpha);
    //double ImRR22 = 0.0;
    
    int matDim = 6*Nstates;
    int i;
    
    
    for (i = 0; i<2*Nstates; i+=2){
        A[i*matDim+i] = ReRR00;
        A[(i+1)*matDim + (i+1)] = ReRR00;
        
        A[(i+2*Nstates)*matDim+i] = ReRR10;
        A[(i+2*Nstates)*matDim+(i+1)] = -ImRR10;
        A[(i+1+2*Nstates)*matDim+(i+1)] = ReRR10;
        A[(i+1+2*Nstates)*matDim+i] = ImRR10;
        
        A[(i+4*Nstates)*matDim+i] = ReRR20;
        A[(i+4*Nstates)*matDim+(i+1)] = -ImRR20;
        A[(i+1+4*Nstates)*matDim+(i+1)] = ReRR20;
        A[(i+1+4*Nstates)*matDim+i] = ImRR20;
        
        A[i*matDim+(i+2*Nstates)] = ReRR01;
        A[i*matDim+(i+1+2*Nstates)] = -ImRR01;
        A[(i+1)*matDim+(i+1+2*Nstates)] = ReRR01;
        A[(i+1)*matDim+(i+2*Nstates)] = ImRR01;
        
        A[(i+2*Nstates)*matDim+(i+2*Nstates)] = ReRR11;
        A[(i+1+2*Nstates)*matDim+(i+1+2*Nstates)] = ReRR11;
        
        A[(i+4*Nstates)*matDim+(i+2*Nstates)] = ReRR21;
        A[(i+4*Nstates)*matDim+(i+1+2*Nstates)] = -ImRR21;
        A[(i+1+4*Nstates)*matDim+(i+1+2*Nstates)] = ReRR21;
        A[(i+1+4*Nstates)*matDim+(i+2*Nstates)] = ImRR21;
        
        A[i*matDim+(i+4*Nstates)] = ReRR02;
        A[i*matDim+(i+1+4*Nstates)] = -ImRR02;
        A[(i+1)*matDim+(i+1+4*Nstates)] = ReRR02;
        A[(i+1)*matDim+(i+4*Nstates)] = ImRR02;
        
        A[(i+2*Nstates)*matDim+(i+4*Nstates)] = ReRR12;
        A[(i+2*Nstates)*matDim+(i+1+4*Nstates)] = -ImRR12;
        A[(i+1+2*Nstates)*matDim+(i+1+4*Nstates)] = ReRR12;
        A[(i+1+2*Nstates)*matDim+(i+4*Nstates)] = ImRR12;
        
        A[(i+4*Nstates)*matDim+(i+4*Nstates)] = ReRR22;
        A[(i+1+4*Nstates)*matDim+(i+1+4*Nstates)] = ReRR22;
    }
    
    memset(B,0,matDim*sizeof(double));
}


/////////////////////////////////////////////////////////////////////////////////
// epg_arr_grelax: Computes the effects of relaxation and diffusion with or without
//		an applied gradient on the state vector, i.e. computes
//		matrices A and B such that the effect of the relaxation and diffusion
//		is described by A*F+B. Assumes use of C++ arrays for matrix computations.
//
// inputs:	T1:		The T1 value at the point (in sec)
//		T2:		The T2 value at the point (in sec)
//		t:		The length of the time interval being computed
//		dk:		The phase change per sec if a gradient is applied
//		D:		The diffusion constant at the point (m^2/sec)
//		G_on:		1 if gradient on, 0 otherwise
//		Nstates:	The number of steady states assumed in EPG computations
//		A:		Variable for storing the matrix A in the operation A*F+B
//		B:		Variable for storing the matrix B in the operation A*F+B
//
//////////////////////////////////////////////////////////////////////////////////
void epg_arr_grelax(double T1,double T2,double t,double dk,double D,int G_on,int Nstates, double* A, double* B){
    
    // Create the matrices Ar,Br,Ad,Bd, where the relaxation effects
    // are described by Ar*F+Br and the diffusion effects are described by
    // Ad*F+Bd
    int matDim = 6*Nstates;
    double Ar [matDim*matDim];
    double Br [matDim];
    double Ad [matDim*matDim];
    double Bd [matDim];
    
    // Set the matrices to zero
    memset(Ar,0,matDim*matDim*sizeof(double));
    memset(Br,0,matDim*sizeof(double));
    memset(Ad,0,matDim*matDim*sizeof(double));
    memset(Bd,0,matDim*sizeof(double));
    
    // Pre-compute the relaxation coefficients to save time
    double E2 = exp(-t/T2);
    double E1 = exp(-t/T1);
    
    // Construct the relaxation matrix Ar
    int i,j;
    for (i=0;i<2*Nstates;i++){
        Ar[i*matDim+i] = E2;
        Ar[(i+2*Nstates)*matDim+(i+2*Nstates)] = E2;
        Ar[(i+4*Nstates)*matDim+(i+4*Nstates)] = E1;
    }
    
    // Construct relaxation matrix Br
    Br[4*Nstates] = 1-E1;
    
    // Construct diffusion matrix Ad
    double eBvFpD, eBvFmD, eBvZD;
    for (i=0; i<2*Nstates; i+=2){
        eBvFpD = exp(-D*(pow((i/2+0.5*G_on)*dk,2)+G_on*pow(dk,2)/12.0)*t);
        eBvFmD = exp(-D*(pow((-i/2+0.5*G_on)*dk,2)+G_on*pow(dk,2)/12.0)*t);
        eBvZD = exp(-D*pow((i/2)*dk,2)*t);
        Ad[i*matDim+i] = eBvFpD;
        Ad[(i+1)*matDim+(i+1)] = eBvFpD;
        Ad[(i+2*Nstates)*matDim+(i+2*Nstates)] = eBvFmD;
        Ad[(i+2*Nstates+1)*matDim+(i+2*Nstates+1)] = eBvFmD;
        Ad[(i+4*Nstates)*matDim+(i+4*Nstates)] = eBvZD;
        Ad[(i+4*Nstates+1)*matDim+(i+4*Nstates+1)] = eBvZD;
    }
    
    // Create diffusion matrix Bd
    // -- Bd is zero, nothing needs to be done --
    
    // Compute A and B from
    // 	A = Ad*Ar
    //	B = Ad*Br + Bd
    matrixMultiply(Ad,matDim,matDim,Ar,matDim,matDim,A);
    matrixMultiply(Ad,matDim,matDim,Br,matDim,1,B);
    
    
    // Create state shifting matrix Ag if gradient is on
    if (G_on){
        double Ag [matDim*matDim];
        
        memset(Ag,0,matDim*matDim*sizeof(double));
        
        for (i=0; i<2*Nstates-2; i++){
            Ag[(i+2)*matDim+i] = 1;
            Ag[(i+2*Nstates)*matDim+(i+2+2*Nstates)] = 1;
        }
        
        for (i=4*Nstates; i<6*Nstates; i++){
            Ag[i*matDim+i] = 1;
        }
        
        Ag[2*Nstates+2] = 1;
        Ag[matDim+2*Nstates+3] = -1;
        
        // The total matrices are then A = Ag*Ad*Ar and B = Ag*(Ad*Br+Bd)
        double tempA [matDim*matDim];
        for (i=0;i<matDim;i++){
            for (j=0; j<matDim; j++){
                tempA[i*matDim+j] = A[i*matDim+j];
            }
        }
        
        double tempB [matDim];
        for (i=0;i<matDim;i++){
            tempB[i] = B[i];
        }
        
        matrixMultiply(Ag,matDim,matDim,tempA,matDim,matDim,A);
        matrixMultiply(Ag,matDim,matDim,tempB,matDim,1,B);
    }
}



//////////////////////////////////////////////////////////////////////////////////////////
// computeEchoesEPG:	Given a description of a DESS sequence applied to a point with given
//			values of T1, T2 and diffusion constant D, estimate the two echo signals
//			obtained by that sequence. This function is very similar to epg_dess_mat, but
//			does not use the TNT and JAMA libraries except for the matrix inversion, which
//			is done using an LU decomposition instead of an SVD decomposition.
//
// Inputs:		T1 value at given point,    (S)
//			T2 value at given point,        (S)
//			TR (repetition time of the sequence),   (S)
//			TE (the echo time, the time from the first RF pulse to the first echo) (S)
//			alpha_deg (the flip angle in degrees)   (deg)
//			G (the gradient strength)               (G/M)
//			Tg (the time of the applied gradient)   (S)
//			D (the diffusion constant at the given point),  (m^2/s)
//			Nstates (the echo values are computed using an EPG procedure with Nstates number of states)
//			ReS1,ImS1,ReS2,ImS2 (variables passed by reference that store the real and imaginary part of
//			the first and second echoes)        (set at 6)
////////////////////////////////////////////////////////////////////////////////////////////
void computeEchoesEPG(double T1, double T2, double TR, double TE, double alpha_deg, double G, double Tg, double D, int Nstates, double *ReS1, double *ImS1, double *ReS2, double *ImS2){
    
    double PI = 3.14159265358979;
    double alpha = PI/180.0*alpha_deg;		// Convert alpha to radians
    double gamma = 2*PI*4258;			// Gamma, rad/(G*s)
    double dk = gamma*G*Tg;				// Step size in k-space
    int matDim = 6*Nstates;
    
    // We create a number of matrices for the EPG computation.
    //	The effect of the RF pulse on the state F is represented by A1*F + B1
    //	The period from the RF pulse to the first echo gives A2*F + B2
    //	The period from the first echo to the gradient gives A3*F + B3
    //	The gradient period gives A4*F + B4
    //	The period from the end of the gradient to the second echo gives A5*F + B5
    //	The period from the second echo to the next RF pulse gives A6*F + B6
    double A1 [matDim*matDim];
    double B1 [matDim];
    double A2 [matDim*matDim];
    double B2 [matDim];
    double A3 [matDim*matDim];
    double B3 [matDim];
    double A4 [matDim*matDim];
    double B4 [matDim];
    double A5 [matDim*matDim];
    double B5 [matDim];
    double A6 [matDim*matDim];
    double B6 [matDim];
    double AA [matDim*matDim];
    double BB [matDim];
    double AA1 [matDim*matDim];
    double BB1 [matDim];
    double F1 [matDim];
    double F2 [matDim];
    
    // Set all the matrices to zero
    memset(A1,0,matDim*matDim*sizeof(double));
    memset(A2,0,matDim*matDim*sizeof(double));
    memset(A3,0,matDim*matDim*sizeof(double));
    memset(A4,0,matDim*matDim*sizeof(double));
    memset(A5,0,matDim*matDim*sizeof(double));
    memset(A6,0,matDim*matDim*sizeof(double));
    memset(AA,0,matDim*matDim*sizeof(double));
    memset(AA1,0,matDim*matDim*sizeof(double));
    memset(B1,0,matDim*sizeof(double));
    memset(B2,0,matDim*sizeof(double));
    memset(B3,0,matDim*sizeof(double));
    memset(B4,0,matDim*sizeof(double));
    memset(B5,0,matDim*sizeof(double));
    memset(B6,0,matDim*sizeof(double));
    memset(BB,0,matDim*sizeof(double));
    memset(BB1,0,matDim*sizeof(double));
    
    // Compute the values of the matrices using the functions epg_arr_rf and epg_arr_grelax
    epg_arr_rf(alpha,0,Nstates,A1,B1);
    epg_arr_grelax(T1,T2,TE,dk,D,0,Nstates,A2,B2);
    epg_arr_grelax(T1,T2,(TR-Tg)/2.0f-TE,dk,D,0,Nstates,A3,B3);
    epg_arr_grelax(T1,T2,Tg,dk,D,1,Nstates,A4,B4);
    epg_arr_grelax(T1,T2,(TR-Tg)/2.0-TE,dk,D,0,Nstates,A5,B5);
    epg_arr_grelax(T1,T2,TE,dk,D,0,Nstates,A6,B6);
    
    
    //printf("A1:\n");
    //int ind1, ind2;
    //for (ind1 = 0; ind1<matDim; ind1++){
    //	for (ind2 = 0; ind2<matDim; ind2++){
    //		printf("%g,", A1[ind1*matDim+ind2]);
    //	}
    //	printf("\n");
    //}
    
    
    // Compute the resulting matrices for an entire gradient period, AA and BB,
    // and the matrices from the first echo to the second echo, AA1 and BB1
    //
    // Note: Ideally this should be done more concisely
    //
    // AA = A2*A1*A6*A5*A4*A3
    // BB = A2*(A1*(A6*(A5*(A4*B3+B4)+B5)+B6)+B1)+B2
    //
    double A4A3 [matDim*matDim];
    double A6A5A4A3 [matDim*matDim];
    double A1A6A5A4A3 [matDim*matDim];
    
    matrixMultiply(A4,matDim,matDim,A3,matDim,matDim,A4A3);
    matrixMultiply(A5,matDim,matDim,A4A3,matDim,matDim,AA1);
    matrixMultiply(A6,matDim,matDim,AA1,matDim,matDim,A6A5A4A3);
    matrixMultiply(A1,matDim,matDim,A6A5A4A3,matDim,matDim,A1A6A5A4A3);
    matrixMultiply(A2,matDim,matDim,A1A6A5A4A3,matDim,matDim,AA);
    
    double A4B3 [matDim];
    double A4B3pB4 [matDim];
    double A5A4B3pB4 [matDim];
    double A6A5A4B3pB4pB5 [matDim];
    double A6A5A4B3pB4pB5pB6 [matDim];
    double A1A6A5A4B3pB4pB5pB6 [matDim];
    double A1A6A5A4B3pB4pB5pB6pB1 [matDim];
    double A2A1A6A5A4B3pB4pB5pB6pB1 [matDim];
    
    matrixMultiply(A4,matDim,matDim,B3,matDim,1,A4B3);
    matrixSum(A4B3,matDim,1,B4,matDim,1,A4B3pB4);
    matrixMultiply(A5,matDim,matDim,A4B3pB4,matDim,1,A5A4B3pB4);
    matrixSum(A5A4B3pB4,matDim,1,B5,matDim,1,BB1);
    matrixMultiply(A6,matDim,matDim,BB1,matDim,1,A6A5A4B3pB4pB5);
    matrixSum(A6A5A4B3pB4pB5,matDim,1,B6,matDim,1,A6A5A4B3pB4pB5pB6);
    matrixMultiply(A1,matDim,matDim,A6A5A4B3pB4pB5pB6,matDim,1,A1A6A5A4B3pB4pB5pB6);
    matrixSum(A1A6A5A4B3pB4pB5pB6,matDim,1,B1,matDim,1,A1A6A5A4B3pB4pB5pB6pB1);
    matrixMultiply(A2,matDim,matDim,A1A6A5A4B3pB4pB5pB6pB1,matDim,1,A2A1A6A5A4B3pB4pB5pB6pB1);
    matrixSum(A2A1A6A5A4B3pB4pB5pB6pB1,matDim,1,B2,matDim,1,BB);
    
    // Construct (I-AA)
    double ImAA[matDim*matDim];
    int i;
    for (i=0; i<matDim*matDim; i++){
        ImAA[i] = -AA[i];
    }
    
    for (i=0; i<matDim; i++){
        ImAA[i*matDim+i] += 1.0f;
    }
    
    // Solve F1 = AA*F1 + BB (or (I-AA)*F1 = BB
    matrixSolve(ImAA, F1, BB, matDim, matDim);
    
    
    // Compute F2 = AA1*F1 + BB1
    double AA1F1 [matDim];
    matrixMultiply(AA1,matDim,matDim,F1,matDim,1,AA1F1);
    matrixSum(AA1F1,matDim,1,BB1,matDim,1,F2);
    
    // Get the corresponding echo values from the state vectors
    *ReS1 = F1[0];
    *ImS1 = F1[1];
    *ReS2 = F2[0];
    *ImS2 = F2[1];
    
}

//////////////////////////////////////////////////////////////////////////////////////////
// computeApproxEchosEPG:	Given a description of a DESS sequence applied to a point with given
//			values of T1, T2 and diffusion constant D, estimate the two echo signals
//			obtained by that sequence. This method does not use inverse matrixes, so it runs much faster in theory
//
// Inputs:		T1 value at given point,    (S)
//			T2 value at given point,        (S)
//			TR (repetition time of the sequence),   (S)
//			TE (the echo time, the time from the first RF pulse to the first echo) (S)
//			alpha_deg (the flip angle in degrees)   (deg)
//			G (the gradient strength)               (G/M)
//			Tg (the time of the applied gradient)   (S)
//			D (the diffusion constant at the given point),  (m^2/s)
//			Nstates (the echo values are computed using an EPG procedure with Nstates number of states)
//			ReS1,ImS1,ReS2,ImS2 (variables passed by reference that store the real and imaginary part of
//			the first and second echoes)        (set at 6)
////////////////////////////////////////////////////////////////////////////////////////////
void computeApproxEchosEPG(double T1, double T2, double TR, double TE, double alpha_deg, double G, double Tg, double D, int Nstates, double *ReS1, double *ImS1, double *ReS2, double *ImS2){
    
    double PI = 3.14159265358979;
    double alpha = PI/180.0*alpha_deg;		// Convert alpha to radians
    double gamma = 2*PI*4258;			// Gamma, rad/(G*s)
    double dk = gamma*G*Tg;				// Step size in k-space
    int matDim = 6*Nstates;
    
    // We create a number of matrices for the EPG computation.
    //	The effect of the RF pulse on the state F is represented by A1*F + B1
    //	The period from the RF pulse to the first echo gives A2*F + B2
    //	The period from the first echo to the gradient gives A3*F + B3
    //	The gradient period gives A4*F + B4
    //	The period from the end of the gradient to the second echo gives A5*F + B5
    //	The period from the second echo to the next RF pulse gives A6*F + B6
    double A1 [matDim*matDim];
    double B1 [matDim];
    double A2 [matDim*matDim];
    double B2 [matDim];
    double A3 [matDim*matDim];
    double B3 [matDim];
    double A4 [matDim*matDim];
    double B4 [matDim];
    double A5 [matDim*matDim];
    double B5 [matDim];
    double A6 [matDim*matDim];
    double B6 [matDim];
    double AA [matDim*matDim];
    double BB [matDim];
    double AA1 [matDim*matDim];
    double BB1 [matDim];
    double F1 [matDim];
    double F2 [matDim];
    
    // Set all the matrices to zero
    memset(A1,0,matDim*matDim*sizeof(double));
    memset(A2,0,matDim*matDim*sizeof(double));
    memset(A3,0,matDim*matDim*sizeof(double));
    memset(A4,0,matDim*matDim*sizeof(double));
    memset(A5,0,matDim*matDim*sizeof(double));
    memset(A6,0,matDim*matDim*sizeof(double));
    memset(AA,0,matDim*matDim*sizeof(double));
    memset(AA1,0,matDim*matDim*sizeof(double));
    memset(B1,0,matDim*sizeof(double));
    memset(B2,0,matDim*sizeof(double));
    memset(B3,0,matDim*sizeof(double));
    memset(B4,0,matDim*sizeof(double));
    memset(B5,0,matDim*sizeof(double));
    memset(B6,0,matDim*sizeof(double));
    memset(BB,0,matDim*sizeof(double));
    memset(BB1,0,matDim*sizeof(double));
    
    // Compute the values of the matrices using the functions epg_arr_rf and epg_arr_grelax
    epg_arr_rf(alpha,0,Nstates,A1,B1);
    epg_arr_grelax(T1,T2,TE,dk,D,0,Nstates,A2,B2);
    epg_arr_grelax(T1,T2,(TR-Tg)/2.0f-TE,dk,D,0,Nstates,A3,B3);
    epg_arr_grelax(T1,T2,Tg,dk,D,1,Nstates,A4,B4);
    epg_arr_grelax(T1,T2,(TR-Tg)/2.0-TE,dk,D,0,Nstates,A5,B5);
    epg_arr_grelax(T1,T2,TE,dk,D,0,Nstates,A6,B6);
    
    
    //printf("A1:\n");
    //int ind1, ind2;
    //for (ind1 = 0; ind1<matDim; ind1++){
    //	for (ind2 = 0; ind2<matDim; ind2++){
    //		printf("%g,", A1[ind1*matDim+ind2]);
    //	}
    //	printf("\n");
    //}
    
    
    // Compute the resulting matrices for an entire gradient period, AA and BB,
    // and the matrices from the first echo to the second echo, AA1 and BB1
    //
    // Note: Ideally this should be done more concisely
    //
    // AA = A2*A1*A6*A5*A4*A3
    // BB = A2*(A1*(A6*(A5*(A4*B3+B4)+B5)+B6)+B1)+B2
    //
    double A4A3 [matDim*matDim];
    double A6A5A4A3 [matDim*matDim];
    double A1A6A5A4A3 [matDim*matDim];
    
    matrixMultiply(A4,matDim,matDim,A3,matDim,matDim,A4A3);
    matrixMultiply(A5,matDim,matDim,A4A3,matDim,matDim,AA1);
    matrixMultiply(A6,matDim,matDim,AA1,matDim,matDim,A6A5A4A3);
    matrixMultiply(A1,matDim,matDim,A6A5A4A3,matDim,matDim,A1A6A5A4A3);
    matrixMultiply(A2,matDim,matDim,A1A6A5A4A3,matDim,matDim,AA);
    
    double A4B3 [matDim];
    double A4B3pB4 [matDim];
    double A5A4B3pB4 [matDim];
    double A6A5A4B3pB4pB5 [matDim];
    double A6A5A4B3pB4pB5pB6 [matDim];
    double A1A6A5A4B3pB4pB5pB6 [matDim];
    double A1A6A5A4B3pB4pB5pB6pB1 [matDim];
    double A2A1A6A5A4B3pB4pB5pB6pB1 [matDim];
    
    matrixMultiply(A4,matDim,matDim,B3,matDim,1,A4B3);
    matrixSum(A4B3,matDim,1,B4,matDim,1,A4B3pB4);
    matrixMultiply(A5,matDim,matDim,A4B3pB4,matDim,1,A5A4B3pB4);
    matrixSum(A5A4B3pB4,matDim,1,B5,matDim,1,BB1);
    matrixMultiply(A6,matDim,matDim,BB1,matDim,1,A6A5A4B3pB4pB5);
    matrixSum(A6A5A4B3pB4pB5,matDim,1,B6,matDim,1,A6A5A4B3pB4pB5pB6);
    matrixMultiply(A1,matDim,matDim,A6A5A4B3pB4pB5pB6,matDim,1,A1A6A5A4B3pB4pB5pB6);
    matrixSum(A1A6A5A4B3pB4pB5pB6,matDim,1,B1,matDim,1,A1A6A5A4B3pB4pB5pB6pB1);
    matrixMultiply(A2,matDim,matDim,A1A6A5A4B3pB4pB5pB6pB1,matDim,1,A2A1A6A5A4B3pB4pB5pB6pB1);
    matrixSum(A2A1A6A5A4B3pB4pB5pB6pB1,matDim,1,B2,matDim,1,BB);
    
    // Construct (I-AA)
    double ImAA[matDim*matDim];
    int i;
    for (i=0; i<matDim*matDim; i++){
        ImAA[i] = -AA[i];
    }
    
    for (i=0; i<matDim; i++){
        ImAA[i*matDim+i] += 1.0f;
    }
    
    // Solve F1 = AA*F1 + BB (or (I-AA)*F1 = BB
    //matrixSolve(ImAA, F1, BB, matDim, matDim);
    
    //compute F1
    // F1 = I+A+A^2+A^3+...+A^nstates
    
    //create New variable for loop
    int pwr;
    //create summ of the A's
    double Apwrs [matDim*matDim];
    double C [matDim*matDim];
    //create increasing powers of A's
    double Summation [matDim*matDim];
    double S [matDim*matDim];
    for (pwr=0; pwr<Nstates; pwr++){
        matrixMultiply(Apwrs, matDim, matDim, AA, matDim, matDim, C);
        matrixSum(C, matDim, matDim, Summation, matDim, matDim, S);
        Apwrs [matDim*matDim] = C [matDim*matDim];
        Summation [matDim*matDim] = S [matDim*matDim];
    }
    
    F1 [matDim] = Summation [matDim*matDim];
    
    // Compute F2 = AA1*F1 + BB1
    double AA1F1 [matDim];
    matrixMultiply(AA1,matDim,matDim,F1,matDim,1,AA1F1);
    matrixSum(AA1F1,matDim,1,BB1,matDim,1,F2);
    
    // Get the corresponding echo values from the state vectors
    *ReS1 = F1[0];
    *ImS1 = F1[1];
    *ReS2 = F2[0];
    *ImS2 = F2[1];
    
}