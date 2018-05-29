/*   No-reference PSNR Estimation Algorithm for H.264 Encoded Video Sequences 
*    Copyright (C) 2018 Leonardo Favario <leonardo.favario@polito.it> 
*
*    If you use this program for research purposes, please cite the paper entitled 
*    M. Siekkinen, T. Kamarainen, L. Favario, E. Masala, "Can You See What I See? 
*    Quality-of-Experience Measurements of Mobile Live Video Broadcasting", ACM 
*    Transactions on Multimedia Computing, Communications, and Applications (TOMM), 
*    Volume 14 Issue 2s, May 2018. 
*    You may find it here: https://dx.doi.org/10.1145/3165279   
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This software is provided by the copyright holders and contributors “AS IS” and
*    any express or implied warranties, including, but not limited to, the implied
*    WARRANTIES of MERCHANTABILITY and FITNESS for a particular purpose are
*    DISCLAIMED. In NO EVENT shall the copyright owner or contributors be liable for
*    any DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, or CONSEQUENTIAL DAMAGES
*    (including, but not limited to, procurement of substitute goods or services;
*    loss of use, data, or profits; or business interruption) however caused and on
*    any theory of LIABILITY, whether in contract, strict liability, or tort
*    (including negligence or otherwise) arising in ANY WAY out of the use of this
*    software, even if advised of the possibility of such damage.
*
*    You should have received a copy of the GNU General Public License
*    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


// Includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Typedefs
typedef enum{false, true} boolean;

typedef struct anotherList{
    double mseBeta;
    double psnrBeta;
    double mseLambda;
    double psnrLambda;
    boolean predicted;
    struct anotherList* next;
} resultsList;

// Reference MSE struct
typedef struct {
    double mseBetaRef;
    double mseLambdaRef;
} mseRefRes;

// Global defines
#define DEBUG false 
// Number of quantized coefficients each line of the matrix
#define COEFFICIENTS 16
// Total number (square matrix)
#define TOTAL_COEFFICIENTS COEFFICIENTS * COEFFICIENTS
// Number of iterations in the Newton approximation
#define LIMIT_NEWTON_ITERATIONS 10000
// Starting value is 0.1. See (11)
#define BETA_START 0.1 
// Seeing the practical results, lambda starts with a lower value
#define LAMBDA_START  0.001
// Input file defines
#define MAX_LEN_FILE_NAME 1500
#define MAX_LEN_LINE 151

// Alpha parameter that controls the width of the dead zone around 0. See (3)
static const double ALPHA_INTER = (double) 5 / 6; 
static const double ALPHA_INTRA = (double)  2 / 3;

// Functions prototypes
double cauchy_distribution(double *a1, double *b1, int zeroCounter, int oneCounter, int qk, double alpha);
double laplace_distribution(double *b1, int zeroCounter, int oneCounter, int qk, double alpha);
double f(double *a1,double *b1, double B, int N0, int N1, int qk, double alpha);
double df(double *a1,double *b1, double B, int N0, int N1, int qk, double alpha);
double lpf(double *b1, double lambda_laplace, int N0, int N1, int qk,double alpha);
double lpdf(double lambda_laplace, int N0, int N1, int qk, double alpha );
// redefine pow(a,2) for better reading
double r(double number);
double check_macroblock_type(char* predmodestring);
double calculate_qk(int qp);
void freeResultsList(resultsList* head);
void insertResultsList(resultsList** head, double mseBeta, double mseLambda, double psnrBeta, double psnrLambda, boolean predicted);
int no_psnr_calculation(resultsList** results_list, int** coefficients, int qp, char* predmodestring, int typeval);
mseRefRes calc_average_mse(resultsList** head);
void mse_prediction(float skipRate, resultsList** results_list, mseRefRes mseRef);
boolean check_all_zeroes(int** coefficients);
