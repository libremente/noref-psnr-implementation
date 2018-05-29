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

/********************************************
 *  General Info: Parser including the NO_REF calculations 
 *  This program parses an XML file, 
 *  extracts the MacroBlocks' quantized coefficients to
 *  perform the MSE and PSNR using a no-ref approach calculation. 
 *  This implementation is based on the Brandao/Queluz paper called: 
 *  "No-ref psnr estimation algorithm". 
 *  Run the makefile or:
 *   $ gcc no_ref_psnr.c no_ref_impl.c -Wall -lm -o no_ref_psnr 
 *   $ ./no_ref_psnr input_file_name.xml output_file_name 
 ********************************************/

// Includes
#include "no_ref_header.h"

// Main
int main(int argc, char **argv){

    // Declare variables
    char fname[MAX_LEN_FILE_NAME], fnameout[MAX_LEN_FILE_NAME], tmp[MAX_LEN_FILE_NAME];
    FILE* fr = NULL;
    FILE* foutMseLambda = NULL;
    FILE* foutMseBeta = NULL;
    FILE* foutPsnrBeta = NULL;
    FILE* foutPsnrLambda = NULL;
    FILE* foutFramePsnrLambda = NULL;
    FILE* foutFramePsnrBeta = NULL;
    FILE* foutFrameMseLambda = NULL;
    FILE* foutFrameMseBeta = NULL;

    // Check arguments provided
    if(argc != 3){
        fprintf(stdout, "Parses the Ghent XML H.264 format to extract QP etc. for each MB - and output some statistics\n");
        fprintf(stdout, "REQUIRES each tag to be on a separate line!!!\n");
        fprintf(stdout, "ERROR: It should be ./executable input_file.xml output_file_name\n");
        exit(1);
    }

    // Open files
    strncpy(fname, argv[1], MAX_LEN_FILE_NAME);
    fr = fopen(fname, "r");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open input file. Check if %s exists in the current path.\n", fname);
        exit(1);
    }

    strncpy(fnameout, argv[2], MAX_LEN_FILE_NAME);
    strcpy(tmp, fnameout);
    strcat(tmp,".mseLambda");
    foutMseLambda = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }

    memset(tmp,0,strlen(tmp));
    strcpy(tmp, fnameout);
    strcat(tmp,".mseBeta");
    foutMseBeta = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }

    memset(tmp,0,strlen(tmp));
    strcpy(tmp, fnameout);
    strcat(tmp,".psnrBeta");
    foutPsnrBeta = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }

    memset(tmp,0,strlen(tmp));
    strcpy(tmp, fnameout);
    strcat(tmp,".framePsnrLambda");
    foutFramePsnrLambda = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }

    memset(tmp,0,strlen(tmp));
    strcpy(tmp, fnameout);
    strcat(tmp,".frameMseLambda");
    foutFrameMseLambda = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }

    memset(tmp,0,strlen(tmp));
    strcpy(tmp, fnameout);
    strcat(tmp,".psnrLambda");
    foutPsnrLambda = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }

    memset(tmp,0,strlen(tmp));
    strcpy(tmp, fnameout);
    strcat(tmp,".framePsnrBeta");
    foutFramePsnrBeta = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }

    memset(tmp,0,strlen(tmp));
    strcpy(tmp, fnameout);
    strcat(tmp,".frameMseBeta");
    foutFrameMseBeta = fopen(tmp, "w");
    if(fr == NULL){
        fprintf(stdout, "Impossible to open output file. Permission error or out of memory problem.\n");
        exit(1);
    }


    if(strstr(fname, "-h") != NULL){
        fprintf(stdout, "N/A\n");
        exit(0);
    }

    // Prints
    if(DEBUG){
        fprintf(foutMseLambda, "MSE calculated using Lambda, one for each MB");
        fprintf(foutMseBeta, "MSE calculated using Beta, one for each MB");
        fprintf(foutPsnrLambda, "PSNR calculated using Lambda, one for each MB");
        fprintf(foutPsnrBeta, "PSNR calculated using Beta, one for each MB");
    }

    // Declare other variables
    boolean inPicture = false, inSlice = false, inMB = false, inRow = false;

    char* line = NULL, *allocLine = NULL, *pos_pic_id = NULL, *pos_pic_id_end = NULL,
          *predmodestring = NULL;
    size_t len = 0;
    ssize_t read = 0;
    int numValues = 0, qp = 0, i = 0, rowCounter = 0, typeval  = 0, skipflag = 0,
        mbcnt = 0, skipCnt = 0, pic_id = 0, picCounter = 0;
    // Reference mse struct 
    mseRefRes mseRef;
    mseRef.mseBetaRef = 0;
    mseRef.mseLambdaRef = 0;

    // Matrix pointer
    int **coefficients = NULL;

    // Array for storing the results
    resultsList* results_list = NULL; 

    // Start parsing line
    read = getline(&line, &len, fr);

    while (read != -1) {
        if (inPicture == false) {
            if (strstr(line, "<Picture") != NULL) {
                // find PIC_ID
                pos_pic_id = strstr(line, "id=");
                pos_pic_id_end = strstr(pos_pic_id+4, "\"");
                numValues = pos_pic_id_end - pos_pic_id - 4;
                allocLine = (char*) malloc((numValues+1) * sizeof(char));
                i = 0;
                while (i < numValues) {
                    allocLine[i] = pos_pic_id[4+i];
                    i++;
                }
                allocLine[i] = '\0';
                pic_id = atoi(allocLine);
                free(allocLine);
                allocLine = NULL;

                // Info print
                fprintf(stdout, "Processing MBs of Pic: %d\n", pic_id);

                // Set flag
                inPicture = true;
            }
        }
        else {
            if (strstr(line, "</Picture") != NULL) {
                inPicture = false;

                // Save MSE of the first I frame which will be taken as the
                // reference. 
                // TODO: is this the only reference?
                if (picCounter == 0) {
                    mseRef = calc_average_mse(&results_list);
                }

                // Increment picture counter
                picCounter++;

                float totFrameMseBeta = 0.0;
                float totFrameMseLambda = 0.0;
                int totCnt = 0;

                // Print if the array contains something
                if (results_list != NULL) {

                    resultsList* cursor  = results_list;
                    while (cursor != NULL) {
                        totFrameMseBeta += cursor->mseBeta;
                        totFrameMseLambda += cursor->mseLambda;
                        totCnt++;

                        fprintf(foutMseLambda, "%lf ", cursor->mseLambda);
                        fprintf(foutMseBeta, "%lf ", cursor->mseBeta);
                        fprintf(foutPsnrBeta, "%lf ", cursor->psnrBeta);
                        fprintf(foutPsnrLambda, "%lf ", cursor->psnrLambda);

                        if (cursor->next == NULL) {
                            totFrameMseBeta += cursor->mseBeta;
                            totFrameMseLambda += cursor->mseLambda;
                            totCnt++;

                            fprintf(foutMseLambda, "%lf\n", cursor->mseLambda);
                            fprintf(foutMseBeta, "%lf\n", cursor->mseBeta);
                            fprintf(foutPsnrBeta, "%lf\n", cursor->psnrBeta);
                            fprintf(foutPsnrLambda, "%lf\n", cursor->psnrLambda);
                            break;
                        }
                        cursor = cursor->next;
                    }

                    totFrameMseBeta /= totCnt;
                    totFrameMseLambda /= totCnt;
                    fprintf(foutFrameMseBeta, "%lf\n", totFrameMseBeta);
                    fprintf(foutFrameMseLambda, "%lf\n", totFrameMseLambda);
                    float psnrFrameBeta = 1000.0;
                    if( mseRef.mseBetaRef > 0.0) {
                        psnrFrameBeta = 10.0 * log10 (255.0*255.0/ totFrameMseBeta);
                    }
                    float psnrFrameLambda = 1000.0;
                    if(mseRef.mseLambdaRef > 0.0) {
                        psnrFrameLambda = 10.0 * log10 (255.0*255.0/ totFrameMseLambda);
                    }
                    fprintf(foutFramePsnrBeta, "%lf\n", psnrFrameBeta);
                    fprintf(foutFramePsnrLambda, "%lf\n", psnrFrameLambda);

                    // Free and reset.  
                    freeResultsList(results_list);
                    results_list = NULL;
                    skipCnt = 0;
                }
            }
            else {
                if (inSlice == false) {
                    if (strstr(line, "<Slice") != NULL) {
                        inSlice = true;
                        mbcnt = 0;
                    }
                }
                else {
                    if (strstr(line, "</Slice") != NULL) {	
                        inSlice = false;
                    }
                    else {
                        if (inMB == false) {
                            if (strstr(line, "<MacroBlock") != NULL) {
                                inMB = true;
                                mbcnt += 1;
                            }
                        }
                        else {
                            if (strstr(line, "</MacroBlock") != NULL) {
                                inMB = false;
                                free(predmodestring);
                                predmodestring = NULL;

                                //Free matrix
                                for (i=0; i<COEFFICIENTS; i++){
                                    free(coefficients[i]);
                                }
                                free(coefficients);
                            }
                            else {
                                if (inRow == false) {
                                    if (strstr(line, "<QP_Y>") != NULL) {
                                        char* pos1;
                                        char* pos2;
                                        pos1 = strstr(line, "<QP_Y>");
                                        pos2 = strstr(line, "</QP_Y>");
                                        numValues = pos2 - pos1 - 6;
                                        allocLine = (char*) malloc((numValues+1) * sizeof(char));
                                        i = 0;
                                        while (i < numValues) {
                                            allocLine[i] = pos1[6+i];
                                            i++;
                                        }
                                        allocLine[i] = '\0';
                                        qp = atoi(allocLine);
                                        free(allocLine);
                                        allocLine = NULL;
                                    }
                                    else if (strstr(line, "<Type>") != NULL) {
                                        char* pos1;
                                        char* pos2;
                                        pos1 = strstr(line, "<Type>");
                                        pos2 = strstr(line, "</Type>");
                                        numValues = pos2 - pos1 - 6;
                                        allocLine = (char*) malloc((numValues+1) * sizeof(char));
                                        i = 0;
                                        while (i < numValues) {
                                            allocLine[i] = pos1[6+i];
                                            i++;
                                        }
                                        allocLine[i] = '\0';
                                        typeval = atoi(allocLine);
                                        free(allocLine);
                                        allocLine = NULL;
                                    }
                                    else if (strstr(line, "<PredModeString>") != NULL) {
                                        char* pos1;
                                        char* pos2;
                                        pos1 = strstr(line, "<PredModeString>");
                                        pos2 = strstr(line, "</PredModeString>");
                                        numValues = pos2 - pos1 - 16;
                                        predmodestring = (char*) malloc((numValues+1) * sizeof(char));
                                        i = 0;
                                        while (i < numValues) {
                                            predmodestring[i] = pos1[16+i];
                                            i++;
                                        }
                                        predmodestring[i] = '\0';
                                    }
                                    else if (strstr(line, "<SkipFlag>") != NULL) {
                                        char* pos1;
                                        char* pos2;
                                        pos1 = strstr(line, "<SkipFlag>");
                                        pos2 = strstr(line, "</SkipFlag>");
                                        numValues = pos2 - pos1 - 10;
                                        allocLine = (char*) malloc((numValues+1) * sizeof(char));
                                        i = 0;
                                        while (i < numValues) {
                                            allocLine[i] = pos1[10+i];
                                            i++;
                                        }
                                        allocLine[i] = '\0';
                                        skipflag = atoi(allocLine);
                                        free(allocLine);
                                        allocLine = NULL;
                                    }
                                    else if (strstr(line, "<Plane type=\"0\">") != NULL) {
                                        inRow = true;
                                        rowCounter = 0;
                                    }
                                }
                                // If inside row
                                else {
                                    if (strstr(line, "</Plane>") == NULL) {
                                        char* pos1;
                                        char* pos2;
                                        pos1 = strstr(line, "<Row>");
                                        pos2 = strstr(line, "</Row>");
                                        numValues = pos2 - pos1 - 5;
                                        // Taking the whole row
                                        allocLine = (char*) malloc((numValues+1) * sizeof(char));
                                        i = 0;
                                        while (i < numValues) {
                                            allocLine[i] = pos1[5+i];
                                            i++;
                                        }
                                        allocLine[i] = '\0';

                                        // Create the coefficient matrix 
                                        if (rowCounter == 0) {
                                            coefficients = (int **) malloc (COEFFICIENTS * sizeof(int*));
                                            for (i = 0; i < COEFFICIENTS; i++) {
                                                coefficients[i] = (int *)malloc(COEFFICIENTS* sizeof(int));
                                            }
                                        }

                                        // Tokenize string
                                        char* token;
                                        token = strtok(allocLine, ",");
                                        i = 0;
                                        while (token != NULL) {
                                            coefficients[rowCounter][i] = atoi(token);
                                            token = strtok(NULL, ",");
                                            i++;
                                        }
                                        rowCounter++;

                                        // Free the line
                                        free(allocLine);
                                        allocLine = NULL;
                                    }
                                    else {
                                        boolean allZeroCoeff = 0;
                                        allZeroCoeff = check_all_zeroes(coefficients);

                                        // If coefficients are *not* all zeros
                                        // and skipFlag is *not* set then
                                        // call the "no_psnr_estimation" function
                                        if ((allZeroCoeff == false) && (skipflag == 0)) {
                                            no_psnr_calculation(&results_list, coefficients, qp, predmodestring, typeval);
                                        }
                                        // If one of the previous conditions
                                        // does not hold, call mse_prediction
                                        else {
                                            float skipRate = .0;
                                            skipCnt++;

                                            if (mbcnt > 0) {
                                                skipRate = (float) skipCnt / (float) mbcnt;
                                                // Call the prediction function
                                                mse_prediction(skipRate, &results_list, mseRef); 
                                            }
                                            // TODO: 'else' branch
                                            if (DEBUG) printf("PIC: %d MB: %d skipped, skipflag = 1.\n", pic_id, mbcnt);
                                        }
                                        inRow = false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        free(line);
        line = NULL;
        read = getline(&line, &len, fr);
    }

    // Free       
    free(line);
    fclose(fr);
    fclose(foutMseLambda);
    fclose(foutMseBeta);
    fclose(foutPsnrBeta);
    fclose(foutPsnrLambda);
    fclose(foutFramePsnrBeta);
    fclose(foutFramePsnrLambda);
    fclose(foutFrameMseBeta);
    fclose(foutFrameMseLambda);

    // Return 0
    exit (0);
}

// Calculate average MSE of all picture of both beta and lambda.
// This is done by making an *average* of all the MSE (both beta and lambda)
// calculated for each MB. 
mseRefRes calc_average_mse (resultsList** head) {
    resultsList* current = *head;
    mseRefRes results;
    results.mseLambdaRef = 0;
    results.mseBetaRef = 0;

    float sumBeta = .0, sumLambda = .0;
    int cnt = 0;

    // Calculate sum of all mseBeta and mseLambda.
    while (current != NULL) {
        sumBeta += current->mseBeta;
        sumLambda += current->mseLambda;
        cnt++;
        current = current->next;
    }

    if (cnt > 0) {
        results.mseBetaRef = (double)(sumBeta / cnt);
        results.mseLambdaRef = (double)(sumLambda/ cnt);
    }
    else {
        results.mseBetaRef = -1;
        results.mseLambdaRef= -1;
    }
    return results;
}

/* Function to handle the prediction
 * This function is called when:
 *     - skipFlag is set to 1
 *     - all coefficients are quantized to 0
 * All the predicted values are inserted in the results_list list.
 * Nothing is returned. */
void mse_prediction (float skipRate, resultsList** results_list, mseRefRes mseRef) {

    float sumBeta = .0, sumLambda = .0, averageBeta = .0, averageLambda = .0,
          estimatedMseBeta = .0, estimatedMseLambda = .0, estimatedpsnrBeta = .0, 
          estimatedpsnrLambda = .0;
    int cnt = 0;
    resultsList* current = *results_list;

    // If the first MB is skipped, the results_list is still NULL (nothing has
    // been inserted in the results) or the skipRate is 1 (1 skip / 1 MB).
    // See equation (22) on the paper for this assumption.
    if ((current == NULL) || (skipRate == 1)) {
        estimatedMseBeta = mseRef.mseBetaRef;
        estimatedMseLambda = mseRef.mseLambdaRef;
    }
    else {
        // Calculate average of the preceeding non-skipped MBs.
        while (current != NULL) {
            if (current->predicted == false ) {
                sumBeta += current->mseBeta;
                sumLambda += current->mseLambda;
                cnt++;
            }
            current = current->next;
        }
        if (cnt > 0) {
            averageBeta = sumBeta / cnt;
            averageLambda = sumLambda / cnt;
            estimatedMseBeta = (skipRate * mseRef.mseBetaRef) + ((1 - skipRate) * averageBeta);
            estimatedMseLambda = (skipRate * mseRef.mseLambdaRef) + ((1 - skipRate) * averageLambda);
        }
        else {
            // It means that until now all the values have been predicted!
            // TODO: refactor this!
            estimatedMseBeta = mseRef.mseBetaRef;
            estimatedMseLambda = mseRef.mseLambdaRef;
        }
    }

    // Calculate PSNR starting from predicted MSE.
    // TODO: add the else branch.
    if ((estimatedMseBeta > 0) && (estimatedMseLambda > 0)) {
        estimatedpsnrBeta = 10 * (log10(r(255) / estimatedMseBeta));
        estimatedpsnrLambda = 10 * (log10(r(255) / estimatedMseLambda));
    }

    if (DEBUG) {
        printf("[BETA]-> Estimated MSE %f, estimated PSNR: %f\n", estimatedMseBeta, estimatedpsnrBeta);
        printf("[LAMBDA]-> Estimated MSE: %f, estimated PSNR: %f\n", estimatedMseLambda, estimatedpsnrLambda);
    }

    // Passing "true" as last argument so that the 'predicted' variable is
    // correctly set.
    insertResultsList(results_list, estimatedMseBeta, estimatedMseLambda, estimatedpsnrBeta, estimatedpsnrLambda, true);
}

// Function to check if all quantized coefficients are set to zero.
// Return a boolean value: true for all zeroes, false otherwise. 
boolean check_all_zeroes (int** coefficients) {
    int i = 0, j = 0; 
    for (i=0; i<COEFFICIENTS; i++) {
        for (j=0; j<COEFFICIENTS; j++) {
            if (coefficients[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

/*** NO Reference psnr Estimation ***/
int no_psnr_calculation (resultsList** results_list, int** coefficients, int qp, char* predmodestring, int typeval) {
    uint i = 0, j = 0, totCntr = 0, zeroCntr = 0, oneCntr = 0;
    int qk = 0, zeroCounter = 0, oneCounter = 0;

    double numIntegralFirstBeta = 0.0, numIntegralSecondBeta = 0.0,
           denumIntegralBeta = 0.0, sumEpsilonSquaredBeta = 0.0, 
           cf = 0.0, psnrBeta = 0.0, mseBeta = 0.0, beta = 0.0, 
           psnrLambda = 0.0, mseLambda = 0.0, ak = 0.0, bk = 0.0,
           alpha = 0.0, lambda_laplace = 0.0, 
           numIntegralFirstLambda = 0.0, numIntegralSecondLambda = 0.0, 
           denumIntegralLambda = 0.0, sumEpsilonSquaredLambda = 0.0;
    
    // Structure for both 0 and !0 coefficients. See (7) equation. 
    //  Arrays {a,b}k0 contain the 0 coefficients
    double a0[TOTAL_COEFFICIENTS];
    double b0[TOTAL_COEFFICIENTS];
    // Arrays {a,b}k1 contain the non 0 coefficients. See (7) equation.
    double a1[TOTAL_COEFFICIENTS];
    double b1[TOTAL_COEFFICIENTS];

    // Style print
    /*printf("[-----------------------------------------\n");*/
    if(DEBUG == 1){
        for(i=0; i<COEFFICIENTS; i++){
            printf("[ ");
            for(j=0; j<COEFFICIENTS; j++){
                printf("%d, ", coefficients[i][j]);
            }
            printf("]\n");
        }
    }

    // Check macroblock type and define which ALPHA to use
    alpha = check_macroblock_type(predmodestring);
    if (alpha == 0) return EXIT_FAILURE;

    // Calculate QK starting from given QP
    qk = calculate_qk(qp);

    // Calculating ak and bk for the case where the coefficient is 0 or !0
    for (i = 0; i < COEFFICIENTS; i++) {
        for (j = 0; j < COEFFICIENTS; j++) {
            cf = coefficients[i][j];
            if (cf == 0) {
                a0[zeroCounter] = alpha * qk * (-1);
                b0[zeroCounter] = alpha * qk;
                zeroCounter++;
            }
            else if (cf != 0) {
                a1[oneCounter] = abs(cf) - ((1 - alpha) * qk);
                b1[oneCounter] = abs(cf) + (alpha * qk);
                oneCounter++;
            }
        }
    }

    // Check whether all coefficients are equal to zero
    if (zeroCounter == TOTAL_COEFFICIENTS) {
        /*printf("All coefficients are zero, skipping block\n");*/
        /*printf("-----------------------------------------]\n");*/
        return EXIT_FAILURE;
    }

    // Retrieve beta from the cauchy distribution calculation
    beta = cauchy_distribution(a1, b1, zeroCounter, oneCounter, qk, alpha);
    lambda_laplace = laplace_distribution(b1, zeroCounter, oneCounter, qk, alpha);

    // Print
    if (DEBUG){
      printf("BETA: %.8lf\n", beta);
      printf("LAMBDA: %.8lf\n", lambda_laplace);
    }

    // Find single epsilon (6)
    // Calculating the integral over ak and bk
    for (i = 0; i < COEFFICIENTS; i++) {
        for (j = 0; j < COEFFICIENTS; j++) {
            cf = coefficients[i][j];
            if (cf == 0) {
                // Use all the "0" related arrays and counters
                ak = a0[zeroCntr];
                bk = b0[zeroCntr];
                zeroCntr++;
            }
            else {
                ak = a1[oneCntr];
                bk = b1[oneCntr];
                oneCntr++;
            }

            cf = fabs(cf);

            // Calculating integral 
            numIntegralFirstBeta = (beta * ( bk - cf * log(r(beta) + r(bk))) + (r(cf) - r(beta)) * atan(bk / beta)) / M_PI;
            numIntegralSecondBeta = (beta * ( ak - cf * log(r(beta) + r(ak))) + (r(cf) - r(beta)) * atan(ak / beta)) / M_PI; 
            denumIntegralBeta = (atan(bk / beta) - (atan(ak / beta))) / M_PI;
            sumEpsilonSquaredBeta += ((numIntegralFirstBeta - numIntegralSecondBeta) / denumIntegralBeta);

            //Using LAMBDA from the LaPlace model
            numIntegralFirstLambda = exp(lambda_laplace* (-1) * bk) * (( lambda_laplace * (cf - bk) * ((lambda_laplace* bk) - (lambda_laplace * cf) + 2)) - 2);
            numIntegralSecondLambda = exp(lambda_laplace* (-1) * ak) * (( lambda_laplace * (cf - ak) * ((lambda_laplace* ak) - (lambda_laplace * cf) + 2)) - 2);
            denumIntegralLambda = 2 * r(lambda_laplace);
            sumEpsilonSquaredLambda += (numIntegralFirstLambda - numIntegralSecondLambda) / denumIntegralLambda;

            totCntr++;
        }
    }

    // MSEBeta (4)
    mseBeta = sumEpsilonSquaredBeta / totCntr;
    if (DEBUG) printf("MSE(Beta): %.8lf\n", mseBeta);

    // MSELambda (4)
    mseLambda = sumEpsilonSquaredLambda / totCntr;
    if (DEBUG) printf("MSE(Lambda): %.8lf\n", mseLambda);

    // PSNRBeta (4)
    psnrBeta = 10 * (log10(r(255) / mseBeta));
    if (DEBUG) printf("PSNR(Beta): %.8lf\n", psnrBeta);

    // PSNRLambda (4)
    psnrLambda = 10 * (log10(r(255) / mseLambda));
    if (DEBUG) printf("PSNR(Lambda): %.8lf\n", psnrLambda);

    // Insert results in structure 
    insertResultsList(results_list, mseBeta, mseLambda, psnrBeta, psnrLambda, false);

    /*printf("-----------------------------------------]\n");*/
    return EXIT_SUCCESS;
}

/***
 * Calculate the Cauchy Distribution using the Newton-Raphson's iterative method
 * Starting beta value defined below
 */
double cauchy_distribution(double *a1, double *b1, int zeroCounter, int oneCounter, int qk, double alpha) {
    double beta = 0.0, newbeta = 0.0, res = 0.0;
    int i = 0;

    // Start the beta approximation
    beta = BETA_START;

    for (i = 0; i < LIMIT_NEWTON_ITERATIONS; i++) {
        // x(n+1) = x(n) - (f(xn) / df(xn))
        newbeta = beta - (f(a1, b1, beta, zeroCounter, oneCounter, qk, alpha) / df(a1, b1, beta, zeroCounter, oneCounter, qk, alpha));

        if (DEBUG) printf("Current beta after %d approximations is: %lf\n", i, beta);

        res = newbeta - beta;
        if (fabs(res) < 0.001f) {
            return newbeta;
        }
        else {
            beta = newbeta;
        }
    }

    if (DEBUG) printf("Beta approximation after %d iterations is: %lf\n", i, beta);

    return beta;
}

// Function f(x0) for cauchy_distribution
double f(double *a1, double *b1, double B, int N0, int N1, int qk, double alpha) {

    double firstSum = 0, secondSum = 0, num = 0, denum = 0;
    int k0 = 0, k1 = 0;

    // Calculating the two sums in different loops
    for (k1 = 0; k1 < N1; k1++) {
        num = (a1[k1] / (r(B) + r(a1[k1]))) - (b1[k1] / (r(B) + r(b1[k1])));
        denum = atan(b1[k1] / B) - atan(a1[k1] / B);
        firstSum += (num / denum);
    }

    for (k0 = 0; k0 < N0; k0++) {
        num = (alpha * qk);
        denum = (atan((alpha * qk) / B)) * (r(alpha * qk) + r(B));
        secondSum += (num / denum);
    }

    return firstSum - secondSum;
}

// Function f'(x0) for cauchy
double df(double *a1, double *b1, double B, int N0, int N1, int qk, double alpha) {

    double firstSum = 0, secondSum = 0;
    double num = 0, num1 = 0, num2 = 0, num3 = 0, denum = 0;
    int k0 = 0, k1 = 0;

    for (k1 = 0; k1 < N1; k1++) {
        num1 = ((2 * b1[k1] * B) / r(r(b1[k1]) + r(B))) - ((2 * a1[k1] * B) / r(r(a1[k1]) + r(B)));
        num2 = atan(b1[k1] / B) - atan(a1[k1] / B);
        num3 = (r(a1[k1] - b1[k1]) * r(r(B) - (a1[k1] * b1[k1]))) / (r(r(a1[k1]) + r(B)) * r(r(b1[k1]) + r(B)));
        denum = r(atan(a1[k1] / B) - atan(b1[k1] / B));
        firstSum += ((num1 * num2) - num3) / denum;
    }

    for (k0 = 0; k0 < N0; k0++) {
        num = (alpha * qk);
        denum = (atan((alpha * qk) / B)) * (r(alpha * qk) + r(B));
        secondSum += (num / denum);
    }

    return firstSum - secondSum;
}

double r(double number) {
    return pow(number, 2);
}

// Function for lambda approximation using laplace distribution
double laplace_distribution(double *b1, int zeroCounter, int oneCounter, int qk, double alpha) {

    double lambda_laplace = 0, new_lambda_laplace = 0, res = 0;
    int i = 0;

    lambda_laplace = LAMBDA_START;

    for (i = 0; i < LIMIT_NEWTON_ITERATIONS; i++){
        // x(n+1) = x(n) - (f(xn) / df(xn))
        new_lambda_laplace = lambda_laplace -  (lpf(b1, lambda_laplace, zeroCounter, oneCounter, qk, alpha) / lpdf(lambda_laplace, zeroCounter, oneCounter, qk, alpha));

        res = new_lambda_laplace - lambda_laplace;
        if( fabs(res) < 0.00001f){
            return new_lambda_laplace;
        }
        else{
            lambda_laplace = new_lambda_laplace;
        }

    }
    if (DEBUG) printf("Lambda approximation after %d iterations is: %lf\n", i, lambda_laplace);

    return lambda_laplace;
}

// Function f(x0) for laplace_distribution 
double lpf(double *b1, double lambda_laplace, int N0, int N1, int qk, double alpha) {

    double firstSum = 0.0, secondSum = 0.0;
    int k0 = 0, k1 = 0;

    // Calculating the two sums in different loops
    for (k1 = 0; k1 < N1; k1++) {
        firstSum += (qk / (exp(lambda_laplace * qk) -1 )) - b1[k1];
    }

    for (k0 = 0; k0 < N0; k0++) {
        secondSum += (alpha * qk)/(exp(alpha * lambda_laplace * qk) -1);
    }

    return firstSum + secondSum;
}

// Function f'(x0) for laplace 
double lpdf(double lambda_laplace, int N0, int N1, int qk, double alpha) {

    double firstSum = 0.0, secondSum = 0.0;
    double num = 0.0, denum = 0.0, coseno = 0.0;
    int k0 = 0, k1 = 0;

    for (k0 = 0; k0 < N1; k0++) {
        num = r(alpha) * r(qk);
        coseno = cosh(alpha* lambda_laplace * qk);
        denum = 2 - (2 * coseno);
        firstSum += ( num / denum );
    }

    for (k1 = 0; k1 < N0; k1++) {
        num = r(qk) * exp(qk * lambda_laplace);
        denum = r(exp(qk * lambda_laplace) -1);
        secondSum += (num / denum) * (-1);
    }

    return firstSum + secondSum;
}

// Function to check the macroblock type and return the proper alpha
// Check MacroBlock type Check Tables 7-11, 7-12, 7-13, 7-14 of ITU H264
double check_macroblock_type(char* predmodestring) {

    double alpha = 0;
    
    if (strcmp(predmodestring, "B_SKIP") == 0) {
        // If it's SKIP I can skip the calculation 
        alpha = 0;
    }
    else if (strcmp(predmodestring, "BLOCK_TYPE_I") == 0 ) {
        // Type P can have only Intra frames
        alpha = ALPHA_INTRA; 
    } 
    else {
        // Treating P and B types as INTERCODED frames.
        if (strcmp(predmodestring , "BLOCK_TYPE_B") == 0) {
            alpha = ALPHA_INTER; 
        } 
        else if ( strcmp(predmodestring , "BLOCK_TYPE_P") == 0 ) {
            alpha = ALPHA_INTER; 
        } 
    }

    return alpha;
}

// Function to calculate QK starting from QP
double calculate_qk(int qp) {
    int mod_return = 0;
    float qb = 0.0;

    mod_return = qp % 6;
    switch (mod_return){
        case 0:
            qb = 0.6250;
            break;
        case 1:
            qb = 0.6875;
            break;
        case 2:
            qb = 0.8125;
            break;
        case 3:
            qb = 0.8750;
            break;
        case 4:
            qb = 1.0000;
            break;
        case 5:
            qb = 1.1250;
            break;
        default:
            break;
    }

    return qb * pow(2, qp / 6);
}

// Insert in list.
void insertResultsList(resultsList** head, double mseBeta, double mseLambda, double psnrBeta, double psnrLambda, boolean predicted) {
    if (*head == NULL) {
        resultsList* newNode = malloc(sizeof(resultsList));
        newNode->mseBeta = mseBeta;
        newNode->mseLambda = mseLambda;
        newNode->psnrBeta = psnrBeta;
        newNode->psnrLambda = psnrLambda;
        newNode->predicted = predicted;
        newNode->next = NULL;
        *head = newNode;
    }
    else {
        resultsList* tmp = *head;
        while(tmp->next != NULL){
                tmp = tmp->next;
        }
        tmp->next = malloc(sizeof(resultsList));
        tmp->next->mseBeta = mseBeta;
        tmp->next->mseLambda = mseLambda;
        tmp->next->psnrBeta = psnrBeta;
        tmp->next->psnrLambda = psnrLambda;
        tmp->next->predicted = predicted;
        tmp->next->next = NULL;
    }
}

// Free the list.
void freeResultsList(resultsList* head){
    while (head != NULL){
        resultsList* toFree = head;
        head = head->next;
        free(toFree);
    }
}
