#include <assert.h>
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>
#include <unistd.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /*
   * Error checking
   */
  if (nlhs == 0){
    mexPrintf("Error: output not assigned\n");
    return;
  }else if(nlhs != 1){
    mexPrintf("Error: there is only one return value\n");
    return;
  }

  if (nrhs != 2){
    mexPrintf("Error: must be two inputs\n");
    return;
  }

  //declare size variables
  const mwSize *dims_L1, *dims_L2;
  int L1_length, L2_length, numdims_L1, numdims_L2;
  
  //figure out dimensions
  dims_L1 = mxGetDimensions(prhs[0]);
  dims_L2 = mxGetDimensions(prhs[1]);
  numdims_L1 = mxGetNumberOfDimensions(prhs[0]);
  numdims_L2 = mxGetNumberOfDimensions(prhs[1]);

  // check that the inputs are matrices
  if (numdims_L1 != 2 || numdims_L2 != 2){
    mexPrintf("Error: Inputs must be matrices\n");  
  }

  // check that there are only two rows
  if ((int) dims_L1[0] != 2 || (int) dims_L2[0] != 2){
    mexPrintf("Error: Inputs must have 2 rows\n");
  }

  /*
   * Start doing real work
   */

  // declare data variables
  mxArray *L1_in_m, *L2_in_m, *P_out_m;

  //associate inputs
  L1_in_m = mxDuplicateArray(prhs[0]);
  L2_in_m = mxDuplicateArray(prhs[1]);

  L1_length = (int)dims_L1[1];
  L2_length = (int)dims_L2[1];

  
  // get normal pointers
  double *L1, *L2, *P;
  L1 = mxGetPr(L1_in_m);
  L2 = mxGetPr(L2_in_m);

  // compute 
  double dx1[L1_length-1];
  double dy1[L1_length-1];

  double dx2[L2_length-1];
  double dy2[L2_length-1];

  double S1[L1_length-1];
  double S2[L2_length-1];

  // int C1[L1_length-1][L2_length-1];

  // int C2[L1_length-1][L2_length-1];
  std::vector<int> ii;
  std::vector<int> jj;
  
  ii.reserve((L2_length-1)*(L1_length-1));
  jj.reserve((L2_length-1)*(L1_length-1));


     // double temp1[L1_length-1][L2_length];
     // double temp2[L2_length-1][L1_length];
     // double D1[L1_length-1][L2_length-1];
     // double D2[L2_length-1][L1_length-1];

 double temp1, temp1_2, temp2, temp2_2, D1, D2;
 int C1, C2;


  for (int i=0; i<L1_length-1; i++){

    dx1[i] = L1[2*(i+1)] - L1[2*i];
    dy1[i] = L1[1 + 2*(i+1)] - L1[1 + 2*i];
    S1[i] = dx1[i]*L1[1 + 2*i] - dy1[i]*L1[2*i];
  }

  for (int j=0; j<L2_length-1; j++){

    dx2[j] = L2[2*(j+1)] - L2[2*j];
    dy2[j] = L2[1 + 2*(j+1)] - L2[1 + 2*j];
    S2[j] = dx2[j]*L2[1 + 2*j] - dy2[j]*L2[2*j];
  }

  // Diff x and Diff y for L1
  for (int i=0; i<L1_length-1; i++){
    for (int j=0; j<L2_length-1; j++){

     temp1 = dx1[i]*L2[1+2*j] - dy1[i]*L2[2*j];
     temp1_2 = dx1[i]*L2[1+2*(j+1)] - dy1[i]*L2[2*(j+1)];
     temp2 = dx2[j]*L1[1+2*i] - dy2[j]*L1[2*i];
     temp2_2 = dx2[j]*L1[1+2*(i+1)] - dy2[j]*L1[2*(i+1)];
     D1 = (temp1-S1[i])*(temp1_2-S1[i]);
     C1 = D1 <= 0;
     D2 = (temp2-S2[j])*(temp2_2-S2[j]);
     C2 = (int) D2 <= 0;
     if ((C1 ==1) && (C2 ==1)){
        ii.push_back(i);
        jj.push_back(j);
     }
     

      // int test;
      // test = ((C1[i][j] == 1) && (C2[i][j] == 1));
      // mexPrintf("test: %d", test);

    }
  }
  // P_out_m = plhs[0] = mxCreateDoubleMatrix(2, 0, mxREAL); return;


  //associate outputs
  if (ii.size() == 0){
    P_out_m = plhs[0] = mxCreateDoubleMatrix(2, 0, mxREAL);
    return;
  }

  // Prepare for output
  double L[ii.size()];
  int num_zeros = 0;

  for (int i = 0; i<ii.size(); i++){
    L[i] = dy2[jj[i]]*dx1[ii[i]] - dy1[ii[i]]*dx2[jj[i]];
    if (L[i] == 0){
      num_zeros++;
    }
    // mexPrintf("L[%d]: %f", i, L[i]);
  }  

  // associate outputs
  P_out_m = plhs[0] = mxCreateDoubleMatrix(2, ii.size() - num_zeros, mxREAL);
  P = mxGetPr(P_out_m);

  // solve system of equations
  int P_idx = 0;
  for (int i=0; i<ii.size(); i++){
    // avoid divisions by 0
    if (L[i] == 0) {
      continue;
    }
    // solve system of equations to get the common points
    P[2*P_idx] = (dx2[jj[i]]*S1[ii[i]] - dx1[ii[i]]*S2[jj[i]])/L[i];// x coordinate
    P[1 + 2*P_idx] = (dy2[jj[i]]*S1[ii[i]] - dy1[ii[i]]*S2[jj[i]])/L[i];// y coordinate
    P_idx++;
  }

}
