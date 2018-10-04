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

  int C1[L1_length-1][L2_length-1];

  int C2[L1_length-1][L2_length-1];
  std::vector<int> ii;
  std::vector<int> jj;
  
  ii.reserve((L2_length-1)*(L1_length-1));
  jj.reserve((L2_length-1)*(L1_length-1));


     double temp1[L1_length-1][L2_length];
     double temp2[L2_length-1][L1_length];
     double D1[L1_length-1][L2_length-1];
     double D2[L2_length-1][L1_length-1];

     
  // Diff x and Diff y for L1
  for (int i=0; i<L1_length-1; i++){
    for (int j=0; j<L2_length-1; j++){


     dx1[i] = L1[2*(i+1)] - L1[2*i];
     dy1[i] = L1[1 + 2*(i+1)] - L1[1 + 2*i];
     dx2[j] = L2[2*(j+1)] - L2[2*j];
     dy2[j] = L2[1 + 2*(j+1)] - L2[1 + 2*j];
     S1[i] = dx1[i]*L1[1 + 2*i] - dy1[i]*L1[2*i];
     S2[j] = dx2[j]*L2[1 + 2*j] - dy2[j]*L2[2*j];
     temp1[i][j] = dx1[i]*L2[1+2*j] - dy1[i]*L2[2*j];
     temp2[j][i] = dx2[j]*L1[1+2*i] - dy2[j]*L1[2*i];
     D1[i][j] = (temp1[i][j]-S1[i])*(temp1[i][j+1]-S1[i]);
     C1[i][j] = D1[i][j] <= 0;
     D2[j][i] = (temp2[j][i]-S2[j])*(temp2[j][i+1]-S2[j]);
     C2[i][j] = (int) D2[j][i] <= 0;
     

      int test;
      test = ((C1[i][j] == 1) && (C2[i][j] == 1));
      mexPrintf("test: %d", test);

           P_out_m = plhs[0] = mxCreateDoubleMatrix(2, 0, mxREAL); return;

      if ((C1[i][j] == 1) && (C2[i][j] == 1)){

        ii.push_back(i);
        jj.push_back(j);
      }
    }
    // mexPrintf("diffx1 %f \n", dx1[i]);
  }


  // Diff x and Diff y for L2
  // for (int j=0; j<L2_length-1; j++){
  //   dx2[j] = L2[2*(j+1)] - L2[2*j];
  //   dy2[j] = L2[1 + 2*(j+1)] - L2[1 + 2*j];
  // }

  // Determine 'signed distances'
 

  // for (int i=0; i<L1_length-1; i++){
  //     S1[i] = dx1[i]*L1[1 + 2*i] - dy1[i]*L1[2*i];
  //     // mexPrintf("S1[%d]: %f \n", i, S1[i]);
  // }

  // for (int j=0; j<L2_length-1; j++){
  //     S2[j] = dx2[j]*L2[1 + 2*j] - dy2[j]*L2[2*j];
  //     // mexPrintf("S2[%d]: %f \n", i, S2[i]);
  // }

  // Some mysterious temporary variables
  // for (int j=0; j<L2_length; j++){
  //   for (int i=0; i<L1_length-1; i++){
  //     temp1[i][j] = dx1[i]*L2[1+2*j] - dy1[i]*L2[2*j];
  //     // mexPrintf("temp1[%d][%d]: %f \n", i, j, temp1[i][j]);
  //   }
  // }

  //////////// This is making Matlab crash, probably because temp1 is too big.
  ////////////   see if you can combine all of these double for loops so we don't
  ////////////   need all of these big temporary variables.
  // mexPrintf("%f", temp1[0][0]);
  ////////////
  
  // for (int i=0; i<L1_length; i++){
  //   for (int j=0; j<L2_length-1; j++){
  //     temp2[j][i] = dx2[j]*L1[1+2*i] - dy2[j]*L1[2*i];
  //     // mexPrintf("temp2[%d][%d]: %f \n", i, j, temp2[i][j]);
  //   }
  // }


  // mexPrintf("D1: %dx%d\n", L1_length-1, L2_length-1);



  // for (int j = 0; j < L2_length-1; j++){
  //   for (int i=0; i<L1_length-1; i++){
  //     D1[i][j] = (temp1[i][j]-S1[i])*(temp1[i][j+1]-S1[i]);
  //     C1[i][j] = D1[i][j] <= 0;
  //     // mexPrintf("first part [%d][%d]: %f \n", i, j, (temp1[i][j]-S1[i]));
  //     // mexPrintf("second part [%d][%d]: %f \n", i, j, (temp1[i][j+1]-S1[i]));
  //     // mexPrintf("D1[%d][%d]: %f \n", i, j, D1[i][j]);
  //     // mexPrintf("C1[%d][%d]: %d \n", i, j, C1[i][j]);
  //   }
  // }
  


  // mexPrintf("D2: %dx%d\n", L2_length-1, L1_length-1);

  // for (int i = 0; i < L1_length-1; i++){
  //   for (int j=0; j<L2_length-1; j++){
  //     D2[j][i] = (temp2[j][i]-S2[j])*(temp2[j][j+1]-S2[j]);
  //     C2[i][j] = (int) D2[j][i] <= 0;
  //     // mexPrintf("first part [%d][%d]: %f \n", i, j, (temp2[i][j]-S2[i]));
  //     // mexPrintf("second part [%d][%d]: %f \n", i, j, (temp1[i][j+1]-S1[i]));
  //     // mexPrintf("D2[%d][%d]: %f \n", i, j, D2[i][j]);
  //     // mexPrintf("C2[%d][%d]: %d \n", j, i, C2[j][i]);
  //   }
  // }

  // find where C1 & C2



  // for (int j = 0; j < L2_length - 1; j++){
  //   for (int i = 0; i < L1_length - 1; i++){
  //     // usleep(1);
  //     // mexPrintf("Checking index: %d, %d\n", i, j);

  //       // mexPrintf("Adding index: %d, %d\n", i, j);
  //     }
  //   }
  // }

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
