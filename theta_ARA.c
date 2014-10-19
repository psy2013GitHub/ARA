








#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/*
  prhs[0], dat

*/
  double *dat, *serial_num, *output;
  int i=0, j=0, k=0, nvox=0, cumul_noxs=0, npoint=0, tmp_point=0;
  double tmp_value=0, avg=0, sum1=0, sum2=0, theta=0;
  
  //check input/output error
  if(nlhs!=1||nrhs!=2)
  	mexErrMsgTxt("Usage: ara=test(dat,serial_num)");

  dat=mxGetPr(prhs[0]);
  serial_num=mxGetPr(prhs[1]);
  nvox=mxGetN(prhs[1]);
  plhs[0]=mxCreateDoubleMatrix(nvox,1,mxREAL);
  output=mxGetPr(plhs[0]); // care !!!
  npoint=mxGetM(prhs[0]);  // dimesion
  while(i++<nvox)
  {
    tmp_point=*serial_num++;
    sum1=0; sum2=0;
    for(j=0;j<tmp_point&&k<npoint;j++,k++)
    	{
    		tmp_value=*dat++;
    		sum1+=tmp_value;
    		sum2+=tmp_value*tmp_value;
    	}
    *output++=sum2/(double)(tmp_point)-sum1*sum1/(double)(tmp_point*tmp_point);
  }
}
