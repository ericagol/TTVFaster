#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ttvfaster.h"


int main(int argc, char **argv)
{
  int i,k,n_planets,npar,m_max,n_events;
  double t0,tf;
  FILE *dynam_param_file; //input param file
  FILE *output_file; //output file for transit times
  CalcTransit* model;

  if(argc!=7){
    printf("Incorrect # of input params: ./predict_formula N_planets Param_File T0 TF J_Max Output_File \n");
    exit(-1);
  }
  sscanf(argv[1],"%d",&n_planets);
  npar = 1+7*n_planets; //dynamics parameters.

  double params[npar]; //dynamics parameter array
  dynam_param_file = fopen(argv[2],"r");
  for(i =0;i<npar;i++){
    fscanf(dynam_param_file, "%lf ",&params[i]);
  }
  fclose(dynam_param_file);
  sscanf(argv[3],"%lf",&t0);
  sscanf(argv[4],"%lf",&tf);
  sscanf(argv[5],"%d",&m_max);

  model = ttvfaster(n_planets, params, t0, tf, m_max, &n_events);

  output_file = fopen(argv[6],"w");
  for(k=0;k<n_events;k++){
    fprintf(output_file,"%d %d %.16le\n",(model+k)->planet,(model+k)->epoch,(model+k)->time);
  }
  fflush(output_file);
  fclose(output_file);

  free(model);
  return 0;
}
