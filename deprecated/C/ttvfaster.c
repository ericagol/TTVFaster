/*Many thanks for Jack Wisdom at MIT for sharing the code "laplace()" below with us. */
/* the longitude of transit = 0 */
/* Please cite Agol & Deck (2015) if you make use of this code in your research.*/
#define LAPLACE_EPS 1.0e-10
#define PI 3.14159265358979323846
#define TWOPI 6.283185307179586476925287
#define MAX_N_PLANETS 9
#define G 2.95994511e-4; /* hey! the units of G are fixed here, to be Msun, Day, AU*/
#include <math.h>
#include <stdlib.h>

#include "ttvfaster.h"

double AJ00,AJ01,AJ02,AJ10,AJ20,AJ11;
double principal_value(double theta);
double laplace(double s, int i, int j, double a);
double gam(int i, int k, int j,double kappa, double beta, double alpha);
double  xi(int i, int k, int j,double kappa, double beta, double alpha);
double C1(int i, int k, int j, double kappa, double beta, double alpha);
double C2(int i, int k, int j, double kappa, double beta, double alpha);
double d1(int i, int k, int j, double kappa, double beta, double alpha);
double d2(int i, int k, int j, double kappa, double beta, double alpha);
double f_j_k_i(double alpha,int j,int k, int i,double kappa_o_m, double beta_o_m,double *a00,double *a01, double *a10, double *a20, double *a02,double *a11);
double u(double gam, double C1v, double C2v);
double v_plus(double z, double d1v, double d2v);
double v_minus(double z, double d1v, double d2v);


CalcTransit* ttvfaster (
    // Inputs:
    int n_planets,
    double* params,
    double t0,
    double tf,
    int m_max,

    // Outputs:
    int* n_events_out
) {
  int i,k,j,count,m;
  double P1,P2,TT1,TT2,n1,n2,m1,m2,mstar,e1,e2,ap1,ap2,trueAnom,EccAnom,MT1,MT2,psi,fac1,fac2,base,TTV,alpha,lambda1,lambda2,Omega1,Omega2;
  double beta_over_m,kappa_over_m;
  int n_transits[n_planets];
  int start[n_planets];
  int n_events = 0;
  CalcTransit *model;

  double b[m_max+2],db[m_max+2],d2b[m_max+2];
  double AJ00_arr[m_max+2],AJ10_arr[m_max+2],AJ01_arr[m_max+2],AJ02_arr[m_max+2],AJ20_arr[m_max+2],AJ11_arr[m_max+2];

  start[0] = 0;
  for(j=0;j<n_planets;j++){
    P1 = params[j*7+2];
    n_transits[j]=(int)((tf-t0)/P1+1.0);
    n_events+=n_transits[j];
    if(j>0){
      start[j] = n_transits[j-1]+start[j-1];
    }
  }
  *n_events_out = n_events;

  model = malloc((sizeof(CalcTransit)*n_events));
  if(model==0){exit(1);}
  mstar = params[0];
  /*"Keplerian" fit to times */
  for(j=0;j<n_planets;j++){
    P1 = params[j*7+2];
    TT1 = params[j*7+7];
    for(count = 0;count<n_transits[j];count++){
      k =count+start[j];
      (model+k)->planet = j;
      (model+k)->epoch = count;
      base =P1*count+TT1;
      (model+k)->time = base+t0;

    }
  }
  /* Now add TTVs; adjacent pairs only*/
  for(j=0;j<n_planets-1;j++){
    i = j+1;
    P1 = params[j*7+2];
    TT1 = params[j*7+7];
    m1 = params[j*7+1];
    e1 = sqrt(params[j*7+6]*params[7*j+6]+params[j*7+3]*params[7*j+3]);
    ap1 = atan2(params[j*7+6],params[j*7+3]);
    ap1 = principal_value(ap1);
    Omega1 = principal_value(params[j*7+5]);
    trueAnom = principal_value(0.0-ap1-Omega1);/* at transit*/
    EccAnom =  2.0*atan(sqrt(1.0-e1)/sqrt(1.0+e1)*tan(trueAnom/2.0));
    MT1= EccAnom-e1*sin(EccAnom);
    n1 = 2.0*PI/P1;

    P2 = params[i*7+2];
    TT2 = params[i*7+7];
    m2 = params[i*7+1];
    e2 = sqrt(params[i*7+6]*params[7*i+6]+params[i*7+3]*params[7*i+3]);
    ap2 = atan2(params[i*7+6],params[i*7+3]);
    ap2 = principal_value(ap2);
    Omega2 = principal_value(params[i*7+5]);
    trueAnom = principal_value(0.0-ap2-Omega2); /* at transit*/
    EccAnom =  2.0*atan(sqrt(1.0-e2)/sqrt(1.0+e2)*tan(trueAnom/2.0));
    MT2= EccAnom-e2*sin(EccAnom);
    n2 = 2.0*PI/P2;

    alpha = pow(P1/P2,2.0/3.0);

    for(count=0;count<m_max+2;count++){
      b[count] = laplace(0.5,count,0,alpha);
      db[count] = laplace(0.5,count,1,alpha);
      d2b[count] = laplace(0.5,count,2,alpha);
    }

    kappa_over_m = (pow(alpha,-1.5)-1.0);
    beta_over_m = (1.0-pow(alpha,1.5));
    for(m=0;m<m_max+2;m++){
      AJ00_arr[m] = b[m];
      AJ10_arr[m] = db[m];
      AJ20_arr[m] = d2b[m];
      AJ01_arr[m] = -db[m]-b[m];
      AJ02_arr[m] = 4.0*db[m]+2.0*b[m]+d2b[m];
      AJ11_arr[m] = -2.0*db[m]-d2b[m];
    }


    /*Inner Planet of Pair First */
    for(count = 0;count<n_transits[j];count++){
      k =count+start[j];
      (model+k)->planet = j;
      (model+k)->epoch = count;
      base =P1*count+TT1;
      lambda1 = n1*(base-TT1)+MT1+ap1+Omega1;
      lambda2 = n2*(base-TT2)+MT2+ap2+Omega2;
      psi = lambda1-lambda2;
      fac1 = lambda1-ap1-Omega1;
      fac2 = lambda1-ap2-Omega2;
      TTV=0;
      for(m=1;m<=m_max;m++){


        TTV+=f_j_k_i(alpha,m,0,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi);
        TTV +=e1*f_j_k_i(alpha,m,1,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi+fac1);
        TTV +=e1*f_j_k_i(alpha,m,-1,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi-fac1);
        TTV +=e2*f_j_k_i(alpha,m+1,2,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi+fac2);
        TTV +=e2*f_j_k_i(alpha,m-1,-2,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi-fac2);
      }

      TTV *=m2/mstar/n1;
      (model+k)->time +=TTV;
    }
    /* Outer Planet of Pair*/
    for(count = 0;count<n_transits[i];count++){
      k =count+start[i];
      (model+k)->planet = i;
      (model+k)->epoch = count;
      base =P2*count+TT2;

      lambda1 = n1*(base-TT1)+MT1+ap1+Omega1;
      lambda2 = n2*(base-TT2)+MT2+ap2+Omega2;
      psi = lambda1-lambda2;
      fac1 = lambda2-ap1-Omega1;
      fac2 = lambda2-ap2-Omega2;
      TTV=0.0;

      for(m=1;m<=m_max;m++){
        TTV+=f_j_k_i(alpha,m,0,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi);
        TTV +=e1*f_j_k_i(alpha,m-1,1,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi+fac1);
        TTV +=e1*f_j_k_i(alpha,m+1,-1,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi-fac1);
        TTV +=e2*f_j_k_i(alpha,m,2,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi+fac2);
        TTV +=e2*f_j_k_i(alpha,m,-2,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(m*psi-fac2);

      }

      TTV *=m1/mstar/n2;
      (model+k)->time +=TTV;

    }
  }

  return model;
}


/* Code due to Jack Wisdom */
/* compute Laplace coefficients and Leverrier derivatives
          j
     j   d     i
    a   ---   b (a)
          j    s
        da

   by series summation */



double laplace(double s, int i, int j, double a)
{
  double as, term, sum, factor1, factor2, factor3, factor4;
  int k,q, q0;

  as = a*a;

  if(i<0) i = -i;

  if(j<=i)     /* compute first term in sum */
    {
      factor4 = 1.0;
      for(k=0; k<j; k++)
        factor4 *= (i - k);
      sum = factor4;
      q0=0;
    }
  else
    {
       q0 = (j + 1 - i) / 2;
      sum = 0.0;
      factor4 = 1.0;
    }

  /* compute factors for terms in sum */

  factor1 = s;
  factor2 = s + i;
  factor3 = i + 1.0;
  for(q=1;q<q0;q++)   /* no contribution for q = 0 */
    {
      factor1 *= s + q;
      factor2 *= s + i + q;
      factor3 *= i + 1.0 + q;
    }

  term = as * factor1 * factor2 / (factor3 * q);

  /* sum series */

  while(term*factor4 > LAPLACE_EPS)
    {
      factor4 = 1.0;
      for(k=0;k<j;k++)
        factor4 *= (2*q + i - k);
      sum += term * factor4;
      factor1 += 1.0;
      factor2 += 1.0;
      factor3 += 1.0;
      q++;
      term *= as * factor1 * factor2 / (factor3 * q);

    }

  /* fix coefficient */

  for(k=0;k<i;k++)
    sum *= (s + ((double) k))/(((double) k)+1.0);

  if(q0 <= 0)
    sum *= 2.0 * pow(a, ((double) i));
  else
    sum *= 2.0 * pow(a, ((double) 2*q0 + i - 2));

  return(sum);
}


double gam(int i, int k, int j,double kappa, double beta, double alpha)
{

  double val =0.0;
  if(i==1){
    if(k ==0){
      val = beta;
    }
    if(k == 1){
      val = beta+1.0;
    }
    if(k == -1){
      val = beta-1.0;
    }
    if(k == 2){
      val = beta+pow(alpha,1.5);;
    }
    if(k == -2){
      val = beta-pow(alpha,1.5);
    }
  }else{
 if(k ==0){
      val = kappa;
    }
    if(k == 2){
      val = kappa+1.0;
    }
    if(k == -2){
      val = kappa-1.0;
    }
    if(k == 1){
      val = kappa+pow(alpha,-1.5);;
    }
    if(k == -1){
      val = kappa-pow(alpha,-1.5);
    }
  }
    return(val);
}


double xi(int i, int k, int j,double kappa, double beta, double alpha)
{
  double val=0.0;
  if(i==1){
    val = beta;
  }else{
    val = kappa;
  }
  return(val);
}

double d1(int i, int k, int j, double kappa, double beta, double alpha)
{
  double  val = 0.0;
  if(i==1){
    if(j==1){
      val = alpha*j*(AJ00-alpha);
    }else{
      val = alpha*j*AJ00;
    }
  }else{
    if(j==1){
      val = -j*(AJ00-pow(alpha,-2.0));
    }else{
      val = -j*AJ00;
    }
  }

 return(val);
}

double d2(int i, int k, int j, double kappa, double beta, double alpha)
{
  double  val = 0.0;
  if(i==1){
    if(j==1){
      val = alpha*(AJ10-alpha);
    }else{
      val = alpha*AJ10;
    }
  }else{
    if(j==1){
      val = AJ01-pow(alpha,-2.0);
    }else{
      val = AJ01;
    }
  }

 return(val);
}


double C1(int i, int k, int j, double kappa, double beta, double alpha)
{

  double val=0.0;
  if(i==1){
    if(k==0){
      val = d1(i,k,j,kappa,beta,alpha);
    }
    if(k ==1){
      if(j==1){
        val = alpha*j*(j*AJ00-0.5*AJ10+0.5*(1.0-2.0)*alpha);
      }else{
        val = alpha*j*(j*AJ00-0.5*AJ10);
      }
    }
    if(k ==-1){
      if(j==1){
        val = alpha*j*(-j*AJ00-0.5*AJ10+0.5*(1.0+2.0)*alpha);
      }else{
        val = alpha*j*(-j*AJ00-0.5*AJ10);
      }
    }


    if(k ==2){
      val = alpha*j*(-j*AJ00-0.5*AJ01);
    }

    if(k ==-2){
      if(j==1){
        val = alpha*j*(j*AJ00-0.5*AJ01-2.0*alpha);
      }else{
        val = alpha*j*(j*AJ00-0.5*AJ01);
      }
    }
  }else{
/* Outer Planet */
    if(k==0){
      val = d1(i,k,j,kappa,beta,alpha);
    }
    if(k ==2){
      if(j==1){
        val = -j*(-j*AJ00-0.5*AJ01+0.5*3.0*pow(alpha,-2.0));
      }else{
        val = -j*(-j*AJ00-0.5*AJ01);
      }
    }
    if(k ==-2){
      if(j==1){
        val = -j*(j*AJ00-0.5*AJ01+0.5*(1.0-2.0)*pow(alpha,-2.0));
      }else{
        val = -j*(j*AJ00-0.5*AJ01);
      }
    }


    if(k ==-1){
      val = -j*(-j*AJ00-0.5*AJ10);
    }

    if(k ==1){
      if(j==1){
        val = -j*(j*AJ00-0.5*AJ10-2.0*pow(alpha,-2.0));

      }else{
        val = -j*(j*AJ00-0.5*AJ10);

      }
    }
  }

  return(val);

}





double C2(int i, int k, int j, double kappa, double beta, double alpha)
{

  double val=0.0;
  if(i==1){
    if(k==0){
      val = d2(i,k,j,kappa,beta,alpha);
    }
    if(k ==1){
      if(j==1){
        val = alpha*(j*AJ10-0.5*AJ20-alpha);
      }else{
        val = alpha*(j*AJ10-0.5*AJ20);
      }
    }
    if(k ==-1){
      if(j==1){
        val = alpha*(-j*AJ10-0.5*AJ20+alpha);
      }else{
        val = alpha*(-j*AJ10-0.5*AJ20);
      }
    }


    if(k ==2){
      val = alpha*(-j*AJ10-0.5*AJ11);
    }

    if(k ==-2){
      if(j==1){
        val = alpha*(j*AJ10-0.5*AJ11-2.0*alpha);
      }else{
        val = alpha*(j*AJ10-0.5*AJ11);
      }
    }
  }else{
/* Outer Planet */
    if(k==0){
      val = d2(i,k,j,kappa,beta,alpha);
    }
    if(k ==2){
      if(j==1){
        val = (-j*AJ01-0.5*AJ02+pow(alpha,-2.0));
      }else{
        val = (-j*AJ01-0.5*AJ02);
      }
    }
    if(k ==-2){
      if(j==1){
        val = (j*AJ01-0.5*AJ02-pow(alpha,-2.0));
      }else{
        val = (j*AJ01-0.5*AJ02);
      }
    }


    if(k ==-1){
      val = (-j*AJ01-0.5*AJ11);
    }

    if(k ==1){
      if(j==1){
        val = (j*AJ01-0.5*AJ11-2.0*pow(alpha,-2.0));

      }else{
        val = (j*AJ01-0.5*AJ11);

      }
    }
  }

  return(val);

}



double principal_value(double theta)
{
  theta -= 2.0*PI*floor(theta/(2.0*PI));
  return(theta);
}
double f_j_k_i(double alpha,int j,int k, int i,double kappa_o_m, double beta_o_m,double *a00,double *a01, double *a10, double *a20, double *a02,double *a11)
{
  double val=0.0;
  double beta,kappa;

  AJ00= a00[j];
  AJ01= a01[j];
  AJ02= a02[j];
  AJ20= a20[j];
  AJ10= a10[j];
  AJ11= a11[j];
  beta = beta_o_m*j;
  kappa = kappa_o_m*j;
  val=u(gam(i,k,j,kappa,beta,alpha),C1(i,k,j,kappa,beta,alpha),C2(i,k,j,kappa,beta,alpha));
  //  printf("%d %d %d %lf %lf %lf %lf\n",i,k,j,val,gam(i,k,j,kappa,beta,alpha),C1(i,k,j,kappa,beta,alpha),C2(i,k,j,kappa,beta,alpha));
  if(i == abs(k)){

    if( k > 0){
      val+= v_plus(xi(i,k,j,kappa,beta,alpha),d1(i,k,j,kappa,beta,alpha),d2(i,k,j,kappa,beta,alpha));
    }
    else{
      val+= v_minus(xi(i,k,j,kappa,beta,alpha),d1(i,k,j,kappa,beta,alpha),d2(i,k,j,kappa,beta,alpha));
    }
  }
  return(val);
}

double u(double gamma, double C1v, double C2v){
  double val=0.0;
  double gsq = gamma*gamma;
  val = ((3.0+gsq)*C1v+2.0*gamma*C2v)/gsq/(1.0-gsq);
  return(val);
}


double v_plus(double z, double d1v, double d2v){
  double val=0.0;
  double zsq = z*z;
  val = (((1.0-zsq)+6.0*z)*d1v+(2.0+zsq)*d2v)/(z*(1.0-zsq)*(z+1.0)*(z+2.0));
  return(val);
}

double v_minus(double z, double d1v, double d2v){
  double val=0.0;
  double zsq = z*z;
  val = ((-(1.0-zsq)+6.0*z)*d1v+(2.0+zsq)*d2v)/(z*(1.0-zsq)*(z-1.0)*(z-2.0));
  return(val);
}
