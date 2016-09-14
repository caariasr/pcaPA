/*-------------------------------------------------------
  polychoric.cpp
  Function for calculating polychoric correlations

  Carlos A. Arias <carias@icfes.gov.co>
  2012-03-27
---------------------------------------------------------*/

/*  C++ Headers */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cfloat>
#include <ctime>
#include <sys/types.h>
#include <unistd.h>

/* This definition is to prevent R headers to use C headers that might 
 * cause conflict with the ones of C++ */

#ifndef NO_C_HEADERS
#define NO_C_HEADERS
#endif

/* R headers */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/Applic.h>

/*  GSL Headers */

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


using namespace std;


#ifndef Pi
const double Pi = 3.1415926535897932384626433832795028841971693993751058;
#endif
 
/*  Function headers */

void bin(int *,int *,int n, int *, int *);
void copy(int *,int *, int);
void copydoub(double *,double *, int);
void getColumn(int *,int *,int , int);
void max_array(int *, int *, int *);
void colcuts(int *, int, int, int, double *);
void rowcuts(int *, int, int, int, double *);
void unique(int *, int, int *);
double function(double, double *, double *, int, int, int *);
double N(double, double, double);
double nc(double);
double ndf(double);
double Ngen(double, double, double, double, double);
double sgn(double);
double f(double, double, double, double, double);
double Brent_fmin(double, double, double (*)(double, double *, double *, int, int, int *), double, double *, double *, int, int, int *);
void correlation(int *, int *, int *, int, int, int, double *);
void probabilities(int *, int, int, double *);
void simulatedBin(int, int, int, gsl_rng *, double *, int *, int *);
void simulatedMulti(int, int, int, int, double *, gsl_rng *, int *, int *);
void eigenvalues(double *, int, double *);
void vecplus (int *, int);
void sum (int *, int, int *);
void tabulate(int *, int *, int *, int *);


/* Main Function */

extern "C" SEXP Cpolychoric(SEXP Rdata, SEXP isBinary, SEXP nReplicates, SEXP nearcor, SEXP rho) {
  
  /* Define environment and seed for the simulations */
  const gsl_rng_type * T; 
  gsl_rng *r_global;

  long seed = time(NULL) * getpid();

  gsl_rng_env_setup ();
  
  T = gsl_rng_default; 

  r_global = gsl_rng_alloc (T);

  gsl_rng_set(r_global, seed);

  
  /* Protect and pointer assignment section */
  PROTECT(Rdata = AS_INTEGER(Rdata));
  int *pRdata = INTEGER_POINTER(Rdata);

  PROTECT(isBinary = AS_INTEGER(isBinary));
  int visBinary = INTEGER_VALUE(isBinary);

  PROTECT(nReplicates = AS_INTEGER(nReplicates));
  int vnReplicates = INTEGER_VALUE(nReplicates);

  /*  Extract number of columns and rows from data */

  SEXP dim;
  PROTECT(dim = NEW_INTEGER(2));
  dim = getAttrib(Rdata, R_DimSymbol);
  int *pdim = INTEGER_POINTER(dim);

  int cols = pdim[1];
  int tot = cols*cols;
  int n; 
  int *pn;
  pn=&n;

  *pn = pdim[0];

  /*  OBSERVED STEP */

  /* Define vec1 and vec2, that will contain columns i and j of the observed data in each iteration */

  int *vec1 =  (int*) calloc (n, sizeof(int));
  int *pvec1=vec1;
  int *vec2 =  (int*) calloc (n, sizeof(int));
  int *pvec2=vec2;

  /* Define correla, array that will contain the polychoric correlation matrix 
   * of the observed data */

  double *correla = (double*) calloc(cols*cols,sizeof(double));
  double *pcorrela = correla;

  /*  Define eigenv, array that will contain the eigenvalues of the polychoric 
   *  correlation matrix */

  double *eigenv = (double*) calloc(cols,sizeof(double));
  double *peigenv = eigenv;

  /*  Calculate the polychoric correlations for the observed data, stored in
   *  correla */

  correlation(pRdata, pvec1, pvec2, cols, n, tot, pcorrela);
  
  /*  Define R objects for input and calculation of the nearcor R function */

  SEXP Rcorrela, R_fcall, R_fcalli;
  PROTECT (Rcorrela = allocMatrix(REALSXP, cols, cols));
  double *pRcorrela = NUMERIC_POINTER(Rcorrela);

  int length = cols * cols;

  /*  Copy the polychoric correlation matrix into the R matrix Rcorrela */

  copydoub(pRcorrela,pcorrela,length);
  
  PROTECT(R_fcall = lang2(nearcor, Rcorrela));

  SEXP correlaAprox, correlaAproxi;
  
  /* Approximate the polychoric correlation matrix to the nearest positive 
   * definite */
  PROTECT(correlaAprox = VECTOR_ELT(eval(R_fcall,rho),0));
  double *pcorrelaAprox = NUMERIC_POINTER(correlaAprox);

  /* Calculate the eigenvalues of the approximated matrix */

  eigenvalues(pcorrelaAprox, cols, peigenv);

  /*  Copy the eigenvalues to an R object */

  SEXP REigen;
  
  PROTECT(REigen = NEW_NUMERIC(cols));
  double *pREigen = NUMERIC_POINTER(REigen);

  copydoub(pREigen,peigenv,cols);

  
  /*  SIMULATION STEP */
  /*  Define vec, that will contain the ith column in the jth iteration of 
   *  the observed data */ 

  int *vec = (int*) calloc(n,sizeof(int));
  int *pvec = vec;
  
  int simtot = vnReplicates * cols;
  
  /* Define EigenMat, the matrix containing the Eigenvalues of each iteration */

  double *EigenMat = (double*) calloc(simtot,sizeof(double));
  double *pEigenMat = EigenMat;
  
  /*  Define data, the matrix that will contain the simulated data in each 
   *  iteration */

  int *data = (int*) calloc(cols * n,sizeof(int));
  int *pdata = data;

  /*  Main loop for the simulation process */

  for (int j = 0; j < vnReplicates; j++) {
    for (int i = 0; i < cols; i++) {
      /*  Extract the ith column of the observed data */
      getColumn(pRdata,pvec,n,i);
      /*  Obtain unique values of the ith column of the observed data */
      int uni=1;
      int *puni;
      puni=&uni;
      unique(pvec,n,puni);
      /*  This function increases the value of all categories in one, to avoid
       *  conflict with function tabulate */
      vecplus(pvec, n);
      /*  Calculate the average of each category, taken as probabilities of 
       *  occurrence */
      double *probs = (double*) calloc(uni,sizeof(double));
      double *pprobs = probs;
      probabilities(pvec, n, uni, pprobs);

      /*  If the data are dichotomous, sample from a Binomial distribution 
       *  using the average of success as parameter */

      if (visBinary == 1) {
        simulatedBin(cols, n, i, r_global, pprobs, pRdata, pdata);
      }  else { /*  If the data are ordered, sample from a Multinomial 
                    distribution using the category averages as parameters */
        simulatedMulti(cols, n, i, uni, pprobs, r_global, pRdata, pdata); 
      }
      free(probs);
    }

    /*  Calculate the polychoric correlation of the simulated data in the jth 
     *  iteration */
    correlation(pdata, pvec1, pvec2, cols, n, tot, pcorrela);

    /*  Approximate the polychoric correlation matrix to the nearest positive
     *  definite matrix */        
    copydoub(pRcorrela,pcorrela,length);
    PROTECT(R_fcalli = lang2(nearcor, Rcorrela));
    PROTECT(correlaAproxi = VECTOR_ELT(eval(R_fcalli,rho),0));
    double *pcorrelaAproxi = NUMERIC_POINTER(correlaAproxi);
    
    /*  Calculate the eigenvalues of the approximated matrix */
    eigenvalues(pcorrelaAproxi, cols, peigenv);
    UNPROTECT(2);
    /* Copy the eigenvalues tho the matrix EigenMat */
    for (int k = 0; k < cols; k++) {
      pEigenMat[k + j*cols]=peigenv[k];
    }
  }
  
  /*  Define REigenMat, R matrix that will contain the eigenvalues for all 
   *  the iterations */
  SEXP REigenMat;

  PROTECT(REigenMat = NEW_NUMERIC(simtot));
  double *pREigenMat = NUMERIC_POINTER(REigenMat);

  /*  Copy the eigenvalues from EigenMat to REigenMat */

  copydoub(pREigenMat,pEigenMat,simtot);
  
 
  /*  Define list, R list that will contain the results to be passed to R,
   *  the observed and simulated eigenvalues */ 
  SEXP list; 

  PROTECT(list = allocVector(VECSXP,2));

  
  
  SET_VECTOR_ELT(list,0,REigen);
  SET_VECTOR_ELT(list,1,REigenMat);
 

  /*  Free section */
  free(data);
  free(vec1);
  free(vec2);
  free(correla);
  free(eigenv);
  free(vec);
  free(EigenMat);

  /*  Unprotect section */
  UNPROTECT(10);


  /*  Return the output list */
  return list;


}
/*  This functions counts the total frequencies of the cross-tabulation
 *  between vec1 and vec2 */

void sum (int *pfreqs, int cross, int *ptotalfreqs) {
  for (int i = 0; i < cross; i++) {
    *ptotalfreqs = *ptotalfreqs + pfreqs[i];
  } 
} 


/*  This function increments in one the value of all the categories found in 
 *  vec, this is to avoid conflicts with tabulate, which only
 *  counts categories different to 0 */

void vecplus (int *pvec, int n) {
  for (int i = 0; i < n; i++) {
    pvec[i]=pvec[i]+1;
  }
}

/*  This function calculates the eigenvalues of a given matrix (in this case 
 *  correla) */

void eigenvalues(double *pcorrela, int cols, double *peigenv) {
  
  int length = cols * cols;
  double *correlacopy = (double*) calloc(length,sizeof(double));
  double *pcorrelacopy = correlacopy; 

  copydoub(pcorrelacopy,pcorrela,length);

  gsl_matrix_view X = gsl_matrix_view_array(pcorrelacopy,cols,cols);

  gsl_vector *eval = gsl_vector_alloc (cols);

  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(cols);

  gsl_eigen_symm(&X.matrix, eval, w);

  gsl_eigen_symm_free (w);

  gsl_sort_vector(eval);

  free(correlacopy);

  for (int k = 0; k < cols; k++) {
    peigenv[k] = gsl_vector_get(eval,k);
  }

  gsl_vector_free(eval);
}

/*  This function calculates the average of each category found in vec */

void probabilities(int *pvec, int n, int uni, double *pprobs) {
  int *pn;
  pn=&n;
  int *freqs = (int*) calloc(uni,sizeof(int));
  int *pfreqs = freqs;
  int max;
  int *pmax;
  pmax=&max;
  *pmax=0; 
  max_array(pvec,pn,pmax);
  tabulate(pvec,pn,pmax,pfreqs);
  for (int i = 0; i < uni; i++) {
    double temp = pfreqs[i];
    pprobs[i] = temp / n;
  }
  free(freqs);
}

/*  This function simulates data from a binomial distribution taking into 
 *  account missing values in the observed data 
 *  so that the simulated data contains this same missing values */

void simulatedBin(int cols, int n,int j, gsl_rng *r_global, double *pprobs, 
    int *pRdata, int *pdata) {
  for (int i = 0; i < n; i++) {
    if (pRdata[i+j*n] == NA_INTEGER)
      pdata[i+j*n] = NA_INTEGER;
    else {
      double prob = pprobs[1];
      unsigned int npar = 1;
      unsigned int temp = gsl_ran_binomial(r_global, prob, npar);
      pdata[i+j*n] = temp;
    }
  }
} 

/*  This function simulates data from a multinomial distribution taking into 
 *  account missing values in the observed data 
 *  so that the simulated data contains this same missing values */

void simulatedMulti(int cols, int n, int j, int uni, double *pprobs,
    gsl_rng *r_global, int *pRdata,  int *pdata) { 
 
  unsigned int *test = (unsigned int *) calloc(uni,sizeof(unsigned int));
  unsigned int *ptest = test;

  int sum_n=0;
  for (int k = 0; k < uni; k++) {
    sum_n += k;
  }
  for (int i = 0; i < n; i++) {
    if (pRdata[i+j*n] == NA_INTEGER)
      pdata[i+j*n] = NA_INTEGER;
    else {
      gsl_ran_multinomial(r_global, uni, sum_n, pprobs, test);
      for ( int l = 0; l < uni; l++) {
        if (ptest[l] == 1) pdata[i+j*n] = l;
      }
    }
  }
}

/*  This function calculates the polychoric correlation between vec1 and vec2 
 *  using two-step estimation  */
void correlation(int *pRdata, int *pvec1, int *pvec2, int cols, 
    int n, int tot, double *pcorrela) {

  int *pn;
  pn=&n;
  /*  Main loop */
  for (int i=0; i < cols - 1; i++) {
    /*  Extract the ith column of the observed data */
    getColumn(pRdata,pvec1,n,i);
    /*  Seconday loop  */
    for(int j=i; j < cols; j++) {
      /* Assign 1 to correlations in the diagonal */
      if (j == i) {
        pcorrela[i+j*cols]=1.0;
        continue;
      }
      else {
        /*  Extract the jth column of the observed data */
        getColumn(pRdata,pvec2,n,j);
        /*  Define constants */
        int uni1=1;
        int uni2=1;
        int *puni2;
        puni2=&uni2;
        int *puni1;
        puni1=&uni1;

        int max;
        int *pmax;
        pmax=&max;
        *pmax=0; 
        
        unique(pvec2,n,puni2);

        unique(pvec1,n,puni1);

        int cross = uni1 * uni2;
        
        /*  Define combi, vector that will containing all the combinations 
         *  of categories in vec1 and vec2 */ 

        int *combi = (int*) calloc (n,sizeof(int));
        int *pcombi = combi;

      
        bin(pvec1,pvec2,n,puni1,pcombi);
        
        /*  Define freqs, array that will contain the frequencies of the
         *  categories found in combi, this will conform the cross-tab of
         *  vec1 and vec2 */ 
 
        int *freqs = (int*) calloc(cross,sizeof(int));
        int *pfreqs = freqs;

        /*  Find the maximum value in combi array */
        max_array(pcombi,pn,pmax);

        /*  Define rc and cc, row cuts and col cuts respectively, see 
         *  function's description */
        double *rc = (double*) calloc(uni1+1,sizeof(double));
        double *prc = rc;
        double *cc = (double*) calloc(uni2+1,sizeof(double));
        double *pcc = cc;
     
        /*  Calculate the frequencies of the categories found in combi */
        tabulate(pcombi,pn,pmax,pfreqs);

        free(combi);

        /* Calculate the total sum of freqs */

        int totalfreqs;
        int *ptotalfreqs;
        ptotalfreqs = & totalfreqs;
        *ptotalfreqs = 0;
        sum(pfreqs, cross, ptotalfreqs);

        
        /*  Calculate the cumulative probabilities for rows and columns */
        

        colcuts(pfreqs,uni1,uni2,totalfreqs,pcc);
        rowcuts(pfreqs,uni1,uni2,totalfreqs,prc);

        double ax  = -1.0;
        double bx  =  1.0;
        double tol =  0.0001220703; 

        pcorrela[i+j*cols]=Brent_fmin(ax,bx,function,tol,prc,pcc,uni1,uni2,pfreqs);
        pcorrela[j+i*cols]=pcorrela[i+j*cols];
        free(freqs);
      }
    }
  }
  pcorrela[tot-1]=1.0;
}

/*  This creates combi as a an array whose categories are all the possible
 *  combinations of the categories in vec1 with the ones of vec2 */
void bin(int *pvec1,int *pvec2,int n, int *puni1, int *pcombi) {
   
  for (int i=0; i < n; i++) {
    *(pcombi+i)=*(pvec1+i);
  }

  for (int k=0; k < n; k++) {
    if ((pvec1[k] != R_NaInt) & (pvec2[k] != R_NaInt))
      pcombi[k]=pvec1[k]+((*puni1)*pvec2[k])+1;
    else pcombi[k] = R_NaInt;
  } 
}

/*  Generic function to copy from one integer array to another */

void copy(int *to,int *from, int len){
  int *toPtr = to;
  int *fromPtr = from;
  for(int i=0;i<len;i++){
    *(toPtr++) = *(fromPtr++);
 }
}

/*  Generic function to copy from one double array to another */

void copydoub(double * to,double *from, int len){
  double *toPtr = to;
  double *fromPtr = from;
  for(int i=0;i<len;i++){
    *(toPtr++) = *(fromPtr++);
 }
}

/*  Generic function to extract the jth column from an array of given
 *  dimensions */

void getColumn(int *D,int *pCol,int m, int j){
int *pc= pCol;
for (int l=0;l<m;l++){
    *(pc++)=  *(D+j*m+l);
}
}

/* Function to fin the maximum value of an array  */
void max_array(int *a, int *pn, int *pmax) {
  for (int i=0; i < *pn; i++) {
    if (a[i] != R_NaInt) {
      if (a[i]>*pmax) {
        *pmax=a[i];
      }
    }
  }
}



/*  This function calculates the normal quantiles taking as parameters the
 *  cumulative sums of the column sums */
void colcuts(int *pfreqs, int uni1, int uni2, int n, double *pcc) {
  int sum=0;
  *(pcc+0)=-1000.0;
  for (int i=1; i < uni2; i++) {
    for (int j=1; j < uni1+1; j++) {
      sum += *(pfreqs+(j-1)+(i-1)*uni1);
    }
    *(pcc+i)=sum;
    sum=0;
  }
  for (int k=1; k < uni2; k++) {
    sum+=*(pcc+k);
    *(pcc+k) = sum;
    *(pcc+k) = *(pcc+k)/n;
    *(pcc+k) = qnorm(*(pcc+k),0,1,1,0);
  }
  *(pcc + uni2) = 1000.0;
}

/*  This function calculates the normal quantiles taking as parameters the
 *  cumulative sums of the row sums */

void rowcuts(int *pfreqs, int uni1, int uni2, int n, double *prc) {   
  int sum=0;
   *(prc+0)=-1000.0;
  for (int i=1; i < uni1; i++) {
    for (int j=1; j < uni2+1; j++) {
      sum += *(pfreqs+(i-1)+(j-1)*uni1);
    }
    *(prc+i) = sum;
    sum=0;
  }
  for (int k=1; k < uni1; k++) {
    sum+=*(prc+k);
    *(prc+k) = sum;
    *(prc+k) = *(prc+k)/n;
    *(prc+k) = qnorm(*(prc+k),0,1,1,0);
  }
  *(prc + uni1) = 1000.0;
}

/*  This functions counts the number of unique values in vec */

void unique(int *pvec, int n, int *puni) {
  int *vecdup = (int*) calloc(n,sizeof(int));
  int *pvecdup = vecdup; 
  for (int i=0; i < n; i++) {
    *(pvecdup+i)=*(pvec+i);
  }
  R_isort(pvecdup,n);
  for (int j=0; j < n - 1; j++) {
    if ((pvecdup[j] != R_NaInt) & (pvecdup[j+1] != R_NaInt) ) {
      if (pvecdup[j]==pvecdup[j+1])
        continue;
      else
        *puni=*puni+1;
    }
  }
  free(vecdup);
}

/*  Function to be minimized */

double function(double rho, double *prc, double *pcc, int uni1, int uni2, int *pfreqs) {
  int cont=0;
  double *P = (double*) calloc(uni1*uni2,sizeof(double));
  double *pP = P;
  for (int j = 0; j < uni2; j++) {
   for (int i = 0; i < uni1; i++) {    
      pP[cont] = pfreqs[cont]*log(Ngen(prc[i],prc[i+1],pcc[j],pcc[j+1],rho));
      cont++;
    }
  }
  double res = 0.0;
  for (int i = 0; i < uni1*uni2; i++) {
    res += pP[i];
  }
  free(P);
  res = 0 -res;
  return res;
}

/*  Function to obtain values from the bivariate normal distribution 
 *  function. Using numerical approximation as found in
 *  http://finance.bi.no/~bernt/gcc_prog/recipes/recipes/node23.html#SECTION002340000000000000000
 *  */


double N(double a, double b, double rho){
  if((a<=0.0) && (b<= 0.0) && (rho<=0.0)) {
    double aprime = a/sqrt(2.0*(1.0-rho*rho));
    double bprime = b/sqrt(2.0*(1.0-rho*rho));
    double A[4] = {0.3253030, 0.4211071, 0.1334425, 0.006374323};
    double B[4] = {0.1337764, 0.6243247, 1.3425378, 2.2626645 };
    double sum = 0;
    for (int i=0;i<4;i++){
       for (int j=0;j<4;j++){
         sum += A[i]*A[j]*f(B[i],B[j],aprime,bprime,rho);
       };
    };
    sum = sum * (sqrt(1.0-rho*rho)/Pi);
    return sum;
  }
  else if (a*b*rho <= 0.0) {
    if ( ( a<=0.0 ) && ( b >=0.0 ) && (rho >=0.0)) {
      return nc(a) - N(a, -b, -rho);
    }
    else if ((a >= 0.0) && (b <= 0.0) && (rho >=0.0)) {
      return nc(b) - N(-a, b, -rho);
    }
    else if ((a >= 0.0) && (b >= 0.0) && (rho <= 0.0)) {
      return nc(a) + nc(b) - 1.0 + N(-a, -b, rho);
    };
  }
  else if (a * b * rho >= 0.0) {
    double denum = sqrt(a*a - 2*rho*a*b + b*b);
    double rho1 = ((rho * a - b) * sgn(a))/denum;
    double rho2 = ((rho * b - a) * sgn(b))/denum;
    double delta = (1.0 - sgn(a)*sgn(b))/4.0;
    return N(a,0.0,rho1) + N(b, 0.0, rho2) - delta;
  };
  return -99.9;
};  

double nc(double x)
{
  double result;
  if (x<-7.)
    result = ndf(x)/sqrt(1.+x*x);
  else if (x>7.)
    result = 1. - nc(-x);
  else
    {
      result = 0.2316419;
      static double a[5] = {0.31938153,-0.356563782,1.781477937,-1.821255978,1.330274429};
      result=1./(1+result*fabs(x));
      result=1-ndf(x)*(result*(a[0]+result*(a[1]+result*(a[2]+result*(a[3]+result*a[4]))))
          );
      if (x<=0.) result=1.-result;
      }
  return result;
}

double ndf(double t)
{
  return 0.398942280401433*exp(-t*t/2);
}



double f(double x, double y, double aprime, double bprime, double rho){
  double r = aprime*(2*x-aprime) + bprime*(2*y-bprime) + 2*rho*(x-aprime)*(y-bprime);
  return exp(r);
};


double sgn(double x){
  if (x>=0.0) return 1.0;
  return -1.0;
};

/*  Function to calculate bivariate normal probabilities for two given 
 *  intervals and a constant correlation parameter rho  */
double Ngen(double a, double b, double c, double d, double rho) {
 if( a== -1000.0 && c == -1000.0) 
   return N(b,d,rho);
 if(a == -1000.0 && d == 1000.0)  
   return nc(b) - N(b,c,rho);
 if(c == -1000.0 && b == 1000.0) 
   return nc(d) - N(a,d,rho);
 if (a == -1000.0 && d != 1000.0)
   return N(b,d,rho) - N(b,c,rho);
 if (c == -1000.0 && b != 1000.0)
   return N(b,d,rho) - N(a,d,rho);
 if (a != -1000.0 && d == 1000.0)
   return nc(b) - nc(a) - N(b,c,rho) + N(a,c,rho);
 if (c != -1000.0 && b == 1000.0)
   return nc(d) - nc(c) - N(a,d,rho) + N(a,c,rho);
 if (b == 1000 && d == 1000)
   return N(-a,-c,rho);
 else return N(b,d,rho) + N(a,c,rho) - N(a,d,rho) - N(b,c,rho);    
}


/*  Brent_fmin taken from the optimize R function and slightly modified to fit this program */
double Brent_fmin(double ax, double bx, 
    double (*f)(double, double *, double *, int, int,  int *), 
    double tol, double *prc, double *pcc, int uni1, int uni2, int *pfreqs)
{

    const double c = (3. - sqrt(5.)) * .5;


    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;


    eps = DBL_EPSILON;
    tol1 = eps + 1.;
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;
    e = 0.;
    fx = (*f)(x,prc,pcc,uni1,uni2,pfreqs);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;


    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;


	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) {

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { 
	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { 

	    d = p / q;
	    u = x + d;

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u,prc,pcc,uni1,uni2,pfreqs);

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    return x;
}

/* Function to tabulate data taken from the former R_tabulate C++
 * function in the R source */

void tabulate(int *pvec, int *pn, int *pmax, int *pfreqs)
{
    int i;
    if(*pn < 1) return;
    for(i = 0 ; i < *pn ; i++)
	if(pvec[i] != R_NaInt && pvec[i] > 0 && pvec[i] <= *pmax)
	    pfreqs[pvec[i] - 1]++;
}

