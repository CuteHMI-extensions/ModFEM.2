/*************************************************************
File contains procedures:
lar_util_dvector - to allocate a double vector: name[0..ncom-1]:
lar_util_ivector - to allocate an integer vector: name[0..ncom-1]:
lar_util_imatrix - to allocate an integer matrix name[0..nrow-1][0..ncol-1]: 
                  name=imatrix(nrow,ncol,error_text) 
lar_util_dmatrix - to allocate a double matrix name[0..nrow-1][0..ncol-1]: 
                  name=imatrix(nrow,ncol,error_text) 
lar_util_chk_list - to check whether a number is on the list
lar_util_put_list - to put Num on the list List with length Ll (obsolete)
lar_util_d_zero - to zero a double vector
lar_util_i_zero - to zero an integer vector
lar_util_sort - to heap-sort an array
lar_util_dgetrf - quasi-LU decomposition of a matrix
lar_util_dgetrs - to perform forward reduction and back substitution
    of the RHS vector for solving a system of linear equations

Required routines:
	BLAS: dnrm2, daxpy, dcopy, dgemv, dscal, ddot, 
              dgetrf, dgetrs, dgetri, dgemm

History:
        05.2001 - Krzysztof Banas, initial version

*************************************************************/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include "./lah_block.h"

#define SMALL 1.e-15 /* approximation of round-off error */

/*---------------------------------------------------------
lar_util_dvector - to allocate a double vector: name[0..ncom-1]:
                  name=lar_util_dvector(ncom,error_text) 
---------------------------------------------------------*/
double *lar_util_dvector( /* return: pointer to allocated vector */
	int ncom,  	/* in: number of components */
	char error_text[]/* in: error text to be printed */
	)
{

  double *v;
  
  v = (double *) malloc (ncom*sizeof(double));
  if(!v){
    printf("Not enough space for allocating vector: %s ! Exiting\n", error_text);
    exit(1);
  }
  return v;
} 


/*---------------------------------------------------------
lar_util_ivector - to allocate an integer vector: name[0..ncom-1]:
                  name=lar_util_ivector(ncom,error_text) 
---------------------------------------------------------*/
int *lar_util_ivector(    /* return: pointer to allocated vector */
	int ncom, 	/* in: number of components */
	char error_text[]/* in: error text to be printed */
	)
{

  int *v;
  
  v = (int *) malloc (ncom*sizeof(int));
  if(!v){
    printf("Not enough space for allocating vector: %s ! Exiting\n", error_text);
    exit(1);
  }
  return v;
} 

/*---------------------------------------------------------
lar_util_imatrix - to allocate an integer matrix name[0..nrow-1][0..ncol-1]: 
                  name=imatrix(nrow,ncol,error_text) 
---------------------------------------------------------*/
int **lar_util_imatrix( /* returns: pointer to array of pointers to integers */
	int Nrow, 	/* in: number of rows */
	int Ncol, 	/* in: number of columns */
	char Error_text[]/* in: text to print in case of error */
	)
{

  int i;
  int **m;
  
  m = (int **) malloc (Nrow*sizeof(int *));
  if(!m){
    printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
    exit(1);
  }
  for(i=0;i<Nrow;i++){
    m[i] = (int *) malloc (Ncol*sizeof(int));
    if(!m[i]){
      printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
      exit(1);
    }
  }
  return m;
} 

/*---------------------------------------------------------
lar_util_dmatrix - to allocate a double matrix name[0..nrow-1][0..ncol-1]: 
                  name=imatrix(nrow,ncol,error_text) 
---------------------------------------------------------*/
double **lar_util_dmatrix( /* returns: pointer to array of pointers to doubles */
	int Nrow, 	/* in: number of rows */
	int Ncol, 	/* in: number of columns */
	char Error_text[]/* in: text to print in case of error */
	)
{

  int i;
  double **m;
  
  m = (double **) malloc (Nrow*sizeof(double *));
  if(!m){
    printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
    exit(1);
  }
  for(i=0;i<Nrow;i++){
    m[i] = (double *) malloc (Ncol*sizeof(double));
    if(!m[i]){
      printf("Not enough space for allocating array: %s ! Exiting\n", Error_text);
      exit(1);
    }
  }
  return m;
} 

/*---------------------------------------------------------
lar_util_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int lar_util_chk_list(	/* returns: */
			/* >0 - position on the list */
            		/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	)
{

  int i, il;
  
  for(i=0;i<Ll;i++){
    if((il=List[i])==0) break;
    /* found on the list on (i+1) position */
    if(Num==il) return(i+1);
  }
  /* not found on the list */
  return(0);
}


/*---------------------------------------------------------
lar_util_d_zero - to zero a double vector
---------------------------------------------------------*/
void lar_util_d_zero(double *Vec, int Num)
{

  int i;
  
  for(i=0;i<Num;i++){
    Vec[i]=0;
  }
  
}

/*---------------------------------------------------------
lar_util_i_zero - to zero an integer vector
---------------------------------------------------------*/
void lar_util_i_zero(int *Vec, int Num)
{

  int i;
  
  for(i=0;i<Num;i++){
    Vec[i]=0;
  }
  
}


/*---------------------------------------------------------
lar_util_sort - to heap-sort an array (code taken from fortran...)
---------------------------------------------------------*/
void lar_util_sort(
   int    *Ind_array,    /* in/out: index array for sorting */
   double *Val_array     /* in: array of values used for sorting */
   )
{

  int i,j,l,ir,index;
  double q;

  l = Ind_array[0]/2+1;
  ir =  Ind_array[0];
  
  s20: {}

  if(l>1){
    l--;
    index=Ind_array[l];
    q=Val_array[index];
  }
  else{
    index=Ind_array[ir];
    q=Val_array[index];
    Ind_array[ir]=Ind_array[1];
    ir--;
    if(ir==1){
      Ind_array[1]=index;
      return;
    }
  }

  i=l;
  j=2*l;

  s30: {}

  if(j>ir) goto s40;
  if((j<ir)&&(Val_array[Ind_array[j]]<Val_array[Ind_array[j+1]])) j++;
  if(q<Val_array[Ind_array[j]]){
    Ind_array[i]=Ind_array[j];
    i=j;
    j*=2;;
  }
  else{
    j=ir+1;
  }

  goto s30;

  s40: {}

  Ind_array[i]=index;

  goto s20;

}

 
/************************************************************
lar_util_dgetrf - quasi-LU decomposition of a matrix
*************************************************************/
void lar_util_dgetrf(double* a, int m, int* ips)
/*
in:
	a - matrix to decompose 
	m - number of rows
out:
	a - decomposed matrix
	ips - partial pivoting information storage
*/
{

  int i,j,k,mm1,kp1,ipiv,kpiv,ipivot;
  double temp,pivot,big;

  for(i=0;i<m;i++){ips[i]=i;}

  mm1=m-1;
  for(k=0;k<mm1;k++){
    big=0.0;
    for(i=k;i<m;i++){
      ipiv=ips[i];
      if(big<=fabs(a[ipiv+m*k])){
	big=fabs(a[ipiv+m*k]);
	ipivot=i;
      }
    }
    if(big<SMALL*SMALL){
      printf("Zero pivot 1 in LU decomposition (%.15lf)\n",big);
      exit(1);
    }
    kpiv=ips[ipivot];
    if(ipivot!=k){
      ips[ipivot]=ips[k];
      ips[k]=kpiv;
    }
    pivot=a[kpiv+m*k];
    kp1=k+1;
    for(i=kp1;i<m;i++){
      ipiv=ips[i];
      temp=a[ipiv+m*k]/pivot;
      a[ipiv+m*k]=temp;
      for(j=kp1;j<m;j++) a[ipiv+m*j] -= temp*a[kpiv+m*j];
    }
  }
  if(fabs(a[ips[mm1]+m*mm1])<SMALL*SMALL){
    printf("Zero pivot 2 in LU decomposition (%.15lf)\n",
	   a[ips[mm1]+m*mm1]);
    exit(1);
  }
}

/************************************************************
lar_util_dgetrs - to perform forward reduction and back substitution
    of the RHS vector for solving a system of linear equations
************************************************************/
void lar_util_dgetrs(double* a, int m, double* b, double* x, int* ips)
/*
in:
	a - LU decomposed matrix
	m - number of rows
	b - right hand side vector
out:
	x - solution of a system of linear equations
	ips - partial pivoting information storage
*/

{
  int i,j,ipiv;
  double sum;
  
  x[0]=b[ips[0]];
  for(i=1;i<m;i++){
    ipiv=ips[i];
    sum=b[ipiv];
    for(j=0;j<i;j++){ sum -= a[ipiv+m*j]*x[j]; }
    x[i]=sum;
  }
  x[m-1]=x[m-1]/a[ipiv+m*(m-1)];
  for(i=m-2;i>=0;i--){
    ipiv=ips[i];
    sum=x[i];
    for(j=i+1;j<m;j++){ sum -= a[ipiv+m*j]*x[j]; }
    x[i]= sum / a[ipiv+m*i];
  }
}


/*---------------------------------------------------------
lar_util_skip_rest_of_line - to allow for comments in input files
---------------------------------------------------------*/
void lar_util_skip_rest_of_line(
			  FILE *Fp  /* in: input file */
			  )
{
  while(fgetc(Fp)!='\n'){}
  return;
}

