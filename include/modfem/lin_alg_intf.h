
#ifndef _lin_alg_intf_
#define _lin_alg_intf_

#ifdef __cplusplus
extern "C" {
#endif 

/** @defgroup LIN_ALG Linear Algebra
 *
 *  @{
 */

/** BLAS and LAPACK names with or without underscore */
#ifdef WITHOUT_
#define daxpy_ daxpy
#define dcopy_ dcopy
#define ddot_ ddot
#define dnrm2_ dnrm2
#define drot_ drot
#define drotg_ drotg
#define dscal_ dscal
#define dgemv_ dgemv
#define dtrsv_ dtrsv
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#define dgemm_ dgemm
#endif

/** BLAS 1 */
void   daxpy_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
void   dcopy_(int *n,double *x,int *incx,double *y,int *incy);
double ddot_(int *n,double *x,int *incx,double *y,int *incy);
double dnrm2_(int *n,double *x,int *incx);
void   drot_(int *n,double *x,int *incx,double *y,int *incy,double *c,double *s);
void   drotg_(double *a,double *b,double *c,double *s);
void   dscal_(int *n,double *a,double *x,int *incx);

/** BLAS 2 */
void dgemv_(char *trans,int *m,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
void dtrsv_(char *uplo,char *trans,char *diag,int *n,double *a,int *lda,double *x,int *incx);
void dgemm_(char *transa, char *transb, int *m, int *
	n, int *k, double *alpha, double *a, int *lda,
	double *b, int *ldb, double *beta, double *c__,
	int *ldc);

/** LAPACK */
void    dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
int	dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
void    dgetrs_(char *trans,int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info);

/** @} */ // end of group


#ifdef __cplusplus
}
#endif 

#endif
