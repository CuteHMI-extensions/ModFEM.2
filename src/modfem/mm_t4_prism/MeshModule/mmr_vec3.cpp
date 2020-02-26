
#include <cmath>
#include "mmh_vec3.h"


void mmr_vec3_add(
    const double vec_a[3],
    const double vec_b[3],
    IN double vec_c[3])
{
    vec_c[0]=vec_a[0]+vec_b[0];
    vec_c[1]=vec_a[1]+vec_b[1];
    vec_c[2]=vec_a[2]+vec_b[2];
}


void mmr_vec3_subst(
    const double vec_a[3],
    const double vec_b[3],
    IN double vec_c[3])
{
    vec_c[0]=vec_a[0]-vec_b[0];
    vec_c[1]=vec_a[1]-vec_b[1];
    vec_c[2]=vec_a[2]-vec_b[2];
}


//----------------------------------------------------------
//mmr_vec3_dot - to compute vector dot product of 3D vectors
//----------------------------------------------------------
double mmr_vec3_dot(
    const double vec_a[3],
    const double vec_b[3])
{
    return vec_a[0]*vec_b[0]+vec_a[1]*vec_b[1]+vec_a[2]*vec_b[2];
}

/*---------------------------------------------------------
mmr_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
void mmr_vec3_prod(
    const double* vec_a, 	/* in: vector a */
    const double* vec_b, 	/* in: vector b */
    double* vec_c	/* out: vector product axb */
) {

    vec_c[0]=vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1];
    vec_c[1]=vec_a[2]*vec_b[0]-vec_a[0]*vec_b[2];
    vec_c[2]=vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0];

}

/*---------------------------------------------------------
mmr_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
double mmr_vec3_mxpr( /* returns: mixed product [a,b,c] */
    const double* vec_a, 	/* in: vector a */
    const double* vec_b, 	/* in: vector b */
    const double* vec_c	/* in: vector c */
) {

    double daux  = vec_c[0]*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1]);
    daux += vec_c[1]*(vec_a[2]*vec_b[0]-vec_a[0]*vec_b[2]);
    daux += vec_c[2]*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);

    return(daux);
}

/*---------------------------------------------------------
mmr_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
double mmr_vec3_length(	/* returns: vector length */
    const double* vec	/* in: vector */
) {
    return(sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]));
}

/*---------------------------------------------------------
mmr_distance3d - to compute distance between points in 3D
---------------------------------------------------------*/
double mmr_distance3d( // returns distance between p1 and p2
					  const double p1[3], // point 1
					  const double p2[3]) // point 2
{
  const double tmpVec[3]={p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
  return mmr_vec3_length(tmpVec);
}


/*---------------------------------------------------------
mmr_point_plane_dist - to compute shortest distance between point and plane
---------------------------------------------------------*/
double mmr_point_plane_dist( // returns distance between p1 and p2
					  const double x[3], // point coords
					  const double A[4]) // plane general equation
{
  assert(x!= NULL);
  assert(A != NULL);
  double dist
	= std::abs(A[0]*x[0] + A[1]*x[1] + A[2]*x[2] + A[3])
	/ sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
  
  assert(dist >= 0.0);
  return dist;
}


