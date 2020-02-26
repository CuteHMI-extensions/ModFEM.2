#ifndef MMH_VEC3_H
#define MMH_VEC3_H

#include "../Common.h"

#ifdef __cplusplus
extern "C"{
#endif

void mmr_vec3_add(
    const double vec_a[3],
    const double vec_b[3],
    IN double vec_c[3]);

void mmr_vec3_subst(
    const double vec_a[3],
    const double vec_b[3],
    IN double vec_c[3]);

//----------------------------------------------------------
//mmr_vec3_dot - to compute vector dot product of 3D vectors
//----------------------------------------------------------
double mmr_vec3_dot(
    const double vec_a[3],
    const double vec_b[3]);
/**--------------------------------------------------------
mmr_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
void mmr_vec3_prod(
    const double* vec_a, 	/* in: vector a */
    const double* vec_b, 	/* in: vector b */
    double* vec_c	/* out: vector product axb */
);
/**--------------------------------------------------------
mmr_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
double mmr_vec3_mxpr( /* returns: mixed product [a,b,c] */
    const double* vec_a, 	/* in: vector a */
    const double* vec_b, 	/* in: vector b */
    const double* vec_c	/* in: vector c */
) ;
/**--------------------------------------------------------
mmr_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
double mmr_vec3_length(	/* returns: vector length */
    const double* vec	/* in: vector */
						);
/**--------------------------------------------------------
mmr_distance3d - to compute distance between points in 3D
---------------------------------------------------------*/
double mmr_distance3d( // returns distance between p1 and p2
					  const double p1[3], // point 1
					  const double p2[3]); // point 2

/**--------------------------------------------------------
mmr_point_plane_dist - to compute shortest distance between point and plane
---------------------------------------------------------*/
double mmr_point_plane_dist( // returns distance between p1 and p2
					  const double point[3], // point coords
					  const double plane[4]); // plane general equation

#ifdef __cplusplus
}
#endif 

#endif
