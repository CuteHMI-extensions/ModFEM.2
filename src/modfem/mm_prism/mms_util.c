/************************************************************************
File mms_util.c - general purpose utilities

Contains routines:
  mmr_chk_list - list manipulation
  mmr_ivector - to allocate space for vectors
  mmr_vec3_prod - to compute vector product of 3D vectors
  mmr_vec3_mxpr - to compute mixed vector product of 3D vectors
  mmr_vec3_length - to compute length of a 3D vector


------------------------------
History:
	  02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

/* interface for the mesh manipulation module */
#include <modfem/mmh_intf.h>

/* mesh manipulation data structure and headers for internal routines */
#include <modfem/mm_prism/mmh_prism.h>


/*---------------------------------------------------------
mmr_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int mmr_chk_list(	/* returns: */
		/* >0 - position on the list */
		/* 0 - not found on the list */
		int Num, 	/* number to be checked */
		int * List, 	/* list of numbers */
		int Ll		/* length of the list */
)
{

	int i, il;

	for (i = 0; i < Ll; i++) {
		if ((il = List[i]) == 0) break;
		/* found on the list on (i+1) position */
		if (Num == il) return (i + 1);
	}
	/* not found on the list */
	return (0);
}

/*---------------------------------------------------------
mmr_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
void mmr_vec3_prod(
		double * vec_a, 	/* in: vector a */
		double * vec_b, 	/* in: vector b */
		double * vec_c	/* out: vector product axb */
)
{

	vec_c[0] = vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
	vec_c[1] = vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2];
	vec_c[2] = vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

	return;
}

/*---------------------------------------------------------
mmr_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
double mmr_vec3_mxpr( /* returns: mixed product [a,b,c] */
		double * vec_a, 	/* in: vector a */
		double * vec_b, 	/* in: vector b */
		double * vec_c	/* in: vector c */
)
{
	double daux;

	daux  = vec_c[0] * (vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1]);
	daux += vec_c[1] * (vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2]);
	daux += vec_c[2] * (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

	return (daux);
}

/*---------------------------------------------------------
mmr_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
double mmr_vec3_length(	/* returns: vector length */
		double * vec	/* in: vector */
)
{
	return (sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
}

/*---------------------------------------------------------
 mmr_ivector - to allocate an integer vector: name[0..ncom-1]:
				  name=mmr_ivector(ncom,error_text)
---------------------------------------------------------*/
int * mmr_ivector( /* returns: pointer to array of integers */
		int ncom, 		/* in: number of components */
		char error_text[]	/* in: text to print in case of error */
)
{
	int * v;

	v = (int *) malloc (ncom * sizeof(int));
	if (!v) {
		printf("Not enough space for allocating vector: %s ! Exiting\n", error_text);
		exit(-1);
	}

	return v;
}

