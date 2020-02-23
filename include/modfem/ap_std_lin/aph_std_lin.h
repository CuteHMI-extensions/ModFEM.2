/************************************************************************
File aph_std_lin.h - internal information for standard linear
			 approximation module for tetrahedral and prismatic elements

Contains:
  - constants
  - data types
  - global variables (for the whole module)
  - headers for internal functions:
  apr_select_field - to select the proper field

------------------------------
History:
		02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#ifndef _aph_std_lin_
#define _aph_std_lin_


#ifdef __cplusplus
extern "C" {
#endif

/*** CONSTANTS ***/
// grab defines from aph interface
// KM 06.2011
#include <modfem/aph_intf.h>


/** \addtogroup APM_STD_LIN Standard Linear Approximation
 *  \ingroup APR
 *  @{
 */




#define APC_MAX_NUM_FIELD   10   /* maximal number of fields */
#define APC_MAXEQ 30 /* maximal number of solution vector components */
#define APC_TRUE 1
#define APC_FALSE 0

/*** DATA TYPES ***/

/* solution degrees of freedom structure (parallel to vertex structure) */
typedef struct {
	double * vec_dof_1;  /* pointer to array of field dofs */
	double * vec_dof_2;  /* pointer to array of field dofs */
	double * vec_dof_3;  /* pointer to array of field dofs */
	int * constr;
	/* for future extensions (for all mesh entities):
	int mesh_ent_type;
	int mesh_ent_id;
	int pdeg;
	*/
} apt_dof_ent;

/* discretization structure, including time and space discretization */
typedef struct {
	int     mesh_id;     /* identifier of the mesh associated with the field */
	int     nreq;        /* number of components in solution vector */
	int     nr_sol;      /* number of solution vectors for each element */
	int     uniform;     /* indicator whether the approximation field is uniform */
	int     constr;      /* indicator whether there are constrained nodes */
	int     nr_dof_ents;   /* number of dof entities */
	int capacity_dof_ents; /* max numer of dof ent. with this malloc*/
	apt_dof_ent * dof_ents;   /* pointer to structure with the vectors of dofs */
} apt_field;

/*** GLOBAL VARIABLES for the whole module ***/

extern int       apv_nr_fields;     /* the number of fields in the problem */
extern int       apv_cur_field_id;              /* ID of the current field */
extern apt_field  apv_fields[APC_MAX_NUM_FIELD];        /* array of fields */


/*** FUNCTIONS DECLARATIONS - headers for internal functions ***/
/**--------------------------------------------------------
  apr_select_field - to select the proper field
---------------------------------------------------------*/
extern apt_field * apr_select_field( /* returns: pointer to the selected field */
		int Field_id  /* in: field ID */
);

/*
In file aps_gauss_util.c - utility routines for Gauss-Legendre integration
*/

/**-------------------------------------------------------------------
  apr_gauss_init - to initialize Gauss-Legendre quadratures
--------------------------------------------------------------------*/
extern void apr_gauss_init();

/**-------------------------------------------------------------------
  apr_gauss_select - to select the proper quadrature
--------------------------------------------------------------------*/
extern void apr_gauss_select(
		int Type,	/* in: type of integration considered */
		/* 	1 - LINE INTEGRATION over [-1,1] */
		/* 	2 - GAUSSIAN QUADRATURES FOR TRIANGLES */
		int Order, 	/* in: order of integration == degree */
		/*     of polynomials integrated exactly */
		int * Ng,	/* out: number of gauss points */
		double ** Xg_p, 	/* out: gauss points */
		double ** Wg_p 	/* out: gauss weights */
);

/** @}*/

#ifdef __cplusplus
}
#endif

#endif //!header guard
