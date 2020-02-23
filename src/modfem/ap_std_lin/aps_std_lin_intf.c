/************************************************************************
File aps_std_prism_intf.c - implementation of the approximation module
				  interface for the standard linear discretization
				  of PDs on 40D meshes with prismatic elements

Contains definitions of routines:
  apr_module_introduce - to return the approximation method's name
  apr_init_field - to initiate new approximation field and read its data
  apr_write_field - to dump-out field data in the standard MOD_FEM format
  apr_check_field - to check approximation field data structure
  apr_get_mesh_id - to return the ID of the associated mesh
  apr_get_mixed_second_field_id - to return ID of second field in mixed approximation
  apr_get_nreq - to return the number of components in solution vector
  apr_get_mixed_nreq - to return the number of components in solution vector for mixed field
  apr_get_nr_sol - to return the number of solution vectors stored
  apr_get_base_type - to return the type of basis functions
  apr_get_ent_pdeg - to return the degree of approximation symbol
					  associated with a given mesh entity
  apr_set_ent_pdeg - to set the degree of approximation index
					  associated with a given mesh entity
  apr_get_ent_numshap - to return the number of shape functions (vector DOFs)
						associated with a given mesh entity
  apr_get_ent_nrdofs - to return the number of dofs associated with
					  a given mesh entity
  apr_get_el_pdeg - to return the degree of approximation vector
					  associated with a given element
  apr_set_el_pdeg - to set the degree of approximation vector
					  associated with a given element
  apr_get_mixed_pdeg - to get approximation pdeg associated with
					   first or second field in mixed approximation
  apr_get_el_geo_pdeg - to return the degree of approximation
										  associated with geometric transformation of element
  apr_set_el_geo_pdeg - to set the degree of approximation
										  associated with geometric transformation of element
  apr_get_el_geo_dofs - to get geometric coordinates for a given element
  apr_move_node - to move real geometric node by vector,
							  node is related to vertex
  apr_node_coor - to return node cordinates related to vertices
  apr_get_el_pdeg_numshap - to return the number of shape functions
							(scalar DOFs) for an element given its
							degree of approximation symbol or vector pdeg
  apr_get_el_pdeg_numshap_related_with_dof_type - return number of shape functions retated with
												  dof type
  apr_get_el_dofs - to return the number and the list of element's degrees
						of freedom (corresponding to standard shape functions)
  apr_get_nrdofs_glob - to return a global dimension of the problem
  apr_read_ent_dofs - to read a vector of dofs associated with a given
				  mesh entity from approximation field data structure
  apr_write_ent_dofs - to write a vector of dofs associated with a given
				   mesh entity to approximation field data structure
  apr_create_ent_dofs - to write a vector of dofs associated with a given
				   mesh entity to approximation field data structure
  apr_set_ini_con - to set an initial condition
  apr_num_int_el - to perform numerical integration for an element
  apr_num_int_fa - to perform numerical integration for a face
  apr_get_stiff_mat_data - to return data on dof entities for an element and
					  to compute or rewrite element's stiffness matrix and RHSV
  apr_proj_dof_ref - to rewrite dofs after modifying the mesh
  apr_rewr_sol - to rewrite solution from one vector to another
  apr_free_field - to free approximation field data structure

  apr_limit_deref - to return whether derefinement is allowed or not
  apr_limit_ref - to return whether refinement is allowed or not
	the last two routines does not use approximation data structures but
	the result depends on the approximation method
  apr_refine - to refine an element or the whole mesh checking mesh irregularity
  apr_derefine - to derefine an element or the whole mesh with irregularity check

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
/* interface of the mesh manipulation module */
#include <modfem/mmh_intf.h>

/* interface for all approximation modules */
// moved to local header where we redifne some constants
//#include "aph_intf.h"

#include <modfem/uth_intf.h>


/* internal header file for the std lin approximation module */
#include <modfem/ap_std_lin/aph_std_lin.h>

#include <modfem/uth_log.h>


#ifdef __cplusplus
extern "C" {
#endif


/** \addtogroup APM_STD_LIN Standard Linear Approximation
 *  \ingroup APR
 *  @{
 */

#define ut_min(x,y) ((x)<(y)?(x):(y))

#ifndef NULL
#define NULL NULL
#endif

#ifndef SAFE_FREE
#define SAFE_FREE(ptr) if(ptr!=NULL) free(ptr); ptr=NULL;
#endif

/* two functions needed from problem dependent module */
/*------------------------------------------------------------
  pdr_select_el_coeff_vect - to select coefficients returned to approximation
						routines for element integrals in weak formulation
		   (the procedure indicates which terms are non-zero in weak form)
------------------------------------------------------------*/
int pdr_select_el_coeff_vect( // returns success indicator
		int Problem_id,
		int * Coeff_vect_ind	/* out: coefficient indicator */
		// input to pdr_select_el_coeff_vect:
		// 0 - perform substitutions in coeff_vect_ind only if coeff_vect_ind[0]==0
		// output from pdr_select_el_coeff_vect:
		// 0 - coeff_vect_ind[0]==1 - always
		// 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz,
		// 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
		// 14 - tx, 15 - ty, 16 - tz, 17 - cval
		// 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval
);


/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
/*------------------------------------------------------------
  pdr_select_el_coeff - to select coefficients returned to approximation
						routines for element integrals
------------------------------------------------------------*/
double * pdr_select_el_coeff( /* returns: pointer !=NULL to indicate selection */
		int Problem_id,
		double ** Mval,	/* out: mass matrix coefficient */
		double ** Axx, double ** Axy, double ** Axz, /* out:diffusion coefficients, e.g.*/
		double ** Ayx, double ** Ayy, double ** Ayz, /* Axy denotes scalar or matrix */
		double ** Azx, double ** Azy, double ** Azz, /* related to terms with dv/dx*du/dy */
		/* second order derivatives in weak formulation (scalar for scalar problems */
		/* matrix for vector problems) */
		double ** Bx, double ** By, double ** Bz,	/* out: convection coefficients */
		/* Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
		double ** Tx, double ** Ty, double ** Tz,	/* out: convection coefficients */
		/* Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
		double ** Cval, /* out: reaction coefficients - for terms without derivatives */
		/*  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
		double ** Lval, /* out: rhs coefficient for time term, Lval denotes scalar */
		/* or matrix corresponding to time derivative - similar as mass matrix but  */
		/* with known solution at the previous time step (usually denoted by u_n) */
		double ** Qx, /* out: rhs coefficients for terms with derivatives */
		double ** Qy, /* Qy denotes scalar or matrix corresponding to terms with dv/dy */
		double ** Qz, /* derivatives in weak formulation */
		double ** Sval	/* out: rhs coefficients without derivatives (source terms) */
);

/*---------------------------------------------------------
  pdr_el_coeff - to return coefficients for element integrals
----------------------------------------------------------*/
int pdr_el_coeff(
		int Problem_id,
		int Elem,	/* in: element number */
		int Mat_num,	/* in: material number */
		double Hsize,	/* in: size of an element */
		int Pdeg,	/* in: local degree of polynomial */
		double * X_loc,     /* in: local coordinates of point within element */
		double * Base_phi,  /* in: basis functions */
		double * Base_dphix, /* in: x-derivatives of basis functions */
		double * Base_dphiy, /* in: y-derivatives of basis functions */
		double * Base_dphiz, /* in: z-derivatives of basis functions */
		double * Xcoor,	/* in: global coordinates of a point */
		double * Uk_val, 	/* in: computed solution from previous iteration */
		double * Uk_x, 	/* in: gradient of computed solution Uk_val */
		double * Uk_y,   	/* in: gradient of computed solution Uk_val */
		double * Uk_z,   	/* in: gradient of computed solution Uk_val */
		double * Un_val, 	/* in: computed solution from previous time step */
		double * Un_x, 	/* in: gradient of computed solution Un_val */
		double * Un_y,   	/* in: gradient of computed solution Un_val */
		double * Un_z,   	/* in: gradient of computed solution Un_val */
		double * Mval,	/* out: mass matrix coefficient */
		double * Axx, double * Axy, double * Axz, /* out:diffusion coefficients */
		double * Ayx, double * Ayy, double * Ayz, /* e.g. Axy denotes scalar or matrix */
		double * Azx, double * Azy, double * Azz, /* related to terms with dv/dx*du/dy */
		/* second order derivatives in weak formulation (scalar for scalar problems */
		/* matrix for vector problems) */
		double * Bx, double * By, double * Bz,	/* out: convection coefficients */
		/* Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
		double * Tx, double * Ty, double * Tz,	/* out: convection coefficients */
		/* Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
		double * Cval,	/* out: reaction coefficients - for terms without derivatives */
		/*  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
		double * Lval,	/* out: rhs coefficient for time term, Lval denotes scalar */
		/* or matrix corresponding to time derivative - similar as mass matrix but  */
		/* with known solution at the previous time step (usually denoted by u_n) */
		double * Qx, /* out: rhs coefficients for terms with derivatives */
		double * Qy, /* Qy denotes scalar or matrix corresponding to terms with dv/dy */
		double * Qz, /* derivatives in weak formulation */
		double * Sval	/* out: rhs coefficients without derivatives (source terms) */
);


/* end of forward declaration of functions from problem dependent module */


/*** GLOBAL VARIABLES for the whole module ***/

int       apv_nr_fields = 0;   /* the number of fields in the problem */
int       apv_cur_field_id;              /* ID of the current field */
apt_field  apv_fields[APC_MAX_NUM_FIELD];        /* array of fields */


/*------------------------------------------------------------
  apr_module_introduce - to return the approximation method's name
------------------------------------------------------------*/
int apr_module_introduce(
		/* returns: >=0 - success code, <0 - error code */
		char * Approx_name /* out: the name of the approximation method */
)
{

	char * string = "STANDARD_LINEAR_WITH_CONSTRAINED_NODES";

	strcpy(Approx_name, string);

	return (1);
}



/*------------------------------------------------------------
  apr_init_field - to initiate new approximation field and
				   read its control parameters
------------------------------------------------------------*/
int apr_init_field(  /* returns: >0 - field ID, <0 - error code */
		char Field_type,
		int Control,	 /* in: control variable: */
		/*      APC_ZERO - to initialize the field to zero */
		/*      APC_READ - to read field values from the file */
		/*      APC_INIT - to initialize the field using function */
		/*                 provided by the problem dependent module */
		int Mesh_id,	 /* in: ID of the corresponding mesh */
		int Nreq,	 /* in: number of equations - solution vector components */
		int Nr_sol,	 /* in: number of solution vectors for each dof entity */
		int Pdeg_in,   /* in: degree of approximating polynomial and geometry in format */
		/*     [psg] p - primary field degree / s - secondary field degree for mixe default 0 / g - geometry */
		//int GeoPdeg_in, /* in: degree of approximating polynomial for geometry */
		char * Filename, /* in: name of the file to read approximation data */
		double (*Fun_p)(int, double *, int) /* pointer to function that provides */
		/* problem dependent initial condition data */
)
{

	double read_v;
	int iter_nodes;
	/* pointer to field structure */

	apt_field * field_p;
	int mxvert, nvert, ivert, nreq, nrvert_glob, nr_dof_ents; //num_dof,

	/* auxiliary variables */
	int i; //iaux,,j,k,jaux
	FILE * fp;
	char mesh_name[100];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* first check that the mesh is compatible */
	mmr_module_introduce(mesh_name);

#ifdef DEBUG_APM
	if (strcmp("3D_PRISM", mesh_name) != 0 && strcmp("3D_Hybrid", mesh_name) != 0 ) {
		// if(monitor>APC_NO_PRINT){
		printf("Wrong mesh %s for standard linear approximation with prisms!\n",
				mesh_name);
		// {
		//exit(-1);
	}
#endif
	/*  */

	/* initiate quadrature data for Gauss-Legendre integration */
	apr_gauss_init();

	/* increase the counter for fields */
	apv_nr_fields++;

	/* set the current field ID */
	apv_cur_field_id = apv_nr_fields;

	field_p = &apv_fields[apv_cur_field_id - 1];

	field_p->mesh_id = Mesh_id;

	/* at the beginning we assume uniform mesh without constrained nodes */
	field_p->uniform = APC_TRUE;
	field_p->constr = APC_FALSE;

	field_p->nreq = Nreq;
	field_p->nr_sol = Nr_sol;

	/* depending on the control variable */
	if (Control == APC_READ) {

		/* open the input file with approximation parameters */
		fp = fopen(Filename, "r");
#ifdef DEBUG_APM
		if (fp == NULL) {
			printf("Not found file '%s' with field data!!! Exiting.\n", Filename);
			exit(-1);
		}
#endif


		/* read the number of solution vectors for each dof entity and */
		/* the number of components in each vector = number of equations */
		fscanf(fp, "%d %d",
				&field_p->nreq, &field_p->nr_sol);

#ifdef DEBUG_APM
		if (field_p->nreq != Nreq) {
			printf("Requested number of equations %d different than read from file %d!!!\n",
					Nreq, field_p->nreq);
			printf("Requested number of components %d different than read from file %d!!!\n",
					Nr_sol, field_p->nr_sol);
			printf("Exiting.\n");
			exit(-1);
		}
#endif
	}
	else {

	}
	nreq = field_p->nreq;

	if (field_p->nreq > APC_MAXEQ) {
		printf("Requested number of equations %d greater than APC_MAXEQ=%d!!!\n",
				field_p->nreq, APC_MAXEQ);
		printf("Change APC_MAXEQ in aph_std_lin_prism.h and recompile the code.\n");
		printf("Exiting.\n");
		exit(-1);
	}

	/* get the maximal admissible number of vertices */
	mxvert = mmr_get_max_node_max(Mesh_id);
	field_p->capacity_dof_ents = mxvert;
	/* allocate the space for degrees of freedom */
	field_p->dof_ents = (apt_dof_ent *) malloc((mxvert + 1) * sizeof(apt_dof_ent));
	if (field_p->dof_ents == NULL) {
		printf("Dofs structures not allocated\n");
		exit(-1);
	}

	for (nvert = 1; nvert <= mxvert; nvert++) {
		field_p->dof_ents[nvert].vec_dof_1 = NULL;
		field_p->dof_ents[nvert].vec_dof_2 = NULL;
		field_p->dof_ents[nvert].vec_dof_3 = NULL;
		field_p->dof_ents[nvert].constr = NULL;
	}

	/* depending on the control variable */
	if (Control == APC_READ) {

		fscanf(fp, "%d", &nrvert_glob);

		assert(nrvert_glob <= mmr_get_max_node_id(Mesh_id));
		nr_dof_ents = nrvert_glob;

	}
	else {

		nrvert_glob = mmr_get_max_node_id(Mesh_id);
		/*kbw
		printf("nrvert_glob from mesh %d\n", nrvert_glob);
		/*kew*/
	}


	/* depending on the control variable */
	if (Control == APC_READ) {

		/* create dofs data structure */
		for (ivert = 1; ivert <= nrvert_glob; ivert++) {

			fscanf(fp, "%d", &nvert);

			field_p->dof_ents[nvert].vec_dof_1 = (double *) malloc(nreq * sizeof(double));

#ifdef DEBUG_APM
			if (field_p->dof_ents[nvert].vec_dof_1 == NULL) {
				printf("Dofs vector not allocated\n");
				exit(-1);
			}
#endif

			if (field_p->nr_sol > 1) {
				field_p->dof_ents[nvert].vec_dof_2 =
						(double *) malloc(nreq * sizeof(double));
			}
			if (field_p->nr_sol > 2) {
				field_p->dof_ents[nvert].vec_dof_3 =
						(double *) malloc(nreq * sizeof(double));

			}


			/* read values from a dump file */
			/*kbw
				  printf("ivert %d, nvert %d, nrvert_glob %d\n",
					 ivert, nvert, nrvert_glob);
			/*kew*/

			for (i = 0; i < nreq; i++) {

				fscanf(fp, "%lg", &read_v);
				field_p->dof_ents[nvert].vec_dof_1[i] = read_v;

				assert(field_p->dof_ents[nvert].vec_dof_1[i] > -10e200 &&
						field_p->dof_ents[nvert].vec_dof_1[i] < 10e200);

				/*kbw
					printf("ieq %d, sol_1 %lf\n", i,
						 field_p->dof_ents[nvert].vec_dof_1[i]);
				/*kew*/

			}
			if (field_p->nr_sol > 1) {
				for (i = 0; i < nreq; i++) {
					fscanf(fp, "%lg", &read_v);
					field_p->dof_ents[nvert].vec_dof_2[i] = read_v;

					assert(field_p->dof_ents[nvert].vec_dof_2[i] > -10e200 &&
							field_p->dof_ents[nvert].vec_dof_2[i] < 10e200);

					/*kbw
						  printf("ieq %d, sol_2 %lf\n", i,
							 field_p->dof_ents[nvert].vec_dof_2[i]);
					/*kew*/

				}
			}
			if (field_p->nr_sol > 2) {
				for (i = 0; i < nreq; i++) {
					fscanf(fp, "%lg", &read_v);
					field_p->dof_ents[nvert].vec_dof_3[i] = read_v;

					assert(field_p->dof_ents[nvert].vec_dof_3[i] > -10e200 &&
							field_p->dof_ents[nvert].vec_dof_3[i] < 10e200);

					/*kbw
						  printf("ieq %d, sol_3 %lf\n", i,
							 field_p->dof_ents[nvert].vec_dof_3[i]);
					/*kew*/

				}
			}


		} // end loop over read nodes

	} /* end if read values from a dump file */
	else {
		/* iniiate the field to zero */

		nr_dof_ents = 0;
		/* create dofs data structure */
		for (nvert = 1; nvert <= nrvert_glob; nvert++) {

			/* for active unconstrained vertices */
			if (mmr_node_status(Mesh_id, nvert) == MMC_ACTIVE) {

				nr_dof_ents++;

				field_p->dof_ents[nvert].vec_dof_1 =
						(double *) malloc(nreq * sizeof(double));

#ifdef DEBUG_APM
				if (field_p->dof_ents[nvert].vec_dof_1 == NULL) {
					printf("Dofs vector not allocated\n");
					exit(-1);
				}
#endif

				if (field_p->nr_sol > 1) {
					field_p->dof_ents[nvert].vec_dof_2 =
							(double *) malloc(nreq * sizeof(double));
				}
				if (field_p->nr_sol > 2) {
					field_p->dof_ents[nvert].vec_dof_3 =
							(double *) malloc(nreq * sizeof(double));
				}


				/*kbw
					  printf("setting ivert %d, nrvert_glob %d\n", nvert, nrvert_glob);
				/*kew*/

				for (i = 0; i < nreq; i++) {
					field_p->dof_ents[nvert].vec_dof_1[i] = 0.0;
				}
				if (field_p->nr_sol > 1) {
					for (i = 0; i < nreq; i++) {
						field_p->dof_ents[nvert].vec_dof_2[i] = 0.0;
					}
				}
				if (field_p->nr_sol > 2) {
					for (i = 0; i < nreq; i++) {
						field_p->dof_ents[nvert].vec_dof_3[i] = 0.0;
					}
				}
			}

		}

	} /* end if not read from file */

	/* number of dof entities does not include constrained nodes ! */
	field_p->nr_dof_ents = nr_dof_ents;

	/*jbw*/
	printf("\n");
	printf("Approximation type: LINEAR\n");
	printf("NR DOF ENTS/VERTs: %d ; nreq = %d ; nr_sol = %d\n", field_p->nr_dof_ents, field_p->nreq, field_p->nr_sol);
	/*jbw*/

	/* close file with approximation data */
	if (Control == APC_READ) {
		fclose(fp);
	}

	/* initiate values using provided function */
	if (Control == APC_INIT) {

		apr_set_ini_con(apv_cur_field_id, Fun_p);

#ifdef DEBUG_APM
		printf("\nSpecified initial data for approximation field %d!\n",
				apv_cur_field_id);
#endif

	}

	#pragma omp critical(creating_constraints)
	{
		/*initiate constraints data*/
		apr_create_constr_data(apv_cur_field_id);
	}

	/*jbw*/
	printf("After create CONSTRAINTS DATA:\n");
	printf("NR DOF ENTS/VERTs: %d ; nreq = %d ; nr_sol = %d\n", field_p->nr_dof_ents, field_p->nreq, field_p->nr_sol);
	printf("\n");
	/*jbw*/

#ifdef DEBUG_APM
	apr_check_field(apv_cur_field_id);
#endif

	return (apv_cur_field_id);
}

/*---------------------------------------------------------
  apr_write_field - to dump-out field data in the standard HP_FEM format
---------------------------------------------------------*/
int apr_write_field( /* returns: >=0 - success code, <0 - error code */
		int Field_id,    /* in: field ID */
		int Nreq,        /* in: number of equations (scalar dofs) */
		int Select,      /* in: parameter to select written vectors */
		double Accuracy, /* in: parameter specyfying accuracy  */
		/* (values below Accuracy are written as 0.0) */
		/* (put 0.0 for full accuracy in "%g" format) */
		char * Filename  /* in: name of the file to write field data */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id, nmvert, nvert, nreq, nr_dof_ents, i, selected[3], nr_sol;
	FILE * fp;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	/* open the input file with approximation parameters */
	fp = fopen(Filename, "w");

	/* check which vectors to write to file */
	if (Select <= 0 || Select >= 7) {
		selected[0] = 1;
		selected[1] = 1;
		selected[2] = 1;
		nr_sol = field_p->nr_sol;
	}
	else {
		// scheme: 1 - 1, 2 - 2, 3 - 4, 1+2 - 3, 1+3 - 5, 2+3 - 6, 1+2+3 - 7
		selected[0] = Select % 2;
		selected[1] = (Select / 2) % 2;
		selected[2] = (Select / 4) % 2;
		nr_sol = selected[0] + selected[1] + selected[2];
	}

	/* write the number of solution vectors for each dof entity */
	mf_log_info("writing field %d %d into file %s \n",
			field_p->nreq, nr_sol, Filename);
	fprintf(fp, "%d %d\n",
			field_p->nreq, nr_sol);
	nreq = apr_get_nreq(Field_id);

	/* loop over dofs data structure including free spaces */
	nmvert = mmr_get_max_node_id(mesh_id);

#ifdef DEBUG_APM
	nr_dof_ents = 0;
	for (nvert = 1; nvert <= nmvert; nvert++) {
		/* for active unconstrained vertices */
		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE &&
				field_p->dof_ents[nvert].vec_dof_1 != NULL) nr_dof_ents++;
	}
	if (nr_dof_ents != field_p->nr_dof_ents) {
		printf("wrong number of dof_entities in apr_write_field! %d %d Exiting!\n", nr_dof_ents, field_p->nr_dof_ents);
		exit(-1);
	}
#endif

	nr_dof_ents = field_p->nr_dof_ents;

	/* write the number of dof entities */
	fprintf(fp, "%d\n", nr_dof_ents);


	// find maximal DOF
	double DOF_max[APC_MAXEQ] = {0.0};
	for (nvert = 1; nvert <= nmvert; nvert++) {

		/* for active vertices */
		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE &&
				field_p->dof_ents[nvert].vec_dof_1 != NULL) {

			for (i = 0; i < nreq; i++) {

				if (fabs(field_p->dof_ents[nvert].vec_dof_1[i]) > DOF_max[i]) {
					DOF_max[i] = fabs(field_p->dof_ents[nvert].vec_dof_1[i]);
				}

			}
		}
	}


	for (nvert = 1; nvert <= nmvert; nvert++) {

		/* for active vertices */
		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE &&
				field_p->dof_ents[nvert].vec_dof_1 != NULL) {

			fprintf(fp, "%d ", nvert);

			if (selected[0] == 1) {
				for (i = 0; i < nreq; i++) {

					if (fabs(field_p->dof_ents[nvert].vec_dof_1[i]) < Accuracy * DOF_max[i]) {
						double temp = 0.0;
						fprintf(fp, "%.1lg ", temp);
					}
					else {
						fprintf(fp, "%.12lg ", field_p->dof_ents[nvert].vec_dof_1[i]);
					}
					/*
					utr_fprintf_double(fp, Accuracy, field_p->dof_ents[nvert].vec_dof_1[i]);
					fprintf(fp," ");
					*/

				}
				fprintf(fp, "\n");
			}

			if (field_p->nr_sol > 1 && selected[1] == 1) {

				for (i = 0; i < nreq; i++) {

					if (fabs(field_p->dof_ents[nvert].vec_dof_2[i]) < Accuracy * DOF_max[i]) {
						double temp = 0.0;
						fprintf(fp, "%.12lg ", temp);
					}
					else {
						fprintf(fp, "%.12lg ", field_p->dof_ents[nvert].vec_dof_2[i]);
					}
					/*
						utr_fprintf_double(fp, Accuracy, field_p->dof_ents[nvert].vec_dof_2[i]);
						fprintf(fp," ");
						*/
				}
				fprintf(fp, "\n");

			}
			if (field_p->nr_sol > 2 && selected[2] == 1) {

				for (i = 0; i < nreq; i++) {

					if (fabs(field_p->dof_ents[nvert].vec_dof_3[i]) < Accuracy * DOF_max[i]) {
						double temp = 0.0;
						fprintf(fp, "%.12lg ", temp);
					}
					else {
						fprintf(fp, "%.12lg ", field_p->dof_ents[nvert].vec_dof_3[i]);
					}
					/*
					utr_fprintf_double(fp, Accuracy, field_p->dof_ents[nvert].vec_dof_3[i]);
					fprintf(fp," ");
					*/

				}
				fprintf(fp, "\n");

			}
		} /* end if vertex active */
	} /* end loop over all dofs structures */


	/* close file with control data */
	fclose(fp);

	return (1);
}


/*------------------------------------------------------------
  apr_check_field - to free approximation field data structure
------------------------------------------------------------*/
int apr_check_field(
		int Field_id    /* in: approximation field ID  */
)
{
	/* pointer to field structure */
	apt_field * field_p;
	int nvert, nmvert, mesh_id, nr_dof_ents; //i,

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	/* loop over dofs data structure including free spaces */
	nmvert = mmr_get_max_node_id(mesh_id);

	/*kbw
	printf("Checking field %d associated with mesh %d, nmvert=%d\n",
	 Field_id, mesh_id, nmvert);
	/*kew*/

	assert(nmvert >= field_p->nr_dof_ents);
	for (nvert = 1; nvert <= nmvert; nvert++) {

		//assert(nvert != 15344);

		if (mmr_node_status(mesh_id, nvert) == MMC_FREE) {

			if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {
				printf("Dofs associated with free space %d in check_field, exiting\n",
						nvert);
				exit(-1);
			}
			if (field_p->dof_ents[nvert].constr != NULL) {
				printf("Constr associated with free space %d in check_field, exiting\n",
						nvert);
				exit(-1);
			}
			if (apr_get_ent_pdeg(Field_id, APC_VERTEX, nvert) != -1) {
				printf("Error 3843 in approximation data structures for %d. Exiting!\n",
						nvert);
				//assert(!"Holding execution before exit!");
				exit(-1);
			}


		}
		else {

			if (apr_get_ent_pdeg(Field_id, APC_VERTEX, nvert) == 1) {

				/*kbw
					printf("Real node %d, dof %lf\n",
					   nvert, field_p->dof_ents[nvert].vec_dof_1[0]);
				/*kew*/

				if (field_p->dof_ents[nvert].vec_dof_1 == NULL) {
					printf("Error 38 in approximation data structures for %d. Exiting!\n",
							nvert);
					exit(-1);
				}

				if (field_p->dof_ents[nvert].constr != NULL) {
					printf("Constr associated with real node %d in check_field, exiting\n",
							nvert);
					exit(-1);
				}

			}
			else if (apr_get_ent_pdeg(Field_id, APC_VERTEX, nvert) == 0) {

				/*kbw
					printf("Constrained node %d, parents %d %d %d %d\n",
					   nvert, nvert, nvert, nvert, nvert);
				/*kew*/


				if (field_p->dof_ents[nvert].constr == NULL) {
					assert(field_p->dof_ents[nvert].constr != NULL);
					assert(field_p->dof_ents[nvert].constr != NULL);
					printf("Error 43 in approximation data structures for %d in check field. Exiting!\n",
							nvert);
					exit(-1);
				}

				if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {
					printf("Dofs associated with constr node %d in check_field, exiting\n",
							nvert);
					exit(-1);
				}
			}

		}

	}

	nr_dof_ents = 0;
	for (nvert = 1; nvert <= nmvert; nvert++) {
		/* for active unconstrained vertices */
		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE &&
				field_p->dof_ents[nvert].vec_dof_1 != NULL) nr_dof_ents++;
	}
	if (nr_dof_ents != field_p->nr_dof_ents) {
		printf("wrong number of dof_entities in apr_check_field! %d %d Exiting!\n", nr_dof_ents, field_p->nr_dof_ents);
		exit(-1);
	}


	return (1);
}

/*------------------------------------------------------------
  apr_get_mesh_id - to return the ID of the associated mesh
------------------------------------------------------------*/
int apr_get_mesh_id( /* returns: >0 - ID of the associated mesh,
										<0 - error code */
		int Field_id     /* in: approximation field ID  */
)
{
	/* pointer to field structure */
	apt_field * field_p;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	mfp_check_debug(field_p->mesh_id > 0, "Wrong mesh id in field!");

	return (field_p->mesh_id);
}

/**-----------------------------------------------------------
apr_get_mixed_second_field_id - to return ID of second field in mixed approximation
------------------------------------------------------------*/
int apr_get_mixed_second_field_id( /** returns: >0 - ID of the associated second field id in mixed,
						  <0 - error code */
		int Field_id     /** in: approximation primary field ID  */
)
{
	/* pointer to field structure */
	apt_field * field_p = NULL;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	return (-1);
}

/*------------------------------------------------------------
  apr_get_nreq - to return the number of components in solution vector
------------------------------------------------------------*/
int apr_get_nreq( /* returns: >0 - number of solution components,
										<0 - error code */
		int Field_id     /* in: approximation field ID  */
)
{
	/* pointer to field structure */
	apt_field * field_p;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	return (field_p->nreq);
}

/*------------------------------------------------------------
  apr_get_mixed_nreq - to return the number of components in solution vector for mixed field
------------------------------------------------------------*/
int apr_get_mixed_nreq( /* returns: >0 - number of solution components,
										<0 - error code */
		int Field_id     /* in: approximation field ID  */
)
{

	/* pointer to field structure */
	apt_field * field_p = NULL;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	return (-1);
}

/*------------------------------------------------------------
  apr_get_nr_sol - to return the number of solution vectors stored
------------------------------------------------------------*/
int apr_get_nr_sol( /* returns: >0 - number of solution vectors,
										<0 - error code */
		int Field_id     /* in: approximation field ID  */
)
{
	/* pointer to field structure */
	apt_field * field_p;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	return (field_p->nr_sol);
}



/*---------------------------------------------------------
  apr_get_base_type - to return the type of basis functions
	 REMARK: type of basis functions differentiates element types as well
		examples for standard linear approximation (from include/aph_intf.h):
		  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements
		  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements
---------------------------------------------------------*/
int apr_get_base_type(/* returns: >0 - type of basis functions,
								  <0 - error code */
		int Field_id,  /* in: field ID */
		int El_id      /* in: element ID */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id, el_type;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* get the corresponding mesh ID */
	mesh_id = field_p->mesh_id;

// base type is used to differentiate element types
// get element type
	el_type = mmr_el_type(mesh_id, El_id);
	if (el_type == MMC_PRISM) {
		return (APC_BASE_PRISM_STD);
	} else {
		return (APC_BASE_TETRA_STD);
	}

	printf("Wrong element type in apr_get_base_type. Exiting!\n");
	exit(-1);
}


/*------------------------------------------------------------
  apr_get_ent_pdeg - to return the degree of approximation index
					  associated with a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_pdeg( /* returns: >0 - approximation index,
								   0 - dof entity inactive (constrained)
								  <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_type,      /* in: type of mesh entity */
		int Ent_id         /* in: mesh entity ID */
)
{
	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;

	/* auxiliary variables */
	int pdeg;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
	/* check input for standard approximation */
	if (Ent_type != APC_VERTEX) {
		assert(!"Holding execution.");
		printf("Wrong dof entity type in get_ent_pdeg !\n");
		exit(-1);
	}
#endif

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* get the corresponding mesh ID */
	mesh_id = field_p->mesh_id;

	if (mmr_node_status(mesh_id, Ent_id) != MMC_ACTIVE) pdeg = -1;
	else {
		if (field_p->dof_ents[Ent_id].vec_dof_1 == NULL) pdeg = 0;
		else pdeg = 1;
	}

	return (pdeg);

}


/*------------------------------------------------------------
  apr_set_ent_pdeg - to set the degree of approximation index
					  associated with a given mesh entity
------------------------------------------------------------*/
int apr_set_ent_pdeg( /* returns: >0 - success code,
										  <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_type,      /* in: type of mesh entity */
		int Ent_id,         /* in: mesh entity ID */
		int Pdeg          /* in: degree of approximation */
)
{
	/* pointer to field structure */
//  apt_field *field_p;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	return (1);

}

/*------------------------------------------------------------
  apr_get_ent_numshap - to return the number of shape functions (vector
						dofs) associated with a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_numshap( /* returns: >0 - the number of shape functions,
										  <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_type,      /* in: type of mesh entity */
		int Ent_id         /* in: mesh entity ID */
)
{
	/* auxiliary variables */
//  int nrdof;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
	/* check input for standard approximation */
	if (abs(Ent_type != APC_VERTEX)) {
		printf("Wrong dof entity type in get_ent_pdeg for STD approximation !\n");
		assert(!"Halting execution!");
		exit(-1);
	}
#endif

	return (1);
}


/*------------------------------------------------------------
  apr_get_ent_nrdofs - to return the number of dofs associated with
					  a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_nrdofs( /* returns: >0 - the number of dofs,
										  <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_type,      /* in: type of mesh entity */
		int Ent_id         /* in: mesh entity ID */
)
{
	/* auxiliary variables */
//  int nrdof;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
	/* check input for standard approximation */
	if (abs(Ent_type != APC_VERTEX)) {
		printf("Wrong dof entity type in get_ent_pdeg for STD approximation !\n");
		assert(!"Holding execution!");
		exit(-1);
	}
#endif

	return (apr_get_nreq(Field_id));
}

/*------------------------------------------------------------
  apr_get_el_pdeg - to return the degree of approximation vector
					  associated with a given element
------------------------------------------------------------*/
int apr_get_el_pdeg( /* returns:  >0 - success code - scalar pdeg
										 <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_id,         /* in: element ID */
		int * Pdeg_vec      /* out: degree of approximation symbol or vector */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (Pdeg_vec != NULL) *Pdeg_vec = 101;

	return (101);

}


/*------------------------------------------------------------
  apr_set_el_pdeg - to set the degree of approximation vector
					  associated with a given element
------------------------------------------------------------*/
int apr_set_el_pdeg( /* returns: >0 - success code,
										  <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_id,         /* in: element ID */
		int * Pdeg_vec      /* in: degree of approximation symbol or vector */
)
{

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	return (1);

}

/**-----------------------------------------------------------
  apr_get_mixed_pdeg - to get approximation pdeg associated with
					   first or second field in mixed approximation
------------------------------------------------------------*/
int apr_get_mixed_pdeg( /** returns: >0 - success code,
					<0 - error code */
		int Field_id,     /** in: approximation field ID  */
		int * Pdeg_vec      /** out: degree of approximation symbol */
)
{
	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (Pdeg_vec != NULL) *Pdeg_vec = 1;

	return (1);
}

/**-----------------------------------------------------------
  apr_get_el_geo_pdeg - to return the degree of approximation
										  associated with geometric transformation of element
------------------------------------------------------------*/
int apr_get_el_geo_pdeg( /** returns: >0 - success code,
					<0 - error code */
		int Field_id,     /** in: approximation field ID  */
		int El_id,         /** in: element ID */
		int * Pdeg_vec      /** out: degree of approximation symbol */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (Pdeg_vec != NULL) *Pdeg_vec = 1;

	return (1);

}

/**-----------------------------------------------------------
  apr_set_el_geo_pdeg - to set the degree of approximation
										  associated with geometric transformation of element
------------------------------------------------------------*/
int apr_set_el_geo_pdeg( /** returns: >0 - success code,
					<0 - error code */
		int Field_id,     /** in: approximation field ID  */
		int El_id,         /** in: element ID */
		int * Pdeg_vec      /** in: degree of approximation symbol */
)
{

	/*++++++++++++++++ executable statements ++++++++++++++++*/
	return (1);
}

/**-----------------------------------------------------------
 apr_move_node - to move real geometric node by vector,
							  node is related to vertex / edge or face
------------------------------------------------------------*/
int apr_move_node( /** returns: >0 - success code,
									  <0 - error code */
		int Field_id, /* in: field ID */
		int Node_id, /* in: node ID */
		int Node_type, /* in: node type - APC_VERTEX / APC_EDGE / APC_FACE */
		double * Shift_vec /* in: shift vector */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;
	int error_code = 0;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* get mesh ID */
	mesh_id = apr_get_mesh_id(Field_id);

	if (Shift_vec != NULL) {

		if (Node_type == APC_VERTEX) {
			error_code = mmr_move_node(mesh_id, Node_id, Shift_vec);
		} else {
			mf_fatal_err("Wrong node type in apr_move_node. Exiting!!!\n");
			exit(-1);
		}

	} else {
		mf_fatal_err("Input shift vector is NULL in apr_move_node. Exiting!!!\n");
		exit(-1);
	}

	return (error_code);

}

/*-----------------------------------------------------------
  apr_node_coor - to return node cordinates related to vertices
------------------------------------------------------------*/
int apr_node_coor( /** returns: >0 - success code,
						   <0 - error code */
		int Field_id, /* in: field ID */
		int Node_id, /* in: node ID */
		int Node_type, /* in: node type - APC_VERTEX / APC_EDGE / APC_FACE */
		double * coord /* out: cordinates of node */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;
	int error_code = 0;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* get mesh ID */
	mesh_id = apr_get_mesh_id(Field_id);

	/* to get coordinates */
	switch (Node_type)
	{

		// Get cordinates for vertex
		case APC_VERTEX:
			mmr_node_coor(mesh_id, Node_id, coord);
			error_code = APC_VERTEX;
			break;

		default:
			mf_fatal_err("Unknown node type in apr_node_coor. Exiting!!!\n");
			exit(-1);
	};

	return error_code;
}

/*---------------------------------------------------------
  apr_get_el_pdeg_numshap - to return the number of shape functions
							(scalar DOFs) for an element given its
							degree of approximation symbol or vector pdeg
---------------------------------------------------------*/
int apr_get_el_pdeg_numshap(
		/* returns: >=0 - success code, <0 - error code*/
		int Field_id,  /* in: field ID */
		int El_id,     /* in: element ID */
		int * Pdeg     /* in: degree of approximation symbol or vector */
)
{

	int base = apr_get_base_type(Field_id, El_id);

	/* REMARK: type of basis functions differentiates element types as well
		 examples for standard linear approximation (from include/aph_intf.h):
			#define APC_BASE_PRISM_STD  3   // for linear prismatic elements
			#define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements
	*/

	if (base == APC_BASE_PRISM_STD) return (6);
	else if (base == APC_BASE_TETRA_STD) return (4);

	printf("Wrong element type in apr_get_base_type. Exiting!\n");
	exit(-1);
}

/*---------------------------------------------------------
apr_get_el_pdeg_numshap_related_with_dof_type - return number of shape functions retated with
												  dof type
---------------------------------------------------------*/
int apr_get_el_pdeg_numshap_related_with_dof_type(
		/* returns: >=0 - success code, <0 - error code*/
		int Field_id,  /* in: field ID */
		int El_id,     /* in: element ID */
		int Dof_type,  /* in: dof type */
		int * Pdeg     /* in: degree of approximation symbol or vector */
) {

	int num_shap = apr_get_el_pdeg_numshap(Field_id, El_id, Pdeg);
	return (num_shap);
}

/*---------------------------------------------------------
  apr_get_el_dofs - to return the number and the list of element's degrees
						of freedom (corresponding to standard shape functions)
!!! The order is the same as for mmr_el_node_coor (loop over ieq inside):
(node1_1, node1_2, ..., node1_NREQ, node2_1, node2_2, ..., node2_NREQ, etc. !!!
---------------------------------------------------------*/
int apr_get_el_dofs( /* returns: >0 - the number of dofs
										 <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int El_id,        /* in: element ID */
		int Vect_id,       /* in: vector ID in case of multiple solution vectors */
		double * El_dofs_std    /* out: the list of values of element dofs */
		/* the size of the table MUST BE >= APC_MAXELSD */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id, nreq, idof, idofent, ieq;
	int ino, icon, nr_dofs;
	int el_nodes[MMC_MAXELVNO + 1];      // list of nodes of El
	int nr_constr[MMC_MAXELVNO + 1];     // list of numbers of constraining nodes
	int constr_id[4 * MMC_MAXELVNO];     // list with IDs of constraining nodes
	double constr_value[4 * MMC_MAXELVNO]; // list with constraint coefficients
	double dofs_loc[APC_MAXEQ];
	int is_constrained;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* get mesh ID */
	mesh_id = apr_get_mesh_id(Field_id);

	nreq = apr_get_nreq(Field_id);

	for (idof = 0; idof < APC_MAXELSD; idof++) El_dofs_std[idof] = 0.0;

	/* if there are no constraints in the mesh */
	if (field_p->constr == APC_FALSE) {

		/* get vertices */
		mmr_el_node_coor(mesh_id, El_id, el_nodes, NULL);

		/* simply rewrite vertex dofs */
		for (ino = 0; ino < el_nodes[0]; ino++) {

			apr_read_ent_dofs(Field_id, APC_VERTEX, el_nodes[ino + 1],
					nreq, Vect_id, &El_dofs_std[ino * nreq]);

		}

	}
	else {

		/* if there are constraints */

		/* get the list of real nodes with constraint coefficients */
		is_constrained = apr_get_el_constr_data(Field_id, El_id, el_nodes,
						nr_constr, constr_id, NULL, constr_value);

		idofent = 0;
		/* rewrite dofs */
		for (ino = 0; ino < el_nodes[0]; ino++) {

			for (icon = 1; icon <= nr_constr[ino + 1]; icon++) {

				apr_read_ent_dofs(Field_id, APC_VERTEX, constr_id[idofent],
						nreq, Vect_id, dofs_loc);

				if (nr_constr[ino + 1] == 1) {

					for (ieq = 0; ieq < nreq; ieq++) {

						El_dofs_std[ino * nreq + ieq] += dofs_loc[ieq];

					}
				}
				else if (nr_constr[ino + 1] == 2) {

					for (ieq = 0; ieq < nreq; ieq++) {

						El_dofs_std[ino * nreq + ieq] += 0.5 * dofs_loc[ieq];

					}
				}
				else if (nr_constr[ino + 1] == 4) {

					for (ieq = 0; ieq < nreq; ieq++) {

						El_dofs_std[ino * nreq + ieq] += 0.25 * dofs_loc[ieq];

					}
				}
				idofent++;
			}
		}
	}

	nr_dofs = nreq * el_nodes[0];

	assert(El_dofs_std == NULL || (El_dofs_std[0] > -10e200 &&
					El_dofs_std[0] < 10e200));

	return (nr_dofs);
}

/**-----------------------------------------------------------
apr_get_el_geo_dofs - to get geometric coordinates for a given element
------------------------------------------------------------*/
int apr_get_el_geo_dofs( /** returns: >0 - the number of dofs
							   <0 - error code */
		int Field_id, /* in: approximation field ID  */
		int El_id,    /* in: element ID */
		int * GeoDofs_id,   /* out: list of dofs in element */
		/*	 (Dofs[0] - number of dofs) */
		int * GeoDofs_type, /* out: type of dofs [APC_VERTEX / APC_EDGE] */
		/*     (Dofs_type[0] - geo_pdeg */
		double * GeoDofs /* out: coordinates of element dofs */
)
{
	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;
	int num_vdofs;
	int i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/
	assert(El_id > 0);

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);
	num_vdofs = 0;

	/* get the corresponding mesh ID */
	mesh_id = field_p->mesh_id;

	/* get coordinates of nodes */
	num_vdofs = mmr_el_node_coor(mesh_id, El_id, GeoDofs_id, GeoDofs); // Get ids and geo coord
	if (GeoDofs_type != NULL) {
		for (i = 1; i <= num_vdofs; i++) {
			GeoDofs_type[i] = APC_VERTEX;    // Write types
		}
		GeoDofs_type[0] = APC_LINEAR_APPROXIMATION_PDEG;
	}

	return num_vdofs;
}


/*---------------------------------------------------------
  apr_get_nrdofs_glob - to return a global dimension of the problem
---------------------------------------------------------*/
int apr_get_nrdofs_glob(	/* returns: global dimension of the problem */
		int Field_id    /* in: field ID */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;

	/* auxiliary variables */
	int glob_dim = 0, nno;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

#ifdef DEBUG_APM
	/* get the corresponding mesh ID */
	mesh_id = field_p->mesh_id;

	/* loop over vertices-nodes */
	nno = 0;
	while ((nno = mmr_get_next_node_all(mesh_id, nno)) != 0)
	{
		if (mmr_node_status(mesh_id, nno) == MMC_ACTIVE
				&& field_p->dof_ents[nno].vec_dof_1 != NULL)
			glob_dim ++;
	}

	if (glob_dim != field_p->nr_dof_ents) {
		printf("wrong number of dof_entities in get_nrdofs_glob! %d %d Exiting!\n", glob_dim, field_p->nr_dof_ents);
		exit(-1);
	}
#endif

	mf_debug("field_p->nr_dof_ents=%d, apr_get_nreq(Field_id)=%d",
			field_p->nr_dof_ents, apr_get_nreq(Field_id));

	return (field_p->nr_dof_ents * apr_get_nreq(Field_id));
}

/*------------------------------------------------------------
  apr_read_ent_dofs - to read a vector of dofs associated with a given
				  mesh entity from approximation field data structure
------------------------------------------------------------*/
int apr_read_ent_dofs(/* returns: >=0 - success code, <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_type,      /* in: type of mesh entity */
		int Ent_id,        /* in: mesh entity ID */
		int Ent_nrdof,     /* in: number of dofs associated with the entity */
		int Vect_id,       /* in: vector ID in case of multiple solution vectors */
		double * Vect_dofs /* out: dofs read from data structure */
)
{
	/* pointer to field structure */
	apt_field * field_p;

	/* auxiliary variables */
	double * glob_dofs;
	int i, mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* get the corresponding mesh ID */
	mesh_id = field_p->mesh_id;


#ifdef DEBUG_APM
	{
		/* check input */
		if (field_p->dof_ents[Ent_id].vec_dof_1 == NULL) {
			printf("Wrong dof entity %d in read_ent_dofs !\n", Ent_id);
			getchar();
			getchar();
			exit(-1);
		}
	}
#endif

	if (Vect_id <= 1) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
	else if (Vect_id == 2) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_2;
	else if (Vect_id == 3) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_3;
	else field_p->dof_ents[Ent_id].vec_dof_1;

	assert(Vect_dofs != NULL);

#ifdef DEBUG_APM
	if (glob_dofs == NULL) {
		printf("no glob_dofs for field %d, ent_type %d, ent_nrdof %d, vect_id %d\n",
				Field_id, Ent_type, Ent_nrdof, Vect_id);
		printf("dof entity %d in read_ent_dofs !\n", Ent_id);
		getchar();
		getchar();
		exit(-1);
	}
#endif

	assert(glob_dofs != NULL);

	for (i = 0; i < Ent_nrdof; i++) Vect_dofs[i] = glob_dofs[i];

	assert((Vect_dofs[0] > -10e200 && Vect_dofs[0] < 10e200));

	return (1);
}

/*------------------------------------------------------------
  apr_write_ent_dofs - to write a vector of dofs associated with a given
				   mesh entity to approximation field data structure
------------------------------------------------------------*/
int apr_write_ent_dofs(/* returns: >=0 - success code, <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_type,      /* in: type of mesh entity */
		int Ent_id,        /* in: mesh entity ID */
		int Ent_nrdof,     /* in: number of dofs associated with the entity */
		int Vect_id,       /* in: vector ID in case of multiple solution vectors */
		double * Vect_dofs /* in: dofs to be written */
)
{
	/* pointer to field structure */
	apt_field * field_p;

	/* auxiliary variables */
	double * glob_dofs;
	int i;//,j

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	assert(Vect_dofs[0] > -10e200 && Vect_dofs[0] < 10e200);

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

#ifdef DEBUG_APM
	{
		/* check input */
		if (field_p->dof_ents[Ent_id].vec_dof_1 == NULL) {
			printf("Wrong dof entity in write_ent_dofs !\n");
			exit(-1);
		}
	}
#endif

//printf("Vect_id %d Ent_nrdof %d ",Vect_id,Ent_nrdof);

	if (Vect_id <= 1) {
		glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
	}
	else if (Vect_id == 2) {
		glob_dofs = field_p->dof_ents[Ent_id].vec_dof_2;
	}
	else if (Vect_id == 3) {
		glob_dofs = field_p->dof_ents[Ent_id].vec_dof_3;
	}
	else {
		field_p->dof_ents[Ent_id].vec_dof_1;
	}

	for (i = 0; i < Ent_nrdof; i++) glob_dofs[i] = Vect_dofs[i];

	return (1);
}


/*------------------------------------------------------------
  apr_create_ent_dofs - to create a vector of dofs associated with a given
				   mesh entity within approximation field data structure
------------------------------------------------------------*/
int apr_create_ent_dofs(/* returns: >=0 - success code, <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Ent_type,      /* in: type of mesh entity */
		int Ent_id,        /* in: mesh entity ID */
		int Ent_nrdof,     /* in: number of dofs associated with the entity */
		int Vect_id,       /* in: vector ID in case of multiple solution vectors */
		double * Vect_dofs /* in: dofs to be written */
)
{
	/* pointer to field structure */
	apt_field * field_p;

	/* auxiliary variables */
	double * glob_dofs;
//  int i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

#ifdef DEBUG_APM
	{
		/* check input */
		if (Ent_nrdof != apr_get_ent_nrdofs(Field_id, Ent_type, Ent_id)) {
			printf("Wrong number of dofs in write_ent_dofs !\n");
			exit(-1);
		}
	}
#endif

	if (Vect_id <= 1) {
		field_p->dof_ents[Ent_id].vec_dof_1 =
				(double *) malloc(Ent_nrdof * sizeof(double));
		glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
	}
	else if (Vect_id == 2) {
		field_p->dof_ents[Ent_id].vec_dof_2 =
				(double *) malloc(Ent_nrdof * sizeof(double));
		glob_dofs = field_p->dof_ents[Ent_id].vec_dof_2;
	}
	else if (Vect_id == 3) {
		field_p->dof_ents[Ent_id].vec_dof_3 =
				(double *) malloc(Ent_nrdof * sizeof(double));
		glob_dofs = field_p->dof_ents[Ent_id].vec_dof_3;
	}
	else  {
		field_p->dof_ents[Ent_id].vec_dof_1 =
				(double *) malloc(Ent_nrdof * sizeof(double));
		glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
	}

	return (1);
}

/*------------------------------------------------------------
  apr_set_ini_con - to set an initial condition
------------------------------------------------------------*/
int apr_set_ini_con(/* returns: >=0 - success code, <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		double (*Fun_p)(int, double *, int) /* pointer to function that provides */
		/* problem dependent initial condition data */
)
{
	double x[3];
	int sol_comp;
	double Uval[APC_MAXEQ];

	int nreq = apr_get_nreq(Field_id);
	int mesh_id = apr_get_mesh_id(Field_id);

	// in a loop over all vertices (nodes)
	int node_id = 0;
	while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
		if (apr_get_ent_pdeg(Field_id, APC_VERTEX, node_id) > 0) {

			// get node coordinates
			mmr_node_coor(mesh_id, node_id, x);

			// for each field vector component
			for (sol_comp = 0; sol_comp < nreq; sol_comp++) {

				// get initial condition value
				Uval[sol_comp] = (*Fun_p)(Field_id, x, sol_comp);

				/*kbw
					printf("in set_ini_con - get value %lf from pd function\n",
						   (*Fun_p)(Field_id, x, sol_comp));
				/*kew*/

				// set value as initial condition for a given node and solution component
				apr_write_ent_dofs(Field_id, APC_VERTEX, node_id, nreq, 0, Uval);

			} // end loop over solution components
		} // end if active node
	} // end loop over vertices (nodes)

	return (1);
}


/*------------------------------------------------------------
  apr_prepare_integration_parameters - used e.g. by apr_num_int_el
------------------------------------------------------------*/
int apr_prepare_integration_parameters(
		int Field_id,     /* in: approximation field ID  */
		int El_id,        /* in: element ID */
		int * Geo_order,  /* out: geometrical order of approximation */
		int * Num_geo_dofs, /* out: number of geometrical degrees of freedom */
		double * Geo_dofs, /* out: geometrical degrees of freedom */
		int * El_mate,    /* out: material index (in materials database) */
		int * Base,       /* out: type of basis functions */
		int * Pdeg_vec,   /* out: degree of approximation symbol or vector */
		int * Num_shap,   /* out: number of shape functions */
		int * Nreq,       /* out: number of components for unknowns */
		int * Num_dofs    /* out: number of scalar dofs (num_dofs = nreq*num_shap) */
)
{
	int mesh_id;
	int pdeg;		/* degree of polynomial */
	int el_nodes[MMC_MAXELVNO + 1];      /* list of nodes of El */

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	*Geo_order = 1; // geometrically (multi)linear elements
	/* get the coordinates of the nodes of El in the right order */
	mmr_el_node_coor(mesh_id, El_id, el_nodes, Geo_dofs);
	/* for geometrically (multi)linear elements number of degrees of freedom */
	/* is equal to the number of vertices - classical FEM nodes */
	*Num_geo_dofs = el_nodes[0];

#ifdef DEBUG_APM
	if (mmr_el_status(mesh_id, El_id) != MMC_ACTIVE) {
		printf("Asking for pdeg of inactive element in apr_num_int_el !\n");
		exit(-1);
	}
#endif

	*El_mate =  mmr_el_groupID(mesh_id, El_id);
	*Base = apr_get_base_type(Field_id, El_id);

	apr_get_el_pdeg(Field_id, El_id, &pdeg);
	Pdeg_vec[0] = pdeg;

	*Num_shap = apr_get_el_pdeg_numshap(Field_id, El_id, &pdeg);
	*Nreq = apr_get_nreq(Field_id);
	*Num_dofs = (*Nreq) * (*Num_shap);

	return (1);

}

/*------------------------------------------------------------
  apr_num_int_el - to perform numerical integration for an element
------------------------------------------------------------*/
int apr_num_int_el(
		int Problem_id,
		int Field_id,    /* in: approximation field ID  */
		int El_id,       /* in: unique identifier of the element */
		int Comp_sm,     /* in: indicator for the scope of computations: */
		/*   APC_NO_COMP  - do not compute anything */
		/*   APC_COMP_SM - compute entries to stiff matrix only */
		/*   APC_COMP_RHS - compute entries to rhs vector only */
		/*   APC_COMP_BOTH - compute entries for sm and rhsv */
		int * Pdeg_vec,       /* in: enforced degree of polynomial (if !=NULL ) */
		double * Sol_dofs_k,  /* in: solution dofs from previous iteration */
		/*     (for nonlinear problems) */
		double * Sol_dofs_n,  /* in: solution dofs from previous time step */
		/*     (for nonlinear problems) */
		/* out: various matrices and vectors containing integrals that appear in */
		/*      the weak form */
		/* 1. matrices for the left hand side with the size num_dofs*num_dofs */
		/*      (num_dofs=nreq*num_shap) */
		int * Diagonal,  /* array of indicators whether matrices are diagonal */
		double * Stiff_mat,	/* out: stiffness matrix stored columnwise */
		double * Rhs_vect	/* out: rhs vector */
		// REMARK:
		//   1. udofs, ueq - correspond to the solution u and go from left to right
		//   2. wdofs, weq - correspond to test functions w and go up and down
		//   3. matrices are stored columnwise in vectors
		//   4. solution indices change in rows (go from left to right)
		//   5. test functions indices change in columns (go up and down)
		//   6. when matrices are stored in a vector the index is computed as:
		//      vector[udofs*nreq*num_dofs+wdofs*nreq+ueq*num_dofs+weq]
		//      (num_dofs=num_shap*nreq)
		//   7. for each pair (wdofs,udofs) there is a small submatrix nreq x nreq
		//      with indices ueq (from left to right) and weq (from top to bottom)
		//   8. when stiffness matrix entries are computed for vector problems with
		//      the same shape functions for each component, then for each pair
		//      (wdofs,udofs) a matrix of coefficients with the size nreq x nreq
		//      is provided by the problem dependent module; its index is computed as
		//      vector[ueq*nreq+weq]
		//   9. when load vector entries are computed for vector problems with
		//      the same shape functions for each component, then for each wdofs index
		//      a vector of coefficients with the size nreq is provided by the problem
		//      dependent module
)
{

	/* pde coefficients */
	static int coeff_ind = 0;
	int coeff_vect_ind[23] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	// input to pdr_select_el_coeff_vect:
	// 0 - perform substitutions in coeff_vect_ind only if coeff_vect_ind[0]==0
	// output from pdr_select_el_coeff_vect:
	// 0 - coeff_vect_ind[0]==1 - always
	// 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz,
	// 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
	// 14 - tx, 15 - ty, 16 - tz, 17 - cval
	// 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval

	double axx[APC_MAXEQ * APC_MAXEQ];
	double axy[APC_MAXEQ * APC_MAXEQ];
	double axz[APC_MAXEQ * APC_MAXEQ];
	double ayx[APC_MAXEQ * APC_MAXEQ];
	double ayy[APC_MAXEQ * APC_MAXEQ];
	double ayz[APC_MAXEQ * APC_MAXEQ];
	double azx[APC_MAXEQ * APC_MAXEQ];
	double azy[APC_MAXEQ * APC_MAXEQ];
	double azz[APC_MAXEQ * APC_MAXEQ];
	double bx[APC_MAXEQ * APC_MAXEQ];
	double by[APC_MAXEQ * APC_MAXEQ];
	double bz[APC_MAXEQ * APC_MAXEQ];
	double tx[APC_MAXEQ * APC_MAXEQ];
	double ty[APC_MAXEQ * APC_MAXEQ];
	double tz[APC_MAXEQ * APC_MAXEQ];
	double cval[APC_MAXEQ * APC_MAXEQ];
	double mval[APC_MAXEQ * APC_MAXEQ];
	double qx[APC_MAXEQ];
	double qy[APC_MAXEQ];
	double qz[APC_MAXEQ];
	double sval[APC_MAXEQ];
	double lval[APC_MAXEQ];

	/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
	/* pde coefficients */
	/*   static int ndofs_max=0;            /\* local dimension of the problem *\/ */
	/*   static double *select_coeff=NULL; */
	/*   static double  *axx=NULL, *axy=NULL, *axz=NULL, *ayx=NULL, *ayy=NULL; */
	/*   static double  *ayz=NULL, *azx=NULL, *azy=NULL, *azz=NULL; */
	/*   static double  *bx=NULL, *by=NULL, *bz=NULL, *tx=NULL, *ty=NULL, *tz=NULL;  */
	/*   static double  *cval=NULL, *mval=NULL, *lval=NULL, *sval=NULL; */
	/*   static double  *qx=NULL, *qy=NULL, *qz=NULL; */

	/* #pragma omp threadprivate (ndofs_max) */
	/* #pragma omp threadprivate (select_coeff) */
	/* #pragma omp threadprivate (axx) */
	/* #pragma omp threadprivate (axy) */
	/* #pragma omp threadprivate (ayx) */
	/* #pragma omp threadprivate (axz) */
	/* #pragma omp threadprivate (ayz) */
	/* #pragma omp threadprivate (ayy) */
	/* #pragma omp threadprivate (ayz) */
	/* #pragma omp threadprivate (azy) */
	/* #pragma omp threadprivate (azz) */
	/* #pragma omp threadprivate (bx) */
	/* #pragma omp threadprivate (by) */
	/* #pragma omp threadprivate (bz) */
	/* #pragma omp threadprivate (tx) */
	/* #pragma omp threadprivate (ty) */
	/* #pragma omp threadprivate (tz) */
	/* #pragma omp threadprivate (cval) */
	/* #pragma omp threadprivate (mval) */
	/* #pragma omp threadprivate (lval) */
	/* #pragma omp threadprivate (sval) */
	/* #pragma omp threadprivate (qx) */
	/* #pragma omp threadprivate (qy) */
	/* #pragma omp threadprivate (qz) */


	/* quadrature rules */
	/*   static int problem_id_old=-1; /\* indicator for recomputing data *\/ */
	/*   int pdeg;		/\* degree of polynomial *\/ */
	/*   static int pdeg_old=-1; /\* indicator for recomputing quadrature data *\/ */
	/*   int base;		/\* type of basis functions *\/ */
	/*   static int base_old=-1; /\* indicator for recomputing quadrature data *\/ */
	/*   static int ngauss;            /\* number of gauss points *\/ */
	/*   static double xg[3000];   	 /\* coordinates of gauss points in 3D *\/ */
	/*   static double wg[1000];       /\* gauss weights *\/ */

	/* #pragma omp threadprivate (problem_id_old) */
	/* #pragma omp threadprivate (pdeg_old) */
	/* #pragma omp threadprivate (base_old) */
	/* #pragma omp threadprivate (ngauss) */
	/* #pragma omp threadprivate (xg) */
	/* #pragma omp threadprivate (wg) */

	// to make old OpenMP compilers working
	int problem_id_old = -1; /* indicator for recomputing data */
	int pdeg;		/* degree of polynomial */
	int pdeg_old = -1; /* indicator for recomputing quadrature data */
	int base;		/* type of basis functions */
	int base_old = -1; /* indicator for recomputing quadrature data */
	int ngauss;            /* number of gauss points */
	double xg[3000];   	 /* coordinates of gauss points in 3D */
	double wg[1000];       /* gauss weights */

	int geo_order, num_geo_dofs;
	double geo_dofs[3 * MMC_MAXELVNO]; /* coord of nodes of El */

	int nreq;	   /* number of equations */
	int el_mate;     /* material index in materials database */
	int num_shap;    /* number of element shape functions */
	int num_dofs;    /* number of element DOFs; usually num_dofs=nreq*num_shap */
	double determ;   /* determinant of jacobi matrix */
	double hsize = 0; /* size of an element */
	double vol;      /* volume for integration rule */
	double xcoor[3] = {0.0}; /* global coord of gauss point */
	double uk_val[APC_MAXEQ]; /* computed solution from previous iteration */
	double uk_x[APC_MAXEQ];   /* x-derivatives of components of uk_val */
	double uk_y[APC_MAXEQ];   /* y-derivatives of components of uk_val */
	double uk_z[APC_MAXEQ];   /* z-derivatives of components of uk_val */
	double un_val[APC_MAXEQ]; /* computed solution from previous time step */
	double un_x[APC_MAXEQ];   /* x-derivatives of components of un_val */
	double un_y[APC_MAXEQ];   /* y-derivatives of components of un_val */
	double un_z[APC_MAXEQ];   /* z-derivatives of components of un_val */
	double base_phi[APC_MAXELVD];    /* basis functions */
	double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
	double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
	double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */

	int i, ki, kk, udofs, wdofs, iaux, weq, ueq, j, k;


	/*++++++++++++++++ executable statements ++++++++++++++++*/

//#define COUNT_OPERATIONS
#ifdef COUNT_OPERATIONS
	int nr_oper = 0;
#endif

	/* prepare data necessary for numerical integration: */
	/* geo_dofs - geometry degrees of freedom - used by apr_elem_calc_3D  */
	/*            (for linear elements these are vertices coordinates) */
	/* pdeg - indicator of element degree of approximation (this may be a single */
	/*        number for simple (linear, quadratic, etc.) elements or a vector */
	/*        with separate degrees of approximation for each mesh entity */
	/*        (vertices, edges, faces, interior) forming an element */
	/* num_shap - number of shape functions */
	/* nreq - number of components for unknowns being vector fields */
	/* num_dofs - number of degrees of freedom (often num_dofs = nreq*num_shap) */
	/* base - type of basis functions (differentiates element types as well)*/
	/*          examples (from include/aph_intf.h): */
	/*          #define APC_BASE_PRISM_STD  3   // for linear prismatic elements */
	/*          #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements */
	apr_prepare_integration_parameters(Field_id, El_id,
			&geo_order, &num_geo_dofs, geo_dofs,
			&el_mate, &base, &pdeg,
			&num_shap, &nreq, &num_dofs);


	/* degree of polynomial passed as argument overrides default element values */
	if (Pdeg_vec != NULL) {
		pdeg = Pdeg_vec[0]; // for simple one number Pdeg_vec
		num_shap = apr_get_el_pdeg_numshap(Field_id, El_id, &pdeg);
		nreq = apr_get_nreq(Field_id);
		num_dofs = nreq * num_shap;
	}


#ifdef DEBUG_APM
	if (nreq > APC_MAXEQ) {
		printf("Number of equations %d greater than the limit %d.\n", nreq, APC_MAXEQ);
		printf("Change APC_MAXEQ in aph_std_lin_prism.h and recompile the code.\n");
		exit(-1);
	}
	if (num_shap > APC_MAXELVD) {
		printf("Number of shape functions %d greater than the limit %d.\n",
				num_shap, APC_MAXELVD);
		printf("Change APC_MAXELVD in include/aph_intf.h and recompile the code.\n");
		exit(-1);
	}
	if (num_dofs > APC_MAXELSD) {
		printf("Number of element degrees of freedom %d greater than the limit %d.\n"
				, num_dofs, APC_MAXELSD);
		printf("Change APC_MAXELSD in include/aph_intf.h and recompile the code.\n");
		exit(-1);
	}
#endif


	// if selection concerned a different problem or has not been done yet
	if (Problem_id != problem_id_old || coeff_ind == 0) {
		pdr_select_el_coeff_vect(Problem_id, coeff_vect_ind);
		// input to pdr_select_el_coeff_vect:
		// 0 - perform substitutions in coeff_vect_ind only if coeff_vect_ind[0]==0
		// output from pdr_select_el_coeff_vect:
		// 0 - coeff_vect_ind[0]==1 - always
		// 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz,
		// 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
		// 14 - tx, 15 - ty, 16 - tz, 17 - cval
		// 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval
	}



	/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
	/* if(Problem_id!=problem_id_old || select_coeff==NULL){ */
	/*   /\* allocate storage for coefficients and select the neeeded ones *\/ */
	/*   select_coeff=pdr_select_el_coeff(Problem_id, &mval, */
	/* 	 &axx,&axy,&axz,&ayx,&ayy,&ayz,&azx,&azy,&azz, */
	/* 	 &bx,&by,&bz,&tx,&ty,&tz,&cval,&lval,&qx,&qy,&qz,&sval); */
	/* } */

	/* prepare data for gaussian integration */
	if (pdeg != pdeg_old || base != base_old) {
		apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
		pdeg_old = pdeg;
		base_old = base;
	}

	/*kbw
	  if(Field_id==2 && El_id==13753)
		{
		printf("In num_int_el: Field_id %d, element %d\n",
		   Field_id, El_id);
		printf("pdeg %d, ngauss %d\n",pdeg, ngauss);
		printf("NREQ %d, ndof %d, local_dim %d\n",nreq, num_shap, num_dofs);
		printf("%d geo_dofs (nodes) with coordinates:\n", num_geo_dofs);
		for(i=0;i<num_geo_dofs;i++){
		  printf("geo_dofs (node) %d: x - %f, y - %f, z - %f\n",
			 i,
			 geo_dofs[3*i],geo_dofs[3*i+1],geo_dofs[3*i+2]);
		}
		printf("DOFS k:\n");
		for(i=0;i<num_dofs;i++) printf("%20.15lf",Sol_dofs_k[i]);
		printf("\n");
		printf("DOFS n:\n");
		for(i=0;i<num_dofs;i++) printf("%20.15lf",Sol_dofs_n[i]);
		printf("\n");
		//getchar();
	  }
	/*kew*/

//#define TIME_TEST_2
#ifdef TIME_TEST_2
	double t00 = 0.0, t01 = 0.0, t02 = 0.0, t03 = 0.0, t04 = 0.0, t05 = 0.0;
#endif

	/*! ----------------------------------------------------------------------! */
	/*! --------------------- LOOP OVER GAUSS POINTS -------------------------! */
	/*! ----------------------------------------------------------------------! */
	for (ki = 0; ki < ngauss; ki++) {

#ifdef TIME_TEST_2
		t00 = time_clock();
#endif

		/* at the gauss point, compute basis functions, determinant etc*/
		iaux = 2; /* calculations with jacobian but not on the boundary */
		determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base,
						&xg[3 * ki], geo_dofs, Sol_dofs_k,
						base_phi, base_dphix, base_dphiy, base_dphiz,
						xcoor, uk_val, uk_x, uk_y, uk_z, NULL);

		vol = determ * wg[ki];

#ifdef TIME_TEST_2
		t01 += time_clock() - t00;
#endif

#ifdef COUNT_OPERATIONS
		// for DG
		//nr_oper += ((4+15)*num_shap + 3*(pdeg+1)*(pdeg+2)/2 + 6*(pdeg+1) + 220);
		// for STD - circa 500; definitely < 1000
		nr_oper += 21 + 24 + 15 * (6 + num_shap) + 12 * 6 + 21 + 8 * num_shap * nreq;
#endif


		/*kbw
			//if(norm_u>0.00010){
			  if(El_id==49031){
				int idofs;
				double hsize = mmr_el_hsize(apr_get_mesh_id(Field_id),
							El_id,NULL,NULL,NULL);
				printf("element %d, num_shap %d, u_x %lf, u_y %lf, u_z %lf\n",
				   El_id, num_shap, uk_val[0], uk_val[1], uk_val[2]);
				for (idofs = 0; idofs < num_shap; idofs++) {
				  printf("size: vol %.12lf, standard %.12lf, dphix %lf, dphiy %lf, dphiz %lf\n",
					 vol, hsize,
					 base_dphix[idofs], base_dphiy[idofs], base_dphiz[idofs]);
			  }
				getchar();  getchar();  getchar();
			  }
		/*kew*/

		// if necessary compute the values of solution at the previous time step
		if (Sol_dofs_n != NULL) {

			int ieq;
			for (ieq = 0; ieq < nreq; ieq++) un_val[ieq] = 0.0;
			for (i = 0; i < num_shap; i++) {
				for (ieq = 0; ieq < nreq; ieq++) {
					un_val[ieq] += Sol_dofs_n[i * nreq + ieq] * base_phi[i];
				}
			}
			for (ieq = 0; ieq < nreq; ieq++) {
				un_x[ieq] = 0.0;
				un_y[ieq] = 0.0;
				un_z[ieq] = 0.0;
			}
			for (i = 0; i < num_shap; i++) {
				for (ieq = 0; ieq < nreq; ieq++) {
					un_x[ieq] += Sol_dofs_n[i * nreq + ieq] * base_dphix[i];
					un_y[ieq] += Sol_dofs_n[i * nreq + ieq] * base_dphiy[i];
					un_z[ieq] += Sol_dofs_n[i * nreq + ieq] * base_dphiz[i];
				}
			}

#ifdef COUNT_OPERATIONS
			// < 333 for std
			nr_oper += nreq * (num_shap * 8 + 4);
			printf("\nNumber of operations in apr_num_int_el for single element %d - before\n", nr_oper);
#endif

		}

		// REMARK:
		//   1. udofs, ueq - correspond to the solution u and go from left to right
		//   2. wdofs, weq - correspond to test functions w and go up and down
		//   3. matrices are stored columnwise in vectors
		//   4. solution indices change in rows (go from left to right)
		//   5. test functions indices change in columns (go up and down)
		//   6. when matrices are stored in a vector the index is computed as:
		//      vector[udofs*nreq*num_dofs+wdofs*nreq+ueq*num_dofs+weq]
		//      (num_dofs=num_shap*nreq)
		//   7. for each pair (wdofs,udofs) there is a small submatrix nreq x nreq
		//      with indices ueq (from left to right) and weq (from top to bottom)
		//   8. when stiffness matrix entries are computed for vector problems with
		//      the same shape functions for each component, then for each pair
		//      (wdofs,udofs) a matrix of coefficients with the size nreq x nreq
		//      is provided by the problem dependent module - the matrix is stored
		//      columnwise in the vector and the index is computed as [ueq*nreq+weq]
		//   9. when load vector entries are computed for vector problems with
		//      the same shape functions for each component, then for each wdofs index
		//      a vector of coefficients with the size nreq is provided by the problem
		//      dependent module

#ifdef TIME_TEST_2
		t02 += time_clock() - t00;
#endif

		/* get coefficients of convection-diffusion-reaction equations */
		pdr_el_coeff(Problem_id, El_id, el_mate, hsize, pdeg, &xg[3 * ki],
				base_phi, base_dphix, base_dphiy, base_dphiz,
				xcoor, uk_val, uk_x, uk_y, uk_z, un_val, un_x, un_y, un_z,
				mval, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz,
				bx, by, bz, tx, ty, tz, cval, lval, qx, qy, qz, sval);

#ifdef TIME_TEST_2
		t03 += time_clock() - t00;
#endif

		/*kbw
		#pragma omp critical(printing_123)
			  {
			  //static int juz_bylo=0;
			  //if(juz_bylo==0){

			  //juz_bylo = 1;
			  //if(Field_id==1 && El_id==123)
			  //if(fabs(xcoor[1]-0.55) < 0.1)
			  extern double pdr_time_d_params(int Problem_id, int Num);
			  double cur_time = pdr_time_d_params(Problem_id,4);
			  if(cur_time > 0.03)
			{

			int i,ieq;

			printf("\n\n\n########################################\n");
			printf("In num_int_el (time %lf): Field_id %d, element %d\n", cur_time, Field_id, El_id);
			printf("\npdeg %d, ngauss %d\n",pdeg, ngauss);
			printf("NREQ %d, ndof %d, local_dim %d\n",nreq, num_shap, num_dofs);
			printf("\n%d geo_dofs (nodes) with coordinates:\n", num_geo_dofs);
			for(i=0;i<num_geo_dofs;i++){
			  printf("geo_dofs (node) %d: x - %f, y - %f, z - %f\n",
				 i,
				 geo_dofs[3*i],geo_dofs[3*i+1],geo_dofs[3*i+2]);
			}

			  printf("\n============================================================\n");
			  printf("at gauss point %d, local coor %lf, %lf, %lf", ki,xg[3*ki],xg[3*ki+1],xg[3*ki+2]);
			  printf("\n============================================================\n");
			  printf("global coor %lf %lf %lf\n",xcoor[0],xcoor[1],xcoor[2]);
			  printf("weight %lf, determ %lf, vol %lf\n",
				 wg[ki],determ,vol);
			  printf("\n%d shape functions and derivatives: \n", num_shap);
			  for(i=0;i<num_shap;i++){
			printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
				   base_phi[i],base_dphix[i],base_dphiy[i],base_dphiz[i]);
			  }

			  printf("\nprevious iteration solution and derivatives: \n");
			  for(ieq=0;ieq<nreq;ieq++){
			printf("component %d: ukx - %lf, der: x - %lf, y - %lf, z - %lf\n",
				   ieq, uk_val[ieq],uk_x[ieq],uk_y[ieq],uk_z[ieq]);
			  }

			  printf("\nprevious timestep solution and derivatives: \n");
			  for(ieq=0;ieq<nreq;ieq++){
			  printf("component %d: unx - %lf, der: x - %lf, y - %lf, z - %lf\n",
				 ieq, un_val[ieq],un_x[ieq],un_y[ieq],un_z[ieq]);
			  }

			  if(mval!=NULL){
			printf("\ntime coeff LHS:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf",mval[ueq]);}
			  }
			  if(axx!=NULL){
			printf("\n\ndiffusion coeff axx:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",axx[ueq]);}
			printf("\n\ndiffusion coeff axy:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",axy[ueq]);}
			printf("\n\ndiffusion coeff axz:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",axz[ueq]);}
			  }
			  if(ayy!=NULL){
			printf("\n\ndiffusion coeff ayx:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",ayx[ueq]);}
			printf("\n\ndiffusion coeff ayy:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",ayy[ueq]);}
			printf("\n\ndiffusion coeff ayz:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",ayz[ueq]);}
			  }
			  if(azz!=NULL){
			printf("\n\ndiffusion coeff azx:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",azx[ueq]);}
			printf("\n\ndiffusion coeff azy:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",azy[ueq]);}
			printf("\n\ndiffusion coeff azz:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",azz[ueq]);}
			  }
			  if(bx!=NULL){
			printf("\n\nconvection coeff bx:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",bx[ueq]);}
			printf("\n\nconvection coeff by:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",by[ueq]);}
			printf("\n\nconvection coeff bz:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",bz[ueq]);}
			  }
			  if(tx!=NULL){
			printf("\n\nconvection coeff tx:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",tx[ueq]);}
			printf("\n\nconvection coeff ty:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",ty[ueq]);}
			printf("\n\nconvection coeff tz:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",tz[ueq]);}
			  }
			  if(cval!=NULL){
			printf("\n\nreaction coeff c:\n");
			for(ueq=0;ueq<nreq*nreq;ueq++) { ueq%4?printf("\t"):printf("\n"); printf("%20.12lf\t",cval[ueq]);}
			  }
			  if(lval!=NULL){
			printf("\n\ntime coeff RHS:\n");
			for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf\t",lval[ueq]);
			  }
			  if(qx!=NULL){
			printf("\n\nqx coeff RHS:\n");
			for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf\t",qx[ueq]);
			printf("\n\nqy coeff RHS:\n");
			for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf\t",qy[ueq]);
			printf("\n\nqz coeff RHS:\n");
			for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf\t",qz[ueq]);
			  }
			  if(sval!=NULL){
			printf("\n\nsource:\n");
			for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf\t",sval[ueq]);
			  }
			printf("\n\ngetchar() -> ");
			getchar();
			}
			  }
		/*kew*/

		if (Comp_sm == APC_COMP_SM || Comp_sm == APC_COMP_BOTH) {

			if (coeff_vect_ind[1] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							Stiff_mat[kk + wdofs] +=
									mval[0] * base_phi[udofs] * base_phi[wdofs] * vol;

						}/* wdofs */
						kk += num_shap;
					} /* udofs */

				}
				else if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[0] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] +=
												mval[ueq * nreq + weq] * base_phi[udofs] * base_phi[wdofs] * vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] +=
											mval[weq * nreq + weq] * base_phi[udofs] * base_phi[wdofs] * vol;

								}

							}

						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}

#ifdef COUNT_OPERATIONS
				if (Diagonal[0] == 0) nr_oper += (2 * nreq * nreq + 2) * num_shap * num_shap;
				else nr_oper += (2 * nreq + 2) * num_shap * num_shap;
				printf("\nNumber of operations in apr_num_int_el for single element %d - m\n", nr_oper);
#endif

			}

			if (coeff_vect_ind[10] == 1) {
				if (coeff_vect_ind[2] == 1) {
					if (nreq == 1) {

						kk = 0;
						for (udofs = 0; udofs < num_shap; udofs++) {
							for (wdofs = 0; wdofs < num_shap; wdofs++) {

								Stiff_mat[kk + wdofs] += (
												axx[0] * base_dphix[udofs] * base_dphix[wdofs] +
												axy[0] * base_dphiy[udofs] * base_dphix[wdofs] +
												axz[0] * base_dphiz[udofs] * base_dphix[wdofs] +
												ayx[0] * base_dphix[udofs] * base_dphiy[wdofs] +
												ayy[0] * base_dphiy[udofs] * base_dphiy[wdofs] +
												ayz[0] * base_dphiz[udofs] * base_dphiy[wdofs] +
												azx[0] * base_dphix[udofs] * base_dphiz[wdofs] +
												azy[0] * base_dphiy[udofs] * base_dphiz[wdofs] +
												azz[0] * base_dphiz[udofs] * base_dphiz[wdofs]
										) * vol;

							}/* wdofs */
							kk += num_shap;
						} /* udofs */
					}
					else if (nreq > 1) {

						kk = 0;
						for (udofs = 0; udofs < num_shap; udofs++) {
							for (wdofs = 0; wdofs < num_shap; wdofs++) {

								if (Diagonal[1] == 0) {

									for (ueq = 0; ueq < nreq; ueq++) {
										for (weq = 0; weq < nreq; weq++) {


											Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += (
															// axx[ueq*nreq+weq] * dw[weq]/dx * du[ueq]/dx
															axx[ueq * nreq + weq] * base_dphix[udofs] * base_dphix[wdofs] +
															// axy[ueq*nreq+weq] * dw[weq]/dx * du[ueq]/dy
															axy[ueq * nreq + weq] * base_dphiy[udofs] * base_dphix[wdofs] +
															// axz[ueq*nreq+weq] * dw[weq]/dx * du[ueq]/dz
															axz[ueq * nreq + weq] * base_dphiz[udofs] * base_dphix[wdofs] +
															// ayx[ueq*nreq+weq] * dw[weq]/dy * du[ueq]/dx
															ayx[ueq * nreq + weq] * base_dphix[udofs] * base_dphiy[wdofs] +
															// ayy[ueq*nreq+weq] * dw[weq]/dy * du[ueq]/dy
															ayy[ueq * nreq + weq] * base_dphiy[udofs] * base_dphiy[wdofs] +
															// ayz[ueq*nreq+weq] * dw[weq]/dy * du[ueq]/dz
															ayz[ueq * nreq + weq] * base_dphiz[udofs] * base_dphiy[wdofs] +
															// azx[ueq*nreq+weq] * dw[weq]/dz * du[ueq]/dx
															azx[ueq * nreq + weq] * base_dphix[udofs] * base_dphiz[wdofs] +
															// azy[ueq*nreq+weq] * dw[weq]/dz * du[ueq]/dy
															azy[ueq * nreq + weq] * base_dphiy[udofs] * base_dphiz[wdofs] +
															// azz[ueq*nreq+weq] * dw[weq]/dz * du[ueq]/dz
															azz[ueq * nreq + weq] * base_dphiz[udofs] * base_dphiz[wdofs]
													) * vol;

										}
									}

								} else {

									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += (
														axx[weq * nreq + weq] * base_dphix[udofs] * base_dphix[wdofs] +
														axy[weq * nreq + weq] * base_dphiy[udofs] * base_dphix[wdofs] +
														axz[weq * nreq + weq] * base_dphiz[udofs] * base_dphix[wdofs] +
														ayx[weq * nreq + weq] * base_dphix[udofs] * base_dphiy[wdofs] +
														ayy[weq * nreq + weq] * base_dphiy[udofs] * base_dphiy[wdofs] +
														ayz[weq * nreq + weq] * base_dphiz[udofs] * base_dphiy[wdofs] +
														azx[weq * nreq + weq] * base_dphix[udofs] * base_dphiz[wdofs] +
														azy[weq * nreq + weq] * base_dphiy[udofs] * base_dphiz[wdofs] +
														azz[weq * nreq + weq] * base_dphiz[udofs] * base_dphiz[wdofs]
												) * vol;

									}
								}

							} /* wdofs */
							kk += nreq * num_dofs;

						} /* udofs */

					}

#ifdef COUNT_OPERATIONS
					if (nreq == 1) nr_oper += 28 * num_shap * num_shap;
					else {
						if (Diagonal[1] == 0) nr_oper += 9 * (2 * nreq * nreq + 2) * num_shap * num_shap;
						else nr_oper += 9 * (2 * nreq + 2) * num_shap * num_shap;
					}
					printf("\nNumber of operations in apr_num_int_el for single element %d - a\n", nr_oper);
#endif

				} else { // if axy == NULL
					if (nreq == 1) {

						kk = 0;
						for (udofs = 0; udofs < num_shap; udofs++) {
							for (wdofs = 0; wdofs < num_shap; wdofs++) {

								Stiff_mat[kk + wdofs] += (
												axx[0] * base_dphix[udofs] * base_dphix[wdofs] +
												ayy[0] * base_dphiy[udofs] * base_dphiy[wdofs] +
												azz[0] * base_dphiz[udofs] * base_dphiz[wdofs]
										) * vol;

							}/* wdofs */
							kk += num_shap;
						} /* udofs */
					}
					else if (nreq > 1) {

						kk = 0;
						for (udofs = 0; udofs < num_shap; udofs++) {
							for (wdofs = 0; wdofs < num_shap; wdofs++) {

								if (Diagonal[1] == 0) {

									for (ueq = 0; ueq < nreq; ueq++) {
										for (weq = 0; weq < nreq; weq++) {

											Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += (
															axx[ueq * nreq + weq] * base_dphix[udofs] * base_dphix[wdofs] +
															ayy[ueq * nreq + weq] * base_dphiy[udofs] * base_dphiy[wdofs] +
															azz[ueq * nreq + weq] * base_dphiz[udofs] * base_dphiz[wdofs]
													) * vol;

										}
									}

								} else {

									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += (
														axx[weq * nreq + weq] * base_dphix[udofs] * base_dphix[wdofs] +
														ayy[weq * nreq + weq] * base_dphiy[udofs] * base_dphiy[wdofs] +
														azz[weq * nreq + weq] * base_dphiz[udofs] * base_dphiz[wdofs]
												) * vol;

									}
								}

							} /* wdofs */
							kk += nreq * num_dofs;

						} /* udofs */

					}

#ifdef COUNT_OPERATIONS
					if (nreq == 1) nr_oper += 10 * num_shap * num_shap;
					else {
						if (Diagonal[1] == 0) nr_oper += 3 * (2 * nreq * nreq + 2) * num_shap * num_shap;
						else nr_oper += 3 * (2 * nreq + 2) * num_shap * num_shap;
					}
#endif

				} // endif if azz!= NULL and axy == NULL

			}
			else if (coeff_vect_ind[1] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							Stiff_mat[kk + wdofs] += (
											axx[0] * base_dphix[udofs] * base_dphix[wdofs] +
											axy[0] * base_dphiy[udofs] * base_dphix[wdofs] +
											ayx[0] * base_dphix[udofs] * base_dphiy[wdofs] +
											ayy[0] * base_dphiy[udofs] * base_dphiy[wdofs]
									) * vol;

						}/* wdofs */
						kk += num_shap;
					} /* udofs */
				}
				else if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[1] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += (
														axx[ueq * nreq + weq] * base_dphix[udofs] * base_dphix[wdofs] +
														axy[ueq * nreq + weq] * base_dphiy[udofs] * base_dphix[wdofs] +
														ayx[ueq * nreq + weq] * base_dphix[udofs] * base_dphiy[wdofs] +
														ayy[ueq * nreq + weq] * base_dphiy[udofs] * base_dphiy[wdofs]
												) * vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += (
													axx[weq * nreq + weq] * base_dphix[udofs] * base_dphix[wdofs] +
													axy[weq * nreq + weq] * base_dphiy[udofs] * base_dphix[wdofs] +
													ayx[weq * nreq + weq] * base_dphix[udofs] * base_dphiy[wdofs] +
													ayy[weq * nreq + weq] * base_dphiy[udofs] * base_dphiy[wdofs]
											) * vol;

								}
							}
						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 13 * num_shap * num_shap;
				else {
					if (Diagonal[1] == 0) nr_oper += 4 * (2 * nreq * nreq + 2) * num_shap * num_shap;
					else nr_oper += 4 * (2 * nreq + 2) * num_shap * num_shap;
				}
#endif

			}

			if (coeff_vect_ind[13] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							Stiff_mat[kk + wdofs] += ( bx[0] * base_dphix[udofs] * base_phi[wdofs] +
											by[0] * base_dphiy[udofs] * base_phi[wdofs] +
											bz[0] * base_dphiz[udofs] * base_phi[wdofs]
									) * vol;

						}/* wdofs */
						kk += num_shap;
					} /* udofs */
				}
				else if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[2] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += (
														//  bx[ueq*nreq+weq] * w[weq] * du[ueq]/dx
														bx[ueq * nreq + weq] * base_dphix[udofs] * base_phi[wdofs] +
														//  by[ueq*nreq+weq] * w[weq] * du[ueq]/dy
														by[ueq * nreq + weq] * base_dphiy[udofs] * base_phi[wdofs] +
														//  bz[ueq*nreq+weq] * w[weq] * du[ueq]/dz
														bz[ueq * nreq + weq] * base_dphiz[udofs] * base_phi[wdofs]
												) * vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += (
													bx[weq * nreq + weq] * base_dphix[udofs] * base_phi[wdofs] +
													by[weq * nreq + weq] * base_dphiy[udofs] * base_phi[wdofs] +
													bz[weq * nreq + weq] * base_dphiz[udofs] * base_phi[wdofs]
											) * vol;

								}
							}
						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 10 * num_shap * num_shap;
				else {
					if (Diagonal[2] == 0) nr_oper += 3 * (2 * nreq * nreq + 2) * num_shap * num_shap;
					else nr_oper += 3 * (2 * nreq + 2) * num_shap * num_shap;
				}
				printf("\nNumber of operations in apr_num_int_el for single element %d - b\n", nr_oper);
#endif

			}
			else if (coeff_vect_ind[11] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							Stiff_mat[kk + wdofs] += ( bx[0] * base_dphix[udofs] * base_phi[wdofs] +
											by[0] * base_dphiy[udofs] * base_phi[wdofs]
									) * vol;

						}/* wdofs */
						kk += num_shap;
					} /* udofs */
				}
				else if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[2] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += (
														bx[ueq * nreq + weq] * base_dphix[udofs] * base_phi[wdofs] +
														by[ueq * nreq + weq] * base_dphiy[udofs] * base_phi[wdofs]
												) * vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += (
													bx[weq * nreq + weq] * base_dphix[udofs] * base_phi[wdofs] +
													by[weq * nreq + weq] * base_dphiy[udofs] * base_phi[wdofs]
											) * vol;

								}
							}
						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 7 * num_shap * num_shap;
				else {
					if (Diagonal[2] == 0) nr_oper += 2 * (2 * nreq * nreq + 2) * num_shap * num_shap;
					else nr_oper += 2 * (2 * nreq + 2) * num_shap * num_shap;
				}
#endif

			}

			if (coeff_vect_ind[16] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							Stiff_mat[kk + wdofs] += (
											tx[0] * base_phi[udofs] * base_dphix[wdofs] +
											ty[0] * base_phi[udofs] * base_dphiy[wdofs] +
											tz[0] * base_phi[udofs] * base_dphiz[wdofs]
									)	* vol;

						}/* wdofs */
						kk += num_shap;
					} /* udofs */
				}
				else if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[3] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += (
														// tx[ueq*nreq+weq] * dw[weq]/dx * u[ueq]
														tx[ueq * nreq + weq] * base_phi[udofs] * base_dphix[wdofs] +
														// ty[ueq*nreq+weq] * dw[weq]/dz * u[ueq]
														ty[ueq * nreq + weq] * base_phi[udofs] * base_dphiy[wdofs] +
														// tz[ueq*nreq+weq] * dw[weq]/dz * u[ueq]
														tz[ueq * nreq + weq] * base_phi[udofs] * base_dphiz[wdofs]
												)	* vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += (
													tx[weq * nreq + weq] * base_phi[udofs] * base_dphix[wdofs] +
													ty[weq * nreq + weq] * base_phi[udofs] * base_dphiy[wdofs] +
													tz[weq * nreq + weq] * base_phi[udofs] * base_dphiz[wdofs]
											)	* vol;

								}
							}
						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 10 * num_shap * num_shap;
				else {
					if (Diagonal[3] == 0) nr_oper += 3 * (2 * nreq * nreq + 2) * num_shap * num_shap;
					else nr_oper += 3 * (2 * nreq + 2) * num_shap * num_shap;
				}
				printf("\nNumber of operations in apr_num_int_el for single element %d - t\n", nr_oper);
#endif

			}
			else if (coeff_vect_ind[14] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							Stiff_mat[kk + wdofs] += ( tx[0] * base_phi[udofs] * base_dphix[wdofs] +
											ty[0] * base_phi[udofs] * base_dphiy[wdofs]
									) * vol;

						}/* wdofs */
						kk += num_shap;
					} /* udofs */
				}
				else if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[3] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += (
														tx[ueq * nreq + weq] * base_phi[udofs] * base_dphix[wdofs] +
														ty[ueq * nreq + weq] * base_phi[udofs] * base_dphiy[wdofs]
												)	* vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += (
													tx[weq * nreq + weq] * base_phi[udofs] * base_dphix[wdofs] +
													ty[weq * nreq + weq] * base_phi[udofs] * base_dphiy[wdofs]
											)	* vol;

								}
							}
						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 7 * num_shap * num_shap;
				else {
					if (Diagonal[3] == 0) nr_oper += 2 * (2 * nreq * nreq + 2) * num_shap * num_shap;
					else nr_oper += 2 * (2 * nreq + 2) * num_shap * num_shap;
				}
#endif

			}


			if (coeff_vect_ind[17] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							Stiff_mat[kk + wdofs] += cval[0] * base_phi[udofs] * base_phi[wdofs] * vol;

						}/* wdofs */
						kk += num_shap;
					} /* udofs */
				}
				else if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[4] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] +=
												cval[ueq * nreq + weq] * base_phi[udofs] * base_phi[wdofs] * vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] +=
											cval[weq * nreq + weq] * base_phi[udofs] * base_phi[wdofs] * vol;

								}
							}
						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 4 * num_shap * num_shap;
				else {
					if (Diagonal[4] == 0) nr_oper += (2 * nreq * nreq + 2) * num_shap * num_shap;
					else nr_oper += (2 * nreq + 2) * num_shap * num_shap;
				}
				printf("\nNumber of operations in apr_num_int_el for single element %d - c\n", nr_oper);
#endif

			}

		} /* end if computing SM */

		if (Comp_sm == APC_COMP_RHS || Comp_sm == APC_COMP_BOTH) {

			if (coeff_vect_ind[18] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						Rhs_vect[kk] += lval[0] * base_phi[wdofs] * vol;

						kk++;
					}/* wdofs */
				}
				else if (nreq > 1) {

					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						for (weq = 0; weq < nreq; weq++) {

							Rhs_vect[wdofs * nreq + weq] += lval[weq] * base_phi[wdofs] * vol;


						}

					} /* wdofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 3 * num_shap;
				else nr_oper += (2 * nreq + 1) * num_shap;
				printf("\nNumber of operations in apr_num_int_el for single element %d - l\n", nr_oper);
#endif

			}


			if (coeff_vect_ind[21] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						Rhs_vect[kk] += ( qx[0] * base_dphix[wdofs] +
										qy[0] * base_dphiy[wdofs] +
										qz[0] * base_dphiz[wdofs] ) * vol;

						kk++;
					}/* wdofs */
				}
				else if (nreq > 1) {

					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						for (weq = 0; weq < nreq; weq++) {

							Rhs_vect[wdofs * nreq + weq] += (
											qx[weq] * base_dphix[wdofs] +
											qy[weq] * base_dphiy[wdofs] +
											qz[weq] * base_dphiz[wdofs] ) * vol;

						}

					} /* wdofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 7 * num_shap;
				else nr_oper += 3 * (2 * nreq + 1) * num_shap;
				printf("\nNumber of operations in apr_num_int_el for single element %d - q\n", nr_oper);
#endif

			}
			else if (coeff_vect_ind[19] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						Rhs_vect[kk] += ( qx[0] * base_dphix[wdofs] +
										qy[0] * base_dphiy[wdofs] ) * vol;

						kk++;
					}/* wdofs */
				}
				else if (nreq > 1) {

					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						for (weq = 0; weq < nreq; weq++) {

							Rhs_vect[wdofs * nreq + weq] += (
											qx[weq] * base_dphix[wdofs] +
											qy[weq] * base_dphiy[wdofs] ) * vol;

						}

					} /* wdofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 5 * num_shap;
				else nr_oper += 2 * (2 * nreq + 1) * num_shap;
#endif

			}


			if (coeff_vect_ind[22] == 1) {

				if (nreq == 1) {

					kk = 0;
					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						Rhs_vect[kk] += sval[0] * base_phi[wdofs] * vol;

						kk++;
					}/* wdofs */
				}

				else if (nreq > 1) {

					for (wdofs = 0; wdofs < num_shap; wdofs++) {

						for (weq = 0; weq < nreq; weq++) {

							Rhs_vect[wdofs * nreq + weq] += sval[weq] * base_phi[wdofs] * vol;

						}

					} /* wdofs */

				}

#ifdef COUNT_OPERATIONS
				if (nreq == 1) nr_oper += 3 * num_shap;
				else nr_oper += (2 * nreq + 1) * num_shap;
				printf("\nNumber of operations in apr_num_int_el for single element %d - s\n", nr_oper);
#endif

			}

		} /* end if computing RHSV */

		if (Comp_sm == APC_COMP_MM) {
			if (coeff_vect_ind[1] == 1) {
				if (nreq > 1) {

					kk = 0;
					for (udofs = 0; udofs < num_shap; udofs++) {
						for (wdofs = 0; wdofs < num_shap; wdofs++) {

							if (Diagonal[0] == 0) {

								for (ueq = 0; ueq < nreq; ueq++) {
									for (weq = 0; weq < nreq; weq++) {

										Stiff_mat[kk + wdofs * nreq + ueq * num_dofs + weq] += base_phi[udofs] * base_phi[wdofs] * vol;

									}
								}

							} else {

								for (weq = 0; weq < nreq; weq++) {

									Stiff_mat[kk + wdofs * nreq + weq * num_dofs + weq] += base_phi[udofs] * base_phi[wdofs] * vol;

								}

							}

						} /* wdofs */
						kk += nreq * num_dofs;

					} /* udofs */

				}
			}
		}
#ifdef TIME_TEST_2
		t04 += time_clock() - t00;
#endif

	} /* end loop over integration points: ki */


#ifdef TIME_TEST_2
	printf("EXECUTION TIME: t01 %lf, t02 %lf, t03 %lf, t04 %lf, t05 %lf\n",
			t01, t02 - t01, t03 - t02, t04 - t03, t05 - t04);
#endif

#ifdef COUNT_OPERATIONS
	printf("\nNumber of operations in apr_num_int_el for single element %d\n", nr_oper);
#endif

	return (1);
}



/*------------------------------------------------------------
  apr_get_stiff_mat_data - to return data on dof entities for an element and
					  to compute or rewrite element's stiffness matrix and RHSV
------------------------------------------------------------*/
int apr_get_stiff_mat_data(
		int Field_id,   /* in: approximation field ID  */
		int El_id,      /* in: unique identifier of the element */
		int Comp_sm,    /* in: indicator for the scope of computations: */
		/*   APC_NO_COMP  - do not compute anything */
		/*   APC_COMP_SM - compute entries to stiff matrix only */
		/*   APC_COMP_RHS - compute entries to rhs vector only */
		/*   APC_COMP_BOTH - compute entries for sm and rhsv */
		/*   APC_REWR_SM - rewrite only entries to stiff matrix only */
		/*   APC_REWR_RHS - rewrite only entries to rhs vector only */
		/*   APC_REWR_BOTH - rewrite only entries for sm and rhsv */
		char  Transpose,/* in: perform transposition while rewriting */
		/*     'y' or 'Y' - yes, otherwise - no */
		int Pdeg_in,    /* in: enforced degree of polynomial (if > 0 ) */
		int Nreq_in,    /* in: enforced nreq (if > 0 ) */
		int * Nr_dof_ent, /* in: size of arrays, */
		/* out: no of filled entries, i.e. number of mesh entities*/
		/* with which dofs and stiffness matrix blocks are associated */
		int * List_dof_ent_type, /* out: list of no of dofs for 'dof' entity */
		int * List_dof_ent_id,  /* out: list of no of dofs for 'dof' entity */
		int * List_dof_ent_nrdofs, /* out: list of no of dofs for 'dof' entity */
		/* below: optional means when Stiff_mat and Rhs_vect are computed */
		int * Nrdofs_loc, /* in(optional): size of Stiff_mat and Rhs_vect */
		/* out(optional): actual number of dofs per integration entity*/
		/* for matrices passing NULL is possible if suitable Comp_sm is passed */
		double * Stiff_mat,     /* out(optional): stiffness matrix stored columnwise*/
		double * Rhs_vect       /* out(optional): rhs vector */
		// REMARK:
		//   1. udofs, ueq - correspond to the solution and go from left to right
		//   2. wdofs, weq - correspond to test functions and go up and down
		//   3. matrices are stored columnwise in vectors
		//   4. solution indices change in rows (go from left to right)
		//   5. test functions indices change in columns (go up and down)
		//   6. when matrices are stored in a vector the index is computed as:
		//      vector[udofs*nreq*num_dofs+wdofs*nreq+ueq*num_dofs+weq]
		//      (num_dofs=num_shap*nreq)
		//   7. for each pair (wdofs,udofs) there is a small submatrix nreq x nreq
		//      with indices ueq (from left to right) and weq (from top to bottom)
)
{

	int mesh_id, nreq,  ino, is_constrained; //idofent,
	int num_shap, num_dofs, constr[5], pdeg, el_nodes_constr[MMC_MAXELVNO + 1];
	double constr_val, constr_el_loc[MMC_MAXELVNO + 1][4];
	int i, j;
	int el_nodes[MMC_MAXELVNO + 1];      /* list of nodes of El */
	int num_shap_constr;  /* number of element shape functions with constraints */
	int num_dofs_constr;  /* number of element degrees of freedom */
	/* with constraints taken into account */

	/* in general APC_MAXELVD should be used instead of 6 */
	double stiff_loc[APC_MAXEQ * 6 * APC_MAXEQ * 6];
	double rhs_loc[APC_MAXEQ * 6];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* checking */
	if (Transpose == 'y' || Transpose == 'Y') {
		printf(
				"Row-major stiffnes matrices not yet supported in apr_get_stiff_mata_data!\n");
		exit(-1);
	}

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	if (Nreq_in > 0) nreq = Nreq_in;
	else {
		nreq = apr_get_nreq(Field_id);
	}

	/* find degree of polynomial and number of element scalar dofs */
	if (Pdeg_in > 0) pdeg = Pdeg_in;
	else {
		apr_get_el_pdeg(Field_id, El_id, &pdeg);
	}
	num_shap = apr_get_el_pdeg_numshap(Field_id, El_id, &pdeg);
	num_dofs = num_shap * nreq;

	mmr_el_node_coor(mesh_id, El_id, el_nodes, NULL);

	is_constrained = 0;
	num_shap_constr = 0;

	/*kbw
	if(El_id==2232){
	  printf("In apr_get_stiff_mat_data: field %d, el %d, comp_sm %d, Pdeg_in %d\n",
		   Field_id, El_id, Comp_sm, Pdeg_in);
	  printf("pdeg %d, num_shap %d, num_nodes %d, el_nodes:",
		 pdeg, num_shap, el_nodes[0]);
	  for(i=1; i<=el_nodes[0]; i++) printf("  %d,",el_nodes[i]);
	  printf("\n");
	}
	/*kew*/

	for (ino = 0; ino < num_shap; ino++)
	{
		if (apr_get_ent_pdeg(Field_id, APC_VERTEX, el_nodes[ino + 1]) > 0 ) {
			/* real node */

			/*kbw
			if(El_id==2232){
				  printf("At real node %d, pdeg %d, constr_id %d\n",
					 el_nodes[ino+1],
					 apr_get_ent_pdeg(Field_id,APC_VERTEX,el_nodes[ino+1]),
					 num_shap_constr);
			 }
			/*kew*/

			List_dof_ent_type[num_shap_constr] = APC_VERTEX;
			List_dof_ent_id[num_shap_constr] = el_nodes[ino + 1];
			List_dof_ent_nrdofs[num_shap_constr] = nreq;
			num_shap_constr++;
			constr_el_loc[ino + 1][0] = 1.0;
			el_nodes_constr[ino + 1] = 1;
		}
		else if (apr_get_ent_pdeg(Field_id, APC_VERTEX, el_nodes[ino + 1]) == 0 ) {
			/* constrained node */
			is_constrained = 1;
			apr_get_constr_data(Field_id, el_nodes[ino + 1], APC_VERTEX, constr, NULL);
			el_nodes_constr[ino + 1] = constr[0];
			if (constr[0] == 2) constr_val = 0.5;
			if (constr[0] == 4) constr_val = 0.25;
			/*kbw
			if(El_id==2232){
				  printf("At constrained node %d, pdeg %d, nr_constr %d\n",
					 el_nodes[ino+1],
					 apr_get_ent_pdeg(Field_id,APC_VERTEX,el_nodes[ino+1]),
					 el_nodes_constr[ino+1]);
			 }
			/*kew*/
			for (j = 1; j <= constr[0]; j++) {
				/*kbw
				if(El_id==2232){
					  printf("At constraining node %d, pdeg %d, constr_id %d, constr_val %lf\n",
						 constr[j],
						 apr_get_ent_pdeg(Field_id,APC_VERTEX,constr[j]),
						 num_shap_constr, constr_val);
				 }
				/*kew*/
				List_dof_ent_type[num_shap_constr] = APC_VERTEX;
				List_dof_ent_id[num_shap_constr] = constr[j];
				List_dof_ent_nrdofs[num_shap_constr] = nreq;
				num_shap_constr++;
				constr_el_loc[ino + 1][j - 1] = constr_val;
			}
		}

	}

	num_dofs_constr = num_shap_constr * nreq;

	/*kbw
	if(El_id==2232){
	printf("In apr_get_stiff_mat_data: field_id %d, El_id %d, Comp_sm %d, Nr_dof_ent %d\n",
		   Field_id, El_id, Comp_sm, *Nr_dof_ent);
	 printf("For each node: \ttype, \tid, \tnrdofs\n");
	 for(i=0;i<num_shap_constr;i++){
	   printf("\t\t%d\t%d\t%d\n",
		  List_dof_ent_type[i],List_dof_ent_id[i],List_dof_ent_nrdofs[i]);
	 }
	 getchar();getchar();
	 }
	/*kew*/

	if (*Nr_dof_ent < num_shap_constr) {
		printf("Too small arrays List_dof_... passed to apr_get_stiff_mat_data\n");
		printf("%d < %d. Exiting!\n", *Nr_dof_ent, num_shap_constr);
		exit(-1);
	}

	*Nr_dof_ent = num_shap_constr;

	if (Comp_sm != APC_NO_COMP) {

		//double constr_value[4*MMC_MAXELVNO];

		if (*Nrdofs_loc < nreq * num_shap_constr) {
			printf(
					"Too small arrays Stiff_mat and Rhs_vect passed to apr_get_stiff_mat_data\n");
			printf("%d < %d. Exiting !!!", *Nrdofs_loc, num_shap_constr * nreq);
			exit(-1);
		}

		*Nrdofs_loc = num_dofs_constr;

		if (Comp_sm == APC_COMP_SM || Comp_sm == APC_COMP_RHS || Comp_sm == APC_COMP_BOTH ) {

			printf("May be some day... (APC_COMP_... not implemented yet\n");
			exit(-1);

		}

		if (is_constrained != 0) {

			int i_constr, j_constr, kk, kk_constr, idofs, jdofs;
			int nci, ncj, ieq, jeq;

			if (Comp_sm == APC_REWR_SM || Comp_sm == APC_REWR_BOTH) {

				/* rewrite Stiff_mat to local array */
				for (i = 0; i < num_dofs * num_dofs; i++) stiff_loc[i] = Stiff_mat[i];
				for (i = 0; i < num_dofs_constr * num_dofs_constr; i++) Stiff_mat[i] = 0;

				j_constr = 0;
				kk = 0;
				kk_constr = 0;
				for (jdofs = 0; jdofs < num_shap; jdofs++) {

					for (ncj = 0; ncj < el_nodes_constr[jdofs + 1]; ncj++) {

						i_constr = 0;
						for (idofs = 0; idofs < num_shap; idofs++) {

							for (nci = 0; nci < el_nodes_constr[idofs + 1]; nci++) {

								for (jeq = 0; jeq < nreq; jeq++) {

									for (ieq = 0; ieq < nreq; ieq++) {

										Stiff_mat[kk_constr + i_constr * nreq + jeq * num_dofs_constr + ieq] +=
												constr_el_loc[jdofs + 1][ncj] * constr_el_loc[idofs + 1][nci] *
												stiff_loc[kk + idofs * nreq + jeq * num_dofs + ieq];

										/*kbw
										printf("rewriting:\nold: \tjdofs %d, idofs %d, kk %d, index %d value %lf\n",
											   jdofs, idofs, kk, kk+idofs*nreq, stiff_loc[kk+idofs*nreq]);
										printf("new: \tncj %d, nci %d, kk_constr %d, i_constr %d, index %d, constr %lf, value %lf\n",
											   ncj, nci, kk_constr, i_constr,
											   kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq,
											   constr_el_loc[jdofs+1][ncj]*constr_el_loc[idofs+1][nci],
											   Stiff_mat[kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq]);
										/*kew*/

									}

								}

								i_constr++;

							}

						}/* idofs */

						j_constr++;
						kk_constr += nreq * num_dofs_constr;

					}

					kk += nreq * num_dofs;

				} /* jdofs */


			} // end if rewriting SM

			if (Comp_sm == APC_REWR_RHS || Comp_sm == APC_REWR_BOTH) {

				/* rewrite Rhs_vect to local array */
				for (i = 0; i < num_dofs; i++) rhs_loc[i] = Rhs_vect[i];
				for (i = 0; i < num_dofs_constr; i++) Rhs_vect[i] = 0;

				i_constr = 0;
				for (idofs = 0; idofs < num_shap; idofs++) {

					for (nci = 0; nci < el_nodes_constr[idofs + 1]; nci++) {

						for (ieq = 0; ieq < nreq; ieq++) {

							Rhs_vect[i_constr * nreq + ieq] +=
									constr_el_loc[idofs + 1][nci] * rhs_loc[idofs * nreq + ieq];

						}

						i_constr++;

					}

				}/* idofs */

			} // end if rewriting RHS

		} // end if constrained nodes exist

	} // end if considering SM and RHS

	return (1);
}


/*------------------------------------------------------------
  apr_proj_dof_ref - to rewrite dofs after modifying the mesh
			  (the procedure also recreates constraints data)
------------------------------------------------------------*/
int apr_proj_dof_ref(
		int Field_id,    /* in: approximation field ID  */
		int El,	   /* in:  >0 - rewrite after one [de]refinement of el */
		/*     <=0 - rewrite after massive [de]refinements */
		int Max_elem_id_before, /* in: maximal element (face, etc.) id */
		int Max_face_id_before, /*     before !!!!!!! refinements */
		int Max_edge_id_before,
		int Max_vert_id_before
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id; //,el_fath;
	int nvert, nreq, nr_dof_ents;
//  int fa_edges[5][5];
	int edge_nodes[3];
	int edge_nodes2[3];
//  double vect_dofs[3+1];
//  double vect_dofs2[3];

	int ned, edge_sons[3], mid_node, i, iaux;
	int num_nodes, fa_nodes[5], fa_sons[5], ino, nfa;
	int * edge_elems = (int *) malloc((MMC_MAX_EDGE_ELEMS + 1) * sizeof(int));

	int node_shift;
	int face_neig[2];	/* list of neighbors */
	int neig_sides[2] = {0, 0}; /* sides of face wrt neigs */

	int max_elem_id, max_face_id, max_edge_id, max_vert_id;

	int repeat;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	/* select the pointer to the approximation field */
	field_p = apr_select_field(Field_id);

	// check if mesh module is ready
	if (mmr_is_ready_for_proj_dof_ref(mesh_id) < 0) {
		printf("Mesh module not ready for dofs projection!");
		exit(-1);
	}

	max_elem_id = mmr_get_max_elem_id(mesh_id);
	if (Max_elem_id_before > max_elem_id) max_elem_id = Max_elem_id_before;
	max_face_id = mmr_get_max_face_id(mesh_id);
	if (Max_face_id_before > max_face_id) max_face_id = Max_face_id_before;
	max_edge_id = mmr_get_max_edge_id(mesh_id);
	if (Max_edge_id_before > max_edge_id) max_edge_id = Max_edge_id_before;
	max_vert_id = mmr_get_max_node_id(mesh_id);
	if (Max_vert_id_before > max_vert_id) max_vert_id = Max_vert_id_before;

	/* may be one pass suffices? */
	repeat = 0;

	// check if there is a need of resizing dofs structure
	if (max_vert_id > field_p->capacity_dof_ents) {

		printf("Sizing dof_ent to %d \n", max_vert_id);

		apt_dof_ent * new_dof_ents
			= (apt_dof_ent *) malloc((2 * max_vert_id + 1) * sizeof(apt_dof_ent));
		if (new_dof_ents == NULL) {
			printf("Dofs structures not reallocated (new size:%d)!", 2 * max_vert_id + 1);
			exit(-1);
		}

		// fast rewriting existing solution
		memcpy(new_dof_ents, field_p->dof_ents,
				(Max_vert_id_before + 1)*sizeof(apt_dof_ent));

		// fast NULL-ing new space for dof_ents
		memset(new_dof_ents + (Max_vert_id_before + 1)*sizeof(apt_dof_ent),
				0, 2 * max_vert_id - Max_vert_id_before);
		// KB: shouldn't it be multiplied by sizeof(apt_dof_ent)
		printf("usage of memset in apr_proj_dof_ref?????\n");
		getchar();
		getchar();
		getchar();

		free(field_p->dof_ents);
		field_p->dof_ents = new_dof_ents;
		field_p->capacity_dof_ents = 2 * max_vert_id;
	}

	nreq = apr_get_nreq(Field_id);

	/* to begin indicate that there are NO constrained nodes in the mesh */
	field_p->constr = APC_FALSE;

	/* allocate the space for elems structure for active and inactive edges */
	mmr_create_edge_elems(mesh_id, max_edge_id);

	/*kbw
	  printf("in  apr_proj_dof_ref: Field_id %d, mesh_id %d\n", Field_id, mesh_id);
	  printf("max_elem_id %d, max_face_id %d, max_edge_id %d, max_vert_id %d\n",
		 max_elem_id, max_face_id, max_edge_id, max_vert_id);
	/*kew*/

	//ned = 0;
	//while((ned = mmr_get_next_edge_all(mesh_id,ned))!=0){
	/* for all edges in the current and the previous mesh */
	for (ned = 1; ned <= max_edge_id; ned++) {

		/*kbw
			if(mmr_edge_status(mesh_id,ned)!=MMC_FREE){
			  int eno[2];
			  mmr_edge_nodes(mesh_id, ned, eno);
			  if(eno[0]==950 || eno[1]==950){
			printf(
				   "In proj_dof_ref node 950 appear: edge %d, nodes %d, %d\n",
				   ned, eno[0], eno[1]);
			if(mmr_edge_status(mesh_id,ned)==MMC_INACTIVE){
			  printf("edge %d is inactive\n",ned);
			}
			if(apr_get_ent_pdeg(mesh_id,APC_VERTEX,eno[0])==0){
			  printf("node %d is constrained\n",eno[0]);
			}
			if(apr_get_ent_pdeg(mesh_id,APC_VERTEX,eno[1])==0){
			  printf("node %d is constrained\n",eno[1]);
			}
			//getchar();getchar();
			  }
			}
		/*kew*/

		if (mmr_edge_status(mesh_id, ned) == MMC_INACTIVE) {

			mmr_edge_sons(mesh_id, ned, edge_sons, NULL);

			/*kbw
			printf("in parent edge %d, sons %d %d\n",
			   ned, edge_sons[0], edge_sons[1]);
			/*kew*/

			/* the node in the middle of a divided edge */
			mmr_edge_nodes(mesh_id, edge_sons[0], edge_nodes2);
			nvert = edge_nodes2[1];
			assert(nvert <= mmr_get_max_node_id(mesh_id));
			assert(nvert <= mmr_get_max_node_max(mesh_id));
			assert(nvert <= field_p->capacity_dof_ents);

			/*kbw
			printf("first son nodes %d %d\n", edge_nodes2[0], edge_nodes2[1]);
			/*kew*/

			/* for inactive edges check whether all surrounding elements are active */
			mmr_edge_elems(mesh_id, ned, edge_elems);

			iaux = 0;
			for (i = 1; i <= edge_elems[0]; i++) {
				if (mmr_el_status(mesh_id, edge_elems[i]) == MMC_INACTIVE) iaux++;
			}

			if (iaux == edge_elems[0]) {

				/* all edge elements are divided - a real node */

				/* check the status in the previous mesh */
				if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {

					/* currently REAL node, previously REAL node  - do nothing */
					/*kbw
						  printf("currently REAL node, previously REAL node %d\n", nvert);
					/*kew*/

				}
				else {

					if (field_p->dof_ents[nvert].constr != NULL) {

						/* currently REAL node, previously CONSTR node */
						/*kbw
							  printf("currently REAL node, previously CONSTR node %d\n", nvert);
						/*kew*/

						/* free constr table */
						free(field_p->dof_ents[nvert].constr);
						field_p->dof_ents[nvert].constr = NULL;

					}
					else {

						/* currently REAL node, previously FREE space */
						/*kbw
							  printf("currently REAL node, previously FREE space %d\n", nvert);
						/*kew*/

					}

					mmr_edge_nodes(mesh_id, ned, edge_nodes);

					/*kbw
					  printf("parent nodes %d %d\n", edge_nodes[0], edge_nodes[1]);
					/*kew*/

					/*kbw
					  if(nvert==2060){
						printf("parent nodes %d %d\n", edge_nodes[0], edge_nodes[1]);
					  }
					/*kew*/


#ifdef DEBUG_APM
					if (field_p->dof_ents[edge_nodes[0]].vec_dof_1 == NULL ||
							field_p->dof_ents[edge_nodes[1]].vec_dof_1 == NULL) {
						printf("False double constrained node %d (parent edge %d, parent nodes %d %d\n",
								nvert, ned, edge_nodes[0], edge_nodes[1]);
						printf("Repeating proj_dof_ref !\n");
						//exit(-1);
					}
#endif
					if (field_p->dof_ents[edge_nodes[0]].vec_dof_1 == NULL ||
							field_p->dof_ents[edge_nodes[1]].vec_dof_1 == NULL) {
						repeat = 1;
						continue;
					}

					/*kbw
						  for(i=0;i<nreq;i++){
							if(fabs(field_p->dof_ents[edge_nodes[0]].vec_dof_1[i]+
							   field_p->dof_ents[edge_nodes[1]].vec_dof_1[i])<1e-9){
					printf("Zero DOFs for node %d (parent edge %d) in parent nodes %d %d\n",
								 nvert, ned, edge_nodes[0], edge_nodes[1]);
							}
						  }
					/*kew*/

					/* allocate dofs vectors */
					field_p->dof_ents[nvert].vec_dof_1 =
							(double *) malloc(nreq * sizeof(double));
					if (field_p->nr_sol > 1) {
						field_p->dof_ents[nvert].vec_dof_2 =
								(double *) malloc(nreq * sizeof(double));
					}
					if (field_p->nr_sol > 2) {
						field_p->dof_ents[nvert].vec_dof_3 =
								(double *) malloc(nreq * sizeof(double));
					}

					for (i = 0; i < nreq; i++) {
						field_p->dof_ents[nvert].vec_dof_1[i] =
								0.5 * field_p->dof_ents[edge_nodes[0]].vec_dof_1[i] +
								0.5 * field_p->dof_ents[edge_nodes[1]].vec_dof_1[i];
					}

					if (field_p->nr_sol > 1) {
						for (i = 0; i < nreq; i++) {
							field_p->dof_ents[nvert].vec_dof_2[i] =
									0.5 * field_p->dof_ents[edge_nodes[0]].vec_dof_2[i] +
									0.5 * field_p->dof_ents[edge_nodes[1]].vec_dof_2[i];
						}
					}
					if (field_p->nr_sol > 2) {
						for (i = 0; i < nreq; i++) {
							field_p->dof_ents[nvert].vec_dof_3[i] =
									0.5 * field_p->dof_ents[edge_nodes[0]].vec_dof_3[i] +
									0.5 * field_p->dof_ents[edge_nodes[1]].vec_dof_3[i];
						}
					}

					/*kbw
					for(i=0;i<nreq;i++){
					  printf("edge nodes' dofs %lf %lf\n",
						 field_p->dof_ents[edge_nodes[0]].vec_dof_1[i],
						 field_p->dof_ents[edge_nodes[1]].vec_dof_1[i]);
					  printf("new node's dofs %lf\n",
						 field_p->dof_ents[nvert].vec_dof_1[i]);
					}
					/*kew*/

				} /* end if REAL node, previously not REAL node */

			} /* end if REAL node */
			else {

				/* at least one edge element not divided - a constrained node */

				/* indicate that there are constrained nodes in the mesh */
				field_p->constr = APC_TRUE;

				/* check the status in the previous mesh */
				if (field_p->dof_ents[nvert].constr != NULL) {

					/* currently CONSTR node, previously CONSTR node - do nothing */
					/*kbw
						  printf("currently CONSTR node, previously CONSTR node %d\n", nvert);
					/*kew*/

				}
				else {

					if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {

						/*kbw
							  printf("currently CONSTR node, previously REAL node %d\n", nvert);
						/*kew*/
						/* currently CONSTR node, previously REAL node */
						free(field_p->dof_ents[nvert].vec_dof_1);
						field_p->dof_ents[nvert].vec_dof_1 = NULL;
						if (field_p->nr_sol > 1) {
							free(field_p->dof_ents[nvert].vec_dof_2);
							field_p->dof_ents[nvert].vec_dof_2 = NULL;
						}
						if (field_p->nr_sol > 2) {
							free(field_p->dof_ents[nvert].vec_dof_3);
							field_p->dof_ents[nvert].vec_dof_3 = NULL;
						}

					}
					else {

						/* currently CONSTR node, previously FREE space */
						/*kbw
							  printf("currently CONSTR node, previously FREE space %d\n", nvert);
						/*kew*/

					}

					field_p->dof_ents[nvert].constr =
							(int *) malloc(3 * sizeof(int));
					field_p->dof_ents[nvert].constr[0] = 2;


					mmr_edge_nodes(mesh_id, ned, edge_nodes);

					/*kbw
						printf("parent nodes %d %d\n", edge_nodes[0], edge_nodes[1]);
					/*kew*/

					field_p->dof_ents[nvert].constr[1] = edge_nodes[0];
					field_p->dof_ents[nvert].constr[2] = edge_nodes[1];


				} /* end if CONSTR node previously not CONSTR node */

			} /* end if CONSTR node, after checking the current status */
		} /* end if edge inactive */
	} /* end loop over edges */


	/* for new mid-face nodes */
	//nfa=0;
	//while((nfa=mmr_get_next_face_all(mesh_id,nfa))!=0){
	for (nfa = 1; nfa <= max_face_id; nfa++) {

		if (mmr_fa_type(mesh_id, nfa) == MMC_QUAD &&
				mmr_fa_status(mesh_id, nfa) == MMC_INACTIVE) {

			//printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

			mmr_fa_fam(mesh_id, nfa, fa_sons, &mid_node);

			/*kbw
			  printf("in big inactive face %d, with mid_node %d without dofs\n",
			  nfa, mid_node);
			/*kew*/

			/* for divided quadrilateral faces */
			if (mid_node > 0 && mmr_node_status(mesh_id, mid_node) == MMC_ACTIVE) {

				/* active node in the middle */
				nvert = mid_node;

				/* check the status of the node */
				mmr_fa_neig(mesh_id, nfa, face_neig, neig_sides, &node_shift, NULL, NULL, NULL);

				/* if face between two elements and only one divided (inactive) */
				if (face_neig[0] != 0 && face_neig[1] != 0 &&
						( ( mmr_el_status(mesh_id, face_neig[0]) == MMC_ACTIVE &&
										mmr_el_status(mesh_id, face_neig[1]) == MMC_INACTIVE ) ||
								( mmr_el_status(mesh_id, face_neig[0]) == MMC_INACTIVE &&
										mmr_el_status(mesh_id, face_neig[1]) == MMC_ACTIVE ) ) ) {

					/* indicate that there are constrained nodes in the mesh */
					field_p->constr = APC_TRUE;

					/* check the status in the previous mesh */
					if (field_p->dof_ents[nvert].constr != NULL) {

						/* currently CONSTR node, previously CONSTR node - do nothing */
						/*kbw
							  printf("currently CONSTR node, previously CONSTR node %d\n", nvert);
						/*kew*/

					}
					else {

						/* currently CONSTR node, previously REAL node  */
						if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {

							/*kbw
								  printf("currently CONSTR node, previously REAL node %d\n", nvert);
							/*kew*/

							free(field_p->dof_ents[nvert].vec_dof_1);
							field_p->dof_ents[nvert].vec_dof_1 = NULL;
							if (field_p->nr_sol > 1) {
								free(field_p->dof_ents[nvert].vec_dof_2);
								field_p->dof_ents[nvert].vec_dof_2 = NULL;
							}
							if (field_p->nr_sol > 2) {
								free(field_p->dof_ents[nvert].vec_dof_3);
								field_p->dof_ents[nvert].vec_dof_3 = NULL;
							}

						}
						else {

							/* currently CONSTR node, previously FREE space */
							/*kbw
								  printf("currently CONSTR node, previously FREE space %d\n", nvert);
							/*kew*/

						}

						field_p->dof_ents[nvert].constr =
								(int *) malloc(5 * sizeof(int));
						field_p->dof_ents[nvert].constr[0] = 4;


						num_nodes = mmr_fa_node_coor(mesh_id, nfa, fa_nodes, NULL);

						for (ino = 1; ino <= num_nodes; ino++) {

							/*kbw
							  printf("%d face node %d, \n", ino,fa_nodes[ino]);
							  /*kew*/

							field_p->dof_ents[nvert].constr[ino] = fa_nodes[ino];

						}

					} /* end if CONSTR node, previously not CONSTR node */

				} /* end if currently constrained active node */
				else {

					/* check the status in the previous mesh */
					if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {

						/* currently REAL node, previously REAL node  - do nothing */
						/*kbw
							  printf("currently REAL node, previously REAL node %d\n", nvert);
						/*kew*/

					}
					else {

						if (field_p->dof_ents[nvert].constr != NULL) {

							/* currently REAL node, previously CONSTR node */
							/*kbw
								  printf("currently REAL node, previously CONSTR node %d\n", nvert);
							/*kew*/

							/* free constr table */
							free(field_p->dof_ents[nvert].constr);
							field_p->dof_ents[nvert].constr = NULL;

						}
						else {

							/* currently REAL node, previously FREE space */
							/*kbw
								  printf("currently REAL node, previously FREE space %d\n", nvert);
							/*kew*/

						}



						field_p->dof_ents[nvert].vec_dof_1 =
								(double *) malloc(nreq * sizeof(double));
						if (field_p->nr_sol > 1) {
							field_p->dof_ents[nvert].vec_dof_2 =
									(double *) malloc(nreq * sizeof(double));
						}
						if (field_p->nr_sol > 2) {
							field_p->dof_ents[nvert].vec_dof_3 =
									(double *) malloc(nreq * sizeof(double));
						}
						for (i = 0; i < nreq; i++) {
							field_p->dof_ents[nvert].vec_dof_1[i] = 0.0;
						}
						if (field_p->nr_sol > 1) {
							for (i = 0; i < nreq; i++) {
								field_p->dof_ents[nvert].vec_dof_2[i] = 0.0;
							}
						}
						if (field_p->nr_sol > 2) {
							for (i = 0; i < nreq; i++) {
								field_p->dof_ents[nvert].vec_dof_3[i] = 0.0;
							}
						}
						num_nodes = mmr_fa_node_coor(mesh_id, nfa, fa_nodes, NULL);

#ifdef DEBUG_APM
						if (field_p->dof_ents[fa_nodes[1]].vec_dof_1 == NULL ||
								field_p->dof_ents[fa_nodes[2]].vec_dof_1 == NULL ||
								field_p->dof_ents[fa_nodes[3]].vec_dof_1 == NULL ||
								field_p->dof_ents[fa_nodes[4]].vec_dof_1 == NULL) {
							printf("False double constrained node %d (parent face %d).\n",
									nvert, nfa);
							printf("Repeating proj_dof_ref !\n");
							//exit(-1);
						}
#endif

						if (field_p->dof_ents[fa_nodes[1]].vec_dof_1 == NULL ||
								field_p->dof_ents[fa_nodes[2]].vec_dof_1 == NULL ||
								field_p->dof_ents[fa_nodes[3]].vec_dof_1 == NULL ||
								field_p->dof_ents[fa_nodes[4]].vec_dof_1 == NULL) {
							repeat = 1;
							continue;
						}

						for (i = 0; i < nreq; i++) {

							for (ino = 1; ino <= num_nodes; ino++) {

								/*kbw
										  printf("face node %d, ieq %d, dofs %lf\n", ino, i,
										  field_p->dof_ents[fa_nodes[ino]].vec_dof_1[i]);
										  printf("new node's dofs before update %lf\n",
										  field_p->dof_ents[nvert].vec_dof_1[i]);
								/*kew*/
								field_p->dof_ents[nvert].vec_dof_1[i] +=
										0.25 * field_p->dof_ents[fa_nodes[ino]].vec_dof_1[i];

								/*kbw
										  printf("new node's dofs after update %lf\n",
										  field_p->dof_ents[nvert].vec_dof_1[i]);
								/*kew*/

								if (field_p->nr_sol > 1) {
									field_p->dof_ents[nvert].vec_dof_2[i] +=
											0.25 * field_p->dof_ents[fa_nodes[ino]].vec_dof_2[i];
								}
								if (field_p->nr_sol > 2) {
									field_p->dof_ents[nvert].vec_dof_3[i] +=
											0.25 * field_p->dof_ents[fa_nodes[ino]].vec_dof_3[i];
								}

							} /* end loop over parent nodes */

						} /* end loop over equations (solution components) */

					} /* and new REAL node previously not REAL node */

				} /* and REAL node */

			} /* end active mid-edge node */

		} /* end inactive quadrilateral face */

	} /* end loop over faces */


	/* deleting dofs of derefined nodes */
	//nvert=0;
	//while((nvert=mmr_get_next_node_all(mesh_id,nvert))!=0)
	for (nvert = 1; nvert <= max_vert_id; nvert++) {

		if (mmr_node_status(mesh_id, nvert) == MMC_FREE) {

			if ( field_p->dof_ents[nvert].vec_dof_1 != NULL) {

				/* currently FREE space, previously REAL node */
				free(field_p->dof_ents[nvert].vec_dof_1);
				field_p->dof_ents[nvert].vec_dof_1 = NULL;
				if (field_p->nr_sol > 1) {
					free(field_p->dof_ents[nvert].vec_dof_2);
					field_p->dof_ents[nvert].vec_dof_2 = NULL;
				}
				if (field_p->nr_sol > 2) {
					free(field_p->dof_ents[nvert].vec_dof_3);
					field_p->dof_ents[nvert].vec_dof_3 = NULL;
				}
			}
			else {

				if ( field_p->dof_ents[nvert].constr != NULL) {
					/* currently FREE space, previously CONSTR node */
					free(field_p->dof_ents[nvert].constr);
					field_p->dof_ents[nvert].constr = NULL;
				}
				else {
					/* currently FREE space, previously FREE space - do nothing */
				}


			}
		}
	}

	nr_dof_ents = 0;
	for (nvert = 1; nvert <= max_vert_id; nvert++) {
		/* for active unconstrained vertices */
		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE &&
				field_p->dof_ents[nvert].vec_dof_1 != NULL) nr_dof_ents++;
	}
	field_p->nr_dof_ents = nr_dof_ents;

	//printf("repeat= %d\n",repeat);
	if (repeat == 1) {
#ifdef DEBUG_APM
		printf("Repeating proj_dof_ref! nr_dof_ents = %d\n", nr_dof_ents);
#endif
		apr_proj_dof_ref(Field_id, El, Max_elem_id_before, Max_face_id_before,
				Max_edge_id_before, Max_vert_id_before);
	}

#ifdef DEBUG_APM
	/* checking */
	//nvert=0;
	//while((nvert=mmr_get_next_node_all(mesh_id,nvert))!=0)
	for (nvert = 1; nvert <= max_vert_id; nvert++) {

		/*kbw
			if(nvert==362){
			  printf("node %d: status %d, dofp %lu, constrp %lu\n",
				 nvert, mmr_node_status(mesh_id,nvert),
				 field_p->dof_ents[nvert].vec_dof_1,
				 field_p->dof_ents[nvert].constr);
			}
		/*kew*/

		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE)
		{
			if ( field_p->dof_ents[nvert].vec_dof_1 == NULL &&
					field_p->dof_ents[nvert].constr == NULL ) {

				printf("Node %d without dofs after proj_dof_ref !\n", nvert);
				exit(-1);

			}
		}
	}

	apr_check_field(Field_id);

#endif

	free(edge_elems);

	return (1);
}


/*------------------------------------------------------------
  apr_rewr_sol - to rewrite solution from one vector to another
------------------------------------------------------------*/
int apr_rewr_sol(
		int Field_id,      /* in: data structure to be used  */
		int Sol_from,      /* in: ID of vector to read solution from */
		int Sol_to         /* in: ID of vector to write solution to */
)
{

	int nno, mesh_id, num_shap, nreq; // Changed on 12.2010 - nreq added
	double dofs_loc[APC_MAXELSD];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	mesh_id = apr_get_mesh_id(Field_id);
	nreq = apr_get_nreq(Field_id); // Changed on 12.2010 - line added

	nno = 0;
	while ((nno = mmr_get_next_node_all(mesh_id, nno)) != 0) {

		if (apr_get_ent_pdeg(Field_id, APC_VERTEX, nno) > 0) {
			// faster: select field_p and next
			// if(field_p->dof_ents[Ent_id].vec_dof_1==NULL){

			num_shap = apr_get_ent_numshap(Field_id, APC_VERTEX, nno);
			apr_read_ent_dofs(Field_id, APC_VERTEX, nno, // Changed on 12.2010
					num_shap * nreq, Sol_from, dofs_loc); // - num_shap*nreq
			apr_write_ent_dofs(Field_id, APC_VERTEX, nno, // Changed on 12.2010
					num_shap * nreq, Sol_to, dofs_loc); // - num_shap*nreq

		}
	}

	return 0;
}


/*------------------------------------------------------------
  apr_free_field - to free approximation field data structure
------------------------------------------------------------*/
int apr_free_field(
		int Field_id    /* in: approximation field ID  */
)
{
	/* pointer to field structure */
	apt_field * field_p;
	int i, nmvert, mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	/* loop over dofs data structure including free spaces */
	nmvert = mmr_get_max_node_id(mesh_id);
	for (i = 1; i <= nmvert; i++) {
		/* clean the space if it has dofs */
		if (field_p->dof_ents[i].vec_dof_1 != NULL) {
			free(field_p->dof_ents[i].vec_dof_1);
			if (field_p->nr_sol > 1) free(field_p->dof_ents[i].vec_dof_2);
			if (field_p->nr_sol > 2) free(field_p->dof_ents[i].vec_dof_3);
		}
		if (field_p->dof_ents[i].constr != NULL) {
			free(field_p->dof_ents[i].constr);
		}
	}
	free(field_p->dof_ents);

	return (1);
}


/*------------------------------------------------------------
  apr_limit_deref - to return whether derefinement is allowed or not
	the routine does not use approximation data structures but the result
	depends on the approximation method
------------------------------------------------------------*/
int apr_limit_deref(
		int Field_id,    /* in: approximation field ID  */
		int El_id      /* in: unique identifier of the element */
)
{
	int ied, ned, el_edges[13], el_edges_2[13];
	int i, j;
	int edge_elems[20 + 1];
	int elsons[MMC_MAXELSONS + 1]; /* family information */
	int mesh_id, father, son, ison, max_gen_diff, ignore_edges, ref_type;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	mesh_id = apr_get_mesh_id(Field_id);
	max_gen_diff = mmr_get_max_gen_diff(mesh_id);

	if (max_gen_diff == 1) {

		father = mmr_el_fam(mesh_id, El_id, NULL, NULL);
		mmr_el_fam(mesh_id, father, elsons, NULL);

		/* we can cluster elements only if all their edges do not contain */
		/* constrained nodes, i.e. are active */
		/* and also edges of all their neighbours across edges... */
		ref_type = mmr_el_type_ref(mesh_id, father);

		if (ref_type == MMC_REF_ANI) ignore_edges = 3;
		else ignore_edges = 0;

		for (ison = 1; ison <= elsons[0]; ison++) {

			son = elsons[ison];

			mmr_el_edges(mesh_id, son, el_edges);

			/* for all edges */
			for (ied = 1; ied <= el_edges[0] - ignore_edges; ied++) {

				ned = el_edges[ied];

				/* find elements containing the edge */
				mmr_edge_elems(mesh_id, ned, edge_elems);

				for (i = 1; i <= edge_elems[0]; i++) {

					mmr_el_edges(mesh_id, edge_elems[i], el_edges_2);

					for (j = 1; j < el_edges_2[0] - ignore_edges; j++) {

						if (mmr_edge_status(mesh_id, el_edges_2[j])
								== MMC_INACTIVE) {

							return (APC_DEREF_DENIED);

						}
					} /* end loop over edges */
				} /* end loop over elements sharing edges with sons */
			} /* end loop over sons' edgse */
		} /* end loop over sons */
	} /* end if(max_gen_diff==1) */
	return (APC_DEREF_ALLOWED);
}

/*------------------------------------------------------------
  apr_limit_ref - to return whether refinement is allowed or not
	the routine does not use approximation data structures but the result
	depends on the approximation method
------------------------------------------------------------*/
int apr_limit_ref(
		int Field_id,    /* in: approximation field ID  */
		int El_id      /* in: unique identifier of the element */
)
{

	int mesh_id, gen_el, max_gen;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	mesh_id = apr_get_mesh_id(Field_id);
	max_gen = mmr_get_max_gen(mesh_id);

	if (max_gen == 0) return (APC_REF_DENIED);
	else if (max_gen > 0) {

		/* check whether max generation is not exceeded */
		gen_el = mmr_el_gen(mesh_id, El_id);

		if (gen_el < max_gen) return (APC_REF_ALLOWED);
		else return (APC_REF_DENIED);

	}

	return (APC_REF_ALLOWED);
}

/*------------------------------------------------------------
  apr_refine - to refine an element or the whole mesh checking mesh irregularity
------------------------------------------------------------*/
int apr_refine( /* returns: >=0 - success code, <0 - error code */
		int Field_id,  /* in: field ID */
		int El   /* in: element ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
)
{

	/* local variables */
	int i, iaux, iel, nmel_old, nrel_old, nrel_div, nelref, nrwait, listwait[100];
	int el_type, gen_el, max_gen, max_gen_diff;
	int ifa, num_face, face, face_sons[5], ifaneig;
	int iprint = MMC_PRINT_INFO + 1;
	int mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh_id = apr_get_mesh_id(Field_id);

	max_gen = mmr_get_max_gen(mesh_id);
	max_gen_diff = mmr_get_max_gen_diff(mesh_id);

	if (El != MMC_DO_UNI_REF && El > 0) {

		/*set number of elements waiting for refinement (they will be stored in listwait)*/
		nrwait = 0;

		/*element to refine*/
		nelref = El;

beginning:
		{}

		/* check whether element is active */
		if (mmr_el_status(mesh_id, nelref) != MMC_ACTIVE) {
#ifdef DEBUG_APM
			if (iprint > MMC_PRINT_INFO) {
				printf("not active element or free space %d in refine!\n", nelref);
			}
#endif
			return (-2);

			/*kb!!!
			if(nrwait==0) return(-2);
			nrwait--;
			nelref=listwait[nrwait];
			goto beginning;
			!!!kb*/
		}

		/* barrier */
		gen_el = mmr_el_gen(mesh_id, nelref);

		if (gen_el >= max_gen) {

#ifdef DEBUG_APM
			if (iprint > MMC_PRINT_INFO) {
				printf("element %d generation (%d) > max_gen (%d) in divide8_p!\n",
						nelref, gen_el, max_gen);
			}
#endif
			return (-3);
			/*kb!!!
			if(nrwait==0) return(-3);
			nrwait--;
			nelref=listwait[nrwait];
			goto beginning;
			!!!kb*/
		}


		/* barrier */
		if (max_gen_diff == 0) {
			El = MMC_DO_UNI_REF;
		}
		else if (max_gen_diff == 1) {

			/* 1-irregular meshes - no double constrained nodes */

			/* proper version: */
			/* 1. find the parent */
			/* 2. check the position within parent */
			/* 3. identify parent's edges concerned  */
			/* 4. check the respective parent's edges do not contain constrained nodes */

			/* simplified version with over-refinements: */
			/* all parent's edges must not contain constrained nodes */

			int parent;
			int ned, ied, nr_edges, el_edges[13];
			int edge_elems[20 + 1];

			parent = mmr_el_fam(mesh_id, nelref, NULL, NULL);

			if (parent != MMC_NO_FATH) {

				mmr_el_edges(mesh_id, parent, el_edges);

				/* for all parent's edges */
				nr_edges = el_edges[0];
				/* for 2D refinements we check only horizontal edges */
				if (mmr_el_type_ref(mesh_id, parent) == MMC_REF_ANI) {
					nr_edges = 6;
				}

				for (ied = 1; ied <= nr_edges; ied++) {

					ned = el_edges[ied];

					mmr_edge_elems(mesh_id, ned, edge_elems);

					for (i = 1; i <= edge_elems[0]; i++) {

						/*kbw
						if(nelref==724||parent==724){
							printf("Edge %d edge_elems: ", ned);
							for(i=1; i<=edge_elems[0]; i++){
							  printf("%d  ", edge_elems[i]);
							}
							printf("\n");
						  }
						/*kew*/
						if ( edge_elems[i] != nelref &&
								mmr_el_status(mesh_id, edge_elems[i]) == MMC_ACTIVE) {

							listwait[nrwait] = nelref;
							nrwait++;
							nelref = edge_elems[i];


							/*kbw
							#ifdef DEBUG_APM
							  if(iprint>MMC_PRINT_INFO){
								printf("adding element %d to refine because of DIFF_GEN\n",nelref);
								printf("nelref %d, parent %d, edge %d, elem_new %d\n",
								   listwait[nrwait-1], parent, ned, edge_elems[i]);

							  }
							#endif
							/*kew*/

							goto beginning;

						}
					}
				}
			} /* end if not initial mesh element */

		} /* end if max_gen_diff == -1 i.e. 1-irregular meshes */
		else {
			printf("Wrong parameter max_gen_diff %d in apr_refine! Exiting.\n"
					, max_gen_diff);
			exit(-1);
		}


		/* divide an element */
		iaux = mmr_refine_el(mesh_id, nelref);

		if (iaux < 0) return (iaux - 10);

		/*get the next element from the list (if any)*/
		if (nrwait == 0) return (iaux);
		nrwait--;
		nelref = listwait[nrwait];
		goto beginning;

	}
	else if (El == MMC_DO_UNI_REF) {
		return (mmr_refine_mesh(mesh_id));
	}
	else {
		printf("Wrong parameter %d in apr_refine! Exiting.\n"
				, El);
		exit(-1);
	}


#ifdef DEBUG_APM
	printf("Error 230956 in apr_refine!\n");
#endif

	/* error condition - that point should not be reached */
	return (-1);
}

/*------------------------------------------------------------
  apr_derefine - to derefine an element or the whole mesh with irregularity check
------------------------------------------------------------*/
int apr_derefine( /* returns: >=0 - success code, <0 - error code */
		int Field_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  /* in: element ID or -2 (MMC_DO_UNI_DEREF) for uniform derefinement */
)
{

	/* local variables */
	int i, iaux, iel, nmel_old, nrel_old, nrel_div;
	int iprint = 5;
	int mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh_id = apr_get_mesh_id(Field_id);

	if (El != MMC_DO_UNI_DEREF && El > 0) {

		if (MMC_ACTIVE == mmr_el_status(mesh_id, El)
				&& APC_DEREF_ALLOWED == apr_limit_deref(Field_id, El) )
		{
			iaux = mmr_derefine_el(mesh_id, El);
		}
	}
	else if (El == MMC_DO_UNI_DEREF) {
		iaux = mmr_derefine_mesh(mesh_id);
	}
	else {
		printf("Wrong parameter %d in apr_derefine! Exiting.\n"
				, El);
		exit(-1);
	}

	return (iaux);
}

/** @}*/

#ifdef __cplusplus
}
#endif

