/************************************************************************
File aps_std_util.c - utilities functions necessary for standard continuous
						 approximation on prismatic and tetrahedral elements

Contains implementation of routines:
  apr_select_field - to select the proper field

  apr_create_constr_data - to create and fill constraint arrays
				   (assuming approximation and mesh data agree)
  apr_get_constr - to return constraints data for a node (dof entity)
  apr_get_el_constr_data - to return the number and the list of element's real
						   nodes, with the corresponding constraint coefficients

  apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
  apr_set_quadr_3D - to prepare quadrature data for a given element
  apr_set_quadr_2D - to prepare quadrature data for a given face
  apr_L2_proj - to L2 project a function onto an element
	given local coordinates within an element of the same family
  apr_sol_xglob - to return the solution at a point with global
	coordinates specified

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
	05.2011 - Kazimierz Michalik, hybrid version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

/* interface for all approximation modules */
// moved to aph_std_lin.h, where we redefine some defines - KM 06.2011
//#include "aph_intf.h"

/* internal header file for the std lin approximation module */
#include <modfem/ap_std_lin/aph_std_lin.h>

/* interface of the mesh manipulation module */
#include <modfem/mmh_intf.h>

/* LAPACK procedures */
#include <modfem/lin_alg_intf.h>

// interface of the vec3 functions
#include <modfem/uth_intf.h>

#include <modfem/uth_log.h>

/** \addtogroup APM_STD_LIN Standard Linear Approximation
 *  \ingroup APR
 *  @{
 */


/* internal utility procedures */
double ut_mat3_inv(
		/* returns: determinant of matrix to invert */
		double * mat,	/* matrix to invert */
		double * mat_inv	/* inverted matrix */
);
void ut_vec3_prod(
		double * vec_a, 	/* in: vector a */
		double * vec_b, 	/* in: vector b */
		double * vec_c	/* out: vector product axb */
);
double ut_vec3_mxpr( /* returns: mixed product [a,b,c] */
		double * vec_a, 	/* in: vector a */
		double * vec_b, 	/* in: vector b */
		double * vec_c	/* in: vector c */
);
double ut_vec3_length(	/* returns: vector length */
		double * vec	/* in: vector */
);
void ut_mat3vec(
		double * m1, 	/* in: matrix (stored by rows as a vector!) */
		double * v1, 	/* in: vector */
		double * v2	/* out: resulting vector */
);
void ut_mat3mat(
		double * m1,	/* in: matrix */
		double * m2,	/* in: matrix */
		double * m3	/* out: matrix m1*m2 */
);
int ut_chk_list(	/* returns: */
		/* >0 - position on the list */
		/* 0 - not found on the list */
		int Num, 	/* number to be checked */
		int * List, 	/* list of numbers */
		int Ll		/* length of the list */
);

/* standard macro for max and min and abs */
#define ut_max(x,y) ((x)>(y)?(x):(y))
#define ut_min(x,y) ((x)<(y)?(x):(y))
#define ut_abs(x)   ((x)<0?-(x):(x))

/*---------------------------------------------------------
  apr_select_field - to select the proper field
---------------------------------------------------------*/
apt_field * apr_select_field( /* returns: pointer to the selected field */
		int Field_id  /* in: field ID */
)
{

	if ( Field_id == APC_CUR_FIELD_ID ) {
		return (&apv_fields[apv_cur_field_id - 1]);
	}
	else if ( Field_id > 0 && Field_id <= apv_nr_fields ) {
		return (&apv_fields[Field_id - 1]);
	}
	else {
		return (&apv_fields[apv_cur_field_id - 1]);
		/* alternative:  return(NULL);   */
	}
}

/*------------------------------------------------------------
  apr_create_constr_data - to create and fill constraint arrays
				   (assuming approximation and mesh data agree)
------------------------------------------------------------*/
int apr_create_constr_data(
		int Field_id    /* in: approximation field ID  */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id; //,el_fath;
	int nvert;
//  int el_faces[5];
//  int fa_edges[5][5];
	int edge_nodes[3];
	int edge_nodes2[3];
	int ned, edge_sons[3], mid_node;
	int num_nodes, fa_nodes[5], fa_sons[5], ino, nfa;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the pointer to the approximation field */
	field_p = apr_select_field(Field_id);

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	// edge_elems data structure is assumed to be created in init_mesh !
	//mmr_create_edge_elems(mesh_id,0);

	/*kbw
	  printf("in create_constr_data: Field_id %d, mesh_id %d\n", Field_id, mesh_id);
	/*kew*/

	/* looking for constrained mid-edge nodes */
	ned = 0;
	while ((ned = mmr_get_next_edge_all(mesh_id, ned)) != 0) {

		if (mmr_edge_status(mesh_id, ned) == MMC_INACTIVE) {

			mmr_edge_sons(mesh_id, ned, edge_sons, NULL);

			/*kbw
				  printf("in parent edge %d, sons %d %d\n",
					 ned, edge_sons[0], edge_sons[1]);
			/*kew*/

			mmr_edge_nodes(mesh_id, edge_sons[0], edge_nodes2);

			/*kbw
				  printf("first son nodes %d %d\n", edge_nodes2[0], edge_nodes2[1]);
			/*kew*/


// if parent edge belongs to at least one active element
			int ned_elems[100];
			mmr_edge_elems(mesh_id, ned, ned_elems);
			int iconstr = 0;
			int iee;
			for (iee = 1; iee <= ned_elems[0]; iee++) {
				if (mmr_el_status(mesh_id, ned_elems[iee]) == MMC_ACTIVE) {
					iconstr = 1;
					break;
				}
			}

			if (iconstr == 1) {
				//if(apr_get_ent_pdeg(Field_id, APC_VERTEX, edge_nodes2[1])==0){

				/* indicate that there are constrained nodes in the mesh */
				field_p->constr = APC_TRUE;

				/* active node with no DOFs vectors */
				nvert = edge_nodes2[1];

				if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {
					free(field_p->dof_ents[nvert].vec_dof_1);
				}
				if (field_p->dof_ents[nvert].vec_dof_2 != NULL) {
					free(field_p->dof_ents[nvert].vec_dof_2);
				}
				if (field_p->dof_ents[nvert].vec_dof_3 != NULL) {
					free(field_p->dof_ents[nvert].vec_dof_3);
				}
				field_p->dof_ents[nvert].vec_dof_1 = NULL;
				field_p->dof_ents[nvert].vec_dof_2 = NULL;
				field_p->dof_ents[nvert].vec_dof_3 = NULL;

#ifdef DEBUG_APM
				if (field_p->dof_ents[nvert].vec_dof_1 != NULL ||
						mmr_node_status(mesh_id, nvert) == MMC_FREE) {

					printf("Error 38784 in create_constr_data. Exiting !\n");

				}
#endif

				field_p->dof_ents[nvert].constr =
						(int *) malloc(3 * sizeof(int));
				field_p->dof_ents[nvert].constr[0] = 2;


				mmr_edge_nodes(mesh_id, ned, edge_nodes);

				/*kbw
					printf("parent nodes %d %d\n", edge_nodes[0], edge_nodes[1]);
				/*kew*/

				field_p->dof_ents[nvert].constr[1] = edge_nodes[0];
				field_p->dof_ents[nvert].constr[2] = edge_nodes[1];

			}
		}
	}

	/* looking for constrained mid-face nodes */
	nfa = 0;
	while ((nfa = mmr_get_next_face_all(mesh_id, nfa)) != 0) {

		if (mmr_fa_type(mesh_id, nfa) == MMC_QUAD &&
				mmr_fa_status(mesh_id, nfa) == MMC_INACTIVE) {

			mmr_fa_fam(mesh_id, nfa, fa_sons, &mid_node);

			int nfa_neig[2];
			mmr_fa_eq_neig(mesh_id, nfa, nfa_neig, NULL, NULL);


			/* // if middle node exists */
			/* if face between two elements and only one divided (inactive) */
			if (mid_node > 0 && nfa_neig[0] != 0 && nfa_neig[1] != 0 &&
					( ( mmr_el_status(mesh_id, nfa_neig[0]) == MMC_ACTIVE &&
									mmr_el_status(mesh_id, nfa_neig[1]) == MMC_INACTIVE ) ||
							( mmr_el_status(mesh_id, nfa_neig[0]) == MMC_INACTIVE &&
									mmr_el_status(mesh_id, nfa_neig[1]) == MMC_ACTIVE ) ) ) {


				/* indicate that there are constrained nodes in the mesh */
				field_p->constr = APC_TRUE;

				/* active node with no DOFs vectors */
				nvert = mid_node;

				if (field_p->dof_ents[nvert].vec_dof_1 != NULL) {
					free(field_p->dof_ents[nvert].vec_dof_1);
				}
				if (field_p->dof_ents[nvert].vec_dof_2 != NULL) {
					free(field_p->dof_ents[nvert].vec_dof_2);
				}
				if (field_p->dof_ents[nvert].vec_dof_3 != NULL) {
					free(field_p->dof_ents[nvert].vec_dof_3);
				}
				field_p->dof_ents[nvert].vec_dof_1 = NULL;
				field_p->dof_ents[nvert].vec_dof_2 = NULL;
				field_p->dof_ents[nvert].vec_dof_3 = NULL;

				/*kbw
				printf("in big inactive face %d, with mid_node %d without dofs\n",
					   nfa, mid_node);
				/*kew*/

				/* indicate that there are constrained nodes in the mesh */
				field_p->constr = APC_TRUE;

				/* active node with no DOFs vectors */
				nvert = mid_node;

#ifdef DEBUG_APM
				if (field_p->dof_ents[nvert].vec_dof_1 != NULL ||
						mmr_node_status(mesh_id, nvert) == MMC_FREE) {

					printf("Error 38884 in create_constr_data. Exiting !\n");

				}
#endif

				field_p->dof_ents[nvert].constr =
						(int *) malloc(5 * sizeof(int));
				field_p->dof_ents[nvert].constr[0] = 4;


				num_nodes = mmr_fa_node_coor(mesh_id, nfa, fa_nodes, NULL);

				for (ino = 1; ino <= num_nodes; ino++) {

					field_p->dof_ents[nvert].constr[ino] = fa_nodes[ino];

				}

			}
		}

	}




#ifdef DEBUG_APM
	/* checking */
	nvert = 0;
	while ((nvert = mmr_get_next_node_all(mesh_id, nvert)) != 0) {

		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE)
		{
			if ( field_p->dof_ents[nvert].vec_dof_1 == NULL &&
					field_p->dof_ents[nvert].constr == NULL ) {

				printf("Node %d without dofs after create_constr_data !\n", nvert);
				exit(-1);

			}
		}
	}
#endif

	// edge_elems data structure is assumed to exist all the time !
	//mmr_delete_edge_elems(mesh_id, 0);

	// set the proper number of dof entities
	int nr_dof_ents = 0;
	int nmvert = mmr_get_max_node_id(mesh_id);
	for (nvert = 1; nvert <= nmvert; nvert++) {
		/* for active unconstrained vertices */
		if (mmr_node_status(mesh_id, nvert) == MMC_ACTIVE &&
				field_p->dof_ents[nvert].vec_dof_1 != NULL) nr_dof_ents++;
	}
	field_p->nr_dof_ents = nr_dof_ents;


	return (1);
}

/*------------------------------------------------------------
  apr_get_constr_data - to return constraints data for a node (dof entity)
------------------------------------------------------------*/
int apr_get_constr_data(/* returns: >=0 - success code, <0 - error code */
		int Field_id,     /* in: approximation field ID  */
		int Node_id,      /* in: id of node (can be vertex or edge) */
		int Node_type,    /* in: type of input node APC_VERTEX or APC_EDGE */

		int * Constr,      /* out: table of constraints data,
			 *      Constr[0] - number of constraints
			 *                  (in linear approximation
			 *                  2 - mid-edge node, 4 - mid-side node)
			 */
		int * Constr_type  /* out: table of constraints element type,
			 *      Constr_type[0] - id of coefficient vector
			 */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int i;
	int Nvert;

	/*++++++++++++++++ executable statements ++++++++++++++++*/
	assert(Constr != NULL);
	assert(Node_type == APC_VERTEX);

	Nvert = Node_id;

	/* select the pointer to the approximation field */
	field_p = apr_select_field(Field_id);
	assert(field_p->dof_ents[Nvert].constr != NULL);
	for (i = 0; i <= field_p->dof_ents[Nvert].constr[0]; i++) {
		Constr[i] = field_p->dof_ents[Nvert].constr[i];
		//printf("Node %d constr=%d\n",Nvert,mesh->node[Nvert].constr[i]);
	}

	return (1);
}

/*------------------------------------------------------------
  apr_get_el_constr_data - to return the number and the list of element's real
						   nodes, with the corresponding constraint coefficients
------------------------------------------------------------*/
extern int apr_get_el_constr_data(
		/* returns: >=0 - success code, <0 - error code */
		int Field_id,   /* in: approximation field ID  */
		int El_id,        /* in: element ID */
		int * Nodes,  	 /* out: list of vertex node IDs */
		/* (Nodes[0] - number of nodes, order from mmr_el_node_coor) */
		int * Nr_constr, /* out: list with the numbers of constraints for each vertex*/
		int * Constr_id, /* out: ID's of parent (constraining) nodes */
		int * Constr_type,  /* out: Unused in linerar approximation */
		double * Constr_val /* out: the corresponding constraint coefficients */
)
{

	/* pointer to field structure */
	apt_field * field_p;
	int ino, j, icount, mesh_id;
	int constr[5];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the pointer to the approximation field */
	field_p = apr_select_field(Field_id);

	/* get mesh ID */
	mesh_id = apr_get_mesh_id(Field_id);

	/* get vertices */
	mmr_el_node_coor(mesh_id, El_id, Nodes, NULL);

	icount = 0;
	for (ino = 1; ino <= Nodes[0]; ino++) {

		if (field_p->dof_ents[Nodes[ino]].constr == NULL) {

			Nr_constr[ino] = 1;
			Constr_id[icount] = Nodes[ino];
			Constr_val[icount] = 1.0;
			icount++;

		}
		else {

			double value;

			apr_get_constr_data(Field_id, Nodes[ino], APC_VERTEX, constr, NULL);
			Nr_constr[ino] = constr[0];
			if (constr[0] == 2) value = 0.5;
			if (constr[0] == 4) value = 0.25;
			for (j = 1; j <= constr[0]; j++) {
				Constr_id[icount] = constr[j];
				Constr_val[icount] = value;
				icount++;
			}
		}

	}

	return (icount);
}

/*---------------------------------------------------------
apr_shape_fun_3D_std_lin_prism - to compute values of shape functions and
their local derivatives at a point within the master 3D prismatic linear element
----------------------------------------------------------*/
int apr_shape_fun_3D_std_lin_prism(
		/* returns: the number of shape functions (<=0 - failure)*/
		double * Eta,	  /* in: local coord of the considered point */
		double * Base_phi, /* out: basis functions */
		double * Base_dphix, /* out: x derivative of basis functions */
		double * Base_dphiy, /* out: y derivative of basis functions */
		double * Base_dphiz /* out: z derivative of basis functions */
)
{


	int num_shap;

	Base_phi[0] = (-Eta[2] + 1) / 2. * (1 - Eta[0] - Eta[1]);
	Base_phi[1] = (-Eta[2] + 1) / 2. * Eta[0];
	Base_phi[2] = (-Eta[2] + 1) / 2. * Eta[1];
	Base_phi[3] = (Eta[2] + 1) / 2. * (1 - Eta[0] - Eta[1]);
	Base_phi[4] = (Eta[2] + 1) / 2. * Eta[0];
	Base_phi[5] = (Eta[2] + 1) / 2. * Eta[1];


	if (Base_dphix != NULL)
	{
		Base_dphix[0] = (Eta[2] - 1) / 2.;
		Base_dphix[1] = (-Eta[2] + 1) / 2.;
		Base_dphix[2] = 0.0;
		Base_dphix[3] = (-Eta[2] - 1) / 2.;
		Base_dphix[4] = (Eta[2] + 1) / 2.;
		Base_dphix[5] = 0.0;

		Base_dphiy[0] = (Eta[2] - 1) / 2.;
		Base_dphiy[1] = 0.0;
		Base_dphiy[2] = (-Eta[2] + 1) / 2.;
		Base_dphiy[3] = (-Eta[2] - 1) / 2.;
		Base_dphiy[4] = 0.0;
		Base_dphiy[5] = (Eta[2] + 1) / 2.;

		Base_dphiz[0] = -0.5 * (1 - Eta[0] - Eta[1]);
		Base_dphiz[1] = -0.5 * Eta[0];
		Base_dphiz[2] = -0.5 * Eta[1];
		Base_dphiz[3] = 0.5 * (1 - Eta[0] - Eta[1]);
		Base_dphiz[4] = 0.5 * Eta[0];
		Base_dphiz[5] = 0.5 * Eta[1];
	}

	num_shap = 6;

	return (num_shap);
}


/*------------------------------------------------------------------
apr_elem_calc_3D_prism - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D prismatic elements
-------------------------------------------------------------------*/
double apr_elem_calc_3D_prism(
		/* returns: Jacobian determinant at a point, either for */
		/* 	volume integration if Vec_norm==NULL,  */
		/* 	or for surface integration otherwise */
		int Control,	    /* in: control parameter (what to compute): */
		/*	1  - shape functions and values */
		/*	2  - derivatives and jacobian */
		/* 	>2 - computations on the (Control-2)-th */
		/*	     element's face */
		int Nreq,	    /* in: number of equations */
		int * Pdeg_vec,     /* in: element degree of polynomial */
		int Base_type,	    /* in: type of basis functions: */
		/* REMARK: type of basis functions differentiates element types as well
			examples for standard linear approximation (from include/aph_intf.h):
			  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements
			  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements  */
		double * Eta,	   /* in: local coordinates of the input point */
		double * Node_coor, /* in: array of coordinates of vertices of element */
		double * Sol_dofs,  /* in: array of element' dofs */
		double * Base_phi,  /* out: basis functions */
		double * Base_dphix, /* out: x-derivatives of basis functions */
		double * Base_dphiy, /* out: y-derivatives of basis functions */
		double * Base_dphiz, /* out: z-derivatives of basis functions */
		double * Xcoor,	   /* out: global coordinates of the point*/
		double * Sol,       /* out: solution at the point */
		double * Dsolx,     /* out: derivatives of solution at the point */
		double * Dsoly,     /* out: derivatives of solution at the point */
		double * Dsolz,     /* out: derivatives of solution at the point */
		double * Vec_nor    /* out: outward unit vector normal to the face */
)
{
	/* local variables */
	int nrgeo;		/* number of element geometry dofs */
	int num_shap;		/* number of shape functions = element solution dofs */
	double det, dxdeta[9], detadx[9]; /* jacobian, jacobian matrix
	and its inverse; for geometrical transformation */

	double ds[3]; /* coordinates of vector normal to a face */

	double geo_phi[8];  	/* geometry shape functions */
	double geo_dphix[8]; 	/* derivatives of geometry shape functions */
	double geo_dphiy[8];  /* derivatives of geometry shape functions */
	double geo_dphiz[8];  /* derivatives of geometry shape functions */

	/* auxiliary variables */
	int i, ieq;
	double daux, faux;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* check parameter Nreq */
#ifdef DEBUG_APM
	if (Nreq > APC_MAXEQ) {
		printf("Number of equations %d greater than the limit %d.\n", Nreq, APC_MAXEQ);
		printf("Change APC_MAXEQ in aph_std_lin_prism.h and recompile the code.\n");
		exit(-1);
	}
#endif


	/* set the number of geometry dofs */
	nrgeo = 6;

	/* get the values of shape functions and their LOCAL derivatives */
	if (Base_phi != NULL) {
		if (Control == 1) {
			num_shap = apr_shape_fun_3D_std_lin_prism(Eta, Base_phi, NULL, NULL, NULL);
		}
		else
			num_shap = apr_shape_fun_3D_std_lin_prism(Eta,
							Base_phi, Base_dphix, Base_dphiy, Base_dphiz);

		/* check parameter num_shap */
#ifdef DEBUG_APM
		if (num_shap != 6) {
			printf("Wrong num_shap from apr_shape_fun_3D in linear approximation %d\n"
					, num_shap);
			exit(-1);
		}
#endif
	}

	num_shap = 6;

	geo_phi[0] = (1.0 - Eta[0] - Eta[1]) * (1.0 - Eta[2]) / 2.0;
	geo_phi[1] = Eta[0] * (1.0 - Eta[2]) / 2.0;
	geo_phi[2] = Eta[1] * (1.0 - Eta[2]) / 2.0;
	geo_phi[3] = (1.0 - Eta[0] - Eta[1]) * (1.0 + Eta[2]) / 2.0;
	geo_phi[4] = Eta[0] * (1.0 + Eta[2]) / 2.0;
	geo_phi[5] = Eta[1] * (1.0 + Eta[2]) / 2.0;

	/*physical coordinates*/
	if (Xcoor != NULL) {

#ifdef DEBUG_APM
		if (Node_coor == NULL) {
			printf("Error 23769483 in shape_fun. Exiting!\n");
			exit(-1);
		}
#endif

		Xcoor[0] = 0;
		Xcoor[1] = 0;
		Xcoor[2] = 0;
		for (i = 0; i < nrgeo; i++) {
			Xcoor[0] += Node_coor[3 * i] * geo_phi[i];
			Xcoor[1] += Node_coor[3 * i + 1] * geo_phi[i];
			Xcoor[2] += Node_coor[3 * i + 2] * geo_phi[i];
		}
	}

	/* function value */
	if (Sol != NULL) {

#ifdef DEBUG_APM
		if (Sol_dofs == NULL) {
			printf("Error 23945483 in shape_fun. Exiting!\n");
			exit(-1);
		}
#endif

		for (ieq = 0; ieq < Nreq; ieq++) Sol[ieq] = 0.0;
		for (i = 0; i < num_shap; i++) {
			for (ieq = 0; ieq < Nreq; ieq++) {
				Sol[ieq] += Sol_dofs[i * Nreq + ieq] * Base_phi[i];
			}
		}

		assert(Sol[0] > -10e200 && Sol[0] < 10e200);

	}

	/*if no computations involving derivatives of shape functions*/
	if (Control == 1) return (0.0);

	/* local derivatives of geometrical shape functions*/

	geo_dphix[0] = -(1.0 - Eta[2]) / 2.0;
	geo_dphix[1] =  (1.0 - Eta[2]) / 2.0;
	geo_dphix[2] =  0.0;
	geo_dphix[3] = -(1.0 + Eta[2]) / 2.0;
	geo_dphix[4] =  (1.0 + Eta[2]) / 2.0;
	geo_dphix[5] =  0.0;
	geo_dphiy[0] = -(1.0 - Eta[2]) / 2.0;
	geo_dphiy[1] =  0.0;
	geo_dphiy[2] =  (1.0 - Eta[2]) / 2.0;
	geo_dphiy[3] = -(1.0 + Eta[2]) / 2.0;
	geo_dphiy[4] =  0.0;
	geo_dphiy[5] =  (1.0 + Eta[2]) / 2.0;
	geo_dphiz[0] = -(1.0 - Eta[0] - Eta[1]) / 2.0;
	geo_dphiz[1] = -Eta[0] / 2.0;
	geo_dphiz[2] = -Eta[1] / 2.0;
	geo_dphiz[3] =  (1.0 - Eta[0] - Eta[1]) / 2.0;
	geo_dphiz[4] =  Eta[0] / 2.0;
	geo_dphiz[5] =  Eta[1] / 2.0;

	/* Jacobian matrix J */
	dxdeta[0] = 0.0;
	dxdeta[1] = 0.0;
	dxdeta[2] = 0.0;
	dxdeta[3] = 0.0;
	dxdeta[4] = 0.0;
	dxdeta[5] = 0.0;
	dxdeta[6] = 0.0;
	dxdeta[7] = 0.0;
	dxdeta[8] = 0.0;
	for (i = 0; i < nrgeo; i++) {
		dxdeta[0] += Node_coor[3 * i]  * geo_dphix[i];
		dxdeta[1] += Node_coor[3 * i]  * geo_dphiy[i];
		dxdeta[2] += Node_coor[3 * i]  * geo_dphiz[i];
		dxdeta[3] += Node_coor[3 * i + 1] * geo_dphix[i];
		dxdeta[4] += Node_coor[3 * i + 1] * geo_dphiy[i];
		dxdeta[5] += Node_coor[3 * i + 1] * geo_dphiz[i];
		dxdeta[6] += Node_coor[3 * i + 2] * geo_dphix[i];
		dxdeta[7] += Node_coor[3 * i + 2] * geo_dphiy[i];
		dxdeta[8] += Node_coor[3 * i + 2] * geo_dphiz[i];
	}

	/* Jacobian |J| and inverse of the Jacobian matrix*/
	det = ut_mat3_inv(dxdeta, detadx);

	/* global derivatives of geometrical shape functions - not used */
	/* for(i=0;i<nrgeo;i++){ */
	/*   daux = geo_dphix[i]*detadx[0] */
	/*        + geo_dphiy[i]*detadx[3] */
	/*        + geo_dphiz[i]*detadx[6]; */
	/*   faux = geo_dphix[i]*detadx[1] */
	/*        + geo_dphiy[i]*detadx[4] */
	/*        + geo_dphiz[i]*detadx[7]; */
	/*   geo_dphiz[i] = geo_dphix[i]*detadx[2] */
	/*        + geo_dphiy[i]*detadx[5] */
	/*        + geo_dphiz[i]*detadx[8]; */
	/*   geo_dphix[i] = daux; */
	/*   geo_dphiy[i] = faux; */
	/* } */

	/* global derivatives of solution shape functions */
	if (Base_dphix != NULL) {

#ifdef DEBUG_APM
		if (Base_dphiy == NULL && Base_dphiz == NULL) {
			printf("Error 23983483 in shape_fun. Exiting!\n");
			exit(-1);
		}
#endif

		for (i = 0; i < num_shap; i++) {
			daux = Base_dphix[i] * detadx[0]
					+   Base_dphiy[i] * detadx[3]
					+   Base_dphiz[i] * detadx[6];
			faux = Base_dphix[i] * detadx[1]
					+   Base_dphiy[i] * detadx[4]
					+   Base_dphiz[i] * detadx[7];
			Base_dphiz[i] =   Base_dphix[i] * detadx[2]
					+   Base_dphiy[i] * detadx[5]
					+   Base_dphiz[i] * detadx[8];
			Base_dphix[i] = daux;
			Base_dphiy[i] = faux;
		}
	}

	/* global derivatives of solution */
	if (Dsolx != NULL) {

#ifdef DEBUG_APM
		if (Dsoly == NULL && Dsolz == NULL) {
			printf("Error 239828483 in shape_fun. Exiting!\n");
			exit(-1);
		}
#endif

		assert(Dsoly != NULL);
		assert(Dsolz != NULL);

		for (ieq = 0; ieq < Nreq; ieq++) {
			Dsolx[ieq] = 0.0;
			Dsoly[ieq] = 0.0;
			Dsolz[ieq] = 0.0;
		}
		for (i = 0; i < num_shap; i++) {
			for (ieq = 0; ieq < Nreq; ieq++) {
				Dsolx[ieq] += Sol_dofs[i * Nreq + ieq] * Base_dphix[i];
				Dsoly[ieq] += Sol_dofs[i * Nreq + ieq] * Base_dphiy[i];
				Dsolz[ieq] += Sol_dofs[i * Nreq + ieq] * Base_dphiz[i];
			}
		}

		assert(Dsolx[0] > -10e200 && Dsolx[0] < 10e200);
		assert(Dsoly[0] > -10e200 && Dsoly[0] < 10e200);
		assert(Dsolz[0] > -10e200 && Dsolz[0] < 10e200);

	}

	switch (Control) {
		case 2:
			return (fabs(det));
			break; /* area element dS = vector normal */
		case 3: {
			ds[0] = - dxdeta[3] * dxdeta[7] + dxdeta[6] * dxdeta[4];
			ds[1] = - dxdeta[6] * dxdeta[1] + dxdeta[0] * dxdeta[7];
			ds[2] = - dxdeta[0] * dxdeta[4] + dxdeta[3] * dxdeta[1];
		}
		break;
		case 4: {
			ds[0] = dxdeta[3] * dxdeta[7] - dxdeta[6] * dxdeta[4];
			ds[1] = dxdeta[6] * dxdeta[1] - dxdeta[0] * dxdeta[7];
			ds[2] = dxdeta[0] * dxdeta[4] - dxdeta[3] * dxdeta[1];
		}
		break;
		case 5: {
			ds[0] = dxdeta[3] * dxdeta[8] - dxdeta[6] * dxdeta[5];
			ds[1] = dxdeta[6] * dxdeta[2] - dxdeta[0] * dxdeta[8];
			ds[2] = dxdeta[0] * dxdeta[5] - dxdeta[3] * dxdeta[2];
		}
		break;
		case 6: {
			ds[0] = (dxdeta[4] - dxdeta[3]) * dxdeta[8]
					- (dxdeta[7] - dxdeta[6]) * dxdeta[5];
			ds[1] = (dxdeta[7] - dxdeta[6]) * dxdeta[2]
					- (dxdeta[1] - dxdeta[0]) * dxdeta[8];
			ds[2] = (dxdeta[1] - dxdeta[0]) * dxdeta[5]
					- (dxdeta[4] - dxdeta[3]) * dxdeta[2];
		}
		break;
		case 7: {
			ds[0] = dxdeta[5] * dxdeta[7] - dxdeta[8] * dxdeta[4];
			ds[1] = dxdeta[8] * dxdeta[1] - dxdeta[2] * dxdeta[7];
			ds[2] = dxdeta[2] * dxdeta[4] - dxdeta[5] * dxdeta[1];
		}
		break;

	} //!switch(Control)


	det = ut_vec3_length(ds);

	if (Vec_nor != NULL) {

		/* normalize vector normal */
		Vec_nor[0] = ds[0] / det;
		Vec_nor[1] = ds[1] / det;
		Vec_nor[2] = ds[2] / det;

	}

	return (det);

}

/*---------------------------------------------------------
apr_shape_fun_3D_std_lin_tetra - to compute values of shape functions and their
local derivatives at a point within the master 3D tetrahedral linear element
----------------------------------------------------------*/
int apr_shape_fun_3D_std_lin_tetra(
		/* returns: the number of shape functions (<=0 - failure)*/
		double * Eta,	  /* in: local coord of the considered point */
		double * Base_phi, /* out: basis functions */
		double * Base_dphix, /* out: x derivative of basis functions */
		double * Base_dphiy, /* out: y derivative of basis functions */
		double * Base_dphiz /* out: z derivative of basis functions */
) {
	const int num_shap = 4;

	assert(Eta[0] >= -APC_SMALL && Eta[0] <= 1.0 + APC_SMALL);
	assert(Eta[1] >= -APC_SMALL && Eta[1] <= 1.0 + APC_SMALL);
	assert(Eta[2] >= -APC_SMALL && Eta[2] <= 1.0 + APC_SMALL);
	assert(Eta[0] + Eta[1] + Eta[2] <= 1.0 + APC_SMALL);

	Base_phi[0] = 1.0 - Eta[0] - Eta[1] - Eta[2];
	Base_phi[1] = Eta[0];
	Base_phi[2] = Eta[1];
	Base_phi[3] = Eta[2];


	if (Base_dphix != NULL) {
		Base_dphix[0] =  -1.0;
		Base_dphix[1] =  1.0;
		Base_dphix[2] =  0.0;
		Base_dphix[3] =  0.0;

		Base_dphiy[0] = -1.0 ;
		Base_dphiy[1] =  0.0;
		Base_dphiy[2] =  1.0;
		Base_dphiy[3] =  0.0;

		Base_dphiz[0] = -1.0;
		Base_dphiz[1] =  0.0;
		Base_dphiz[2] =  0.0;
		Base_dphiz[3] =  1.0;


	}

	return (num_shap);
}

double apr_elem_calc_3D_tetra(
		/* returns: Jacobian determinant at a point, either for */
		/* 	volume integration if Vec_norm==NULL,  */
		/* 	or for surface integration otherwise */
		int Control,	    /* in: control parameter (what to compute): */
		/*	1  - shape functions and values */
		/*	2  - derivatives and jacobian */
		/* 	>2 - computations on the (Control-2)-th */
		/*	     element's face */
		int Nreq,	    /* in: number of equations */
		int * Pdeg_vec,	   /* in: element degree of polynomial */
		int Base_type,	    /* in: type of basis functions: */
		/* REMARK: type of basis functions differentiates element types as well
			examples for standard linear approximation (from include/aph_intf.h):
			  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements
			  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements  */
		double * Eta,	   /* in: local coordinates of the input point */
		double * Node_coor, /* in: array of coordinates of vertices of element */
		double * Sol_dofs,  /* in: array of element' dofs */
		double * Base_phi,  /* out: basis functions */
		double * Base_dphix, /* out: x-derivatives of basis functions */
		double * Base_dphiy, /* out: y-derivatives of basis functions */
		double * Base_dphiz, /* out: z-derivatives of basis functions */
		double * Xcoor,	   /* out: global coordinates of the point*/
		double * Sol,       /* out: solution at the point */
		double * Dsolx,     /* out: derivatives of solution at the point */
		double * Dsoly,     /* out: derivatives of solution at the point */
		double * Dsolz,     /* out: derivatives of solution at the point */
		double * Vec_nor    /* out: outward unit vector normal to the face */
) {
	/* local variables */
	//int nrgeo=4;		/* number of element geometry dofs */
	int num_shap;	/* number of shape functions = element solution dofs */
	double det = 0, dxdeta[9] = {0}, detadx[9] = {0}; /* jacobian, jacobian matrix
						   and its inverse; for geometrical transformation */

	double ds[3] = {0}, v1[3] = {0}, v2[3] = {0};
	/* coordinates of vector normal to a face */

	//    double geo_phi[8];  	/* geometry shape functions */
	//    double geo_dphix[8]; 	/* derivatives of geometry shape functions */
	//    double geo_dphiy[8];  /* derivatives of geometry shape functions */
	//    double geo_dphiz[8];  /* derivatives of geometry shape functions */

	/* auxiliary variables */
	int i, ieq, n0 = 0, n1 = 0, n2 = 0;
	double daux, faux;
	// below is 1.0/6.0 which is fraction of normal volume for tetrahedron
	static const double one_sixth = 0.16666666666666666666666666666667;
	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* check parameter Nreq */
#ifdef DEBUG_APM
	if (Nreq > APC_MAXEQ) {
		printf("Number of equations %d greater than the limit %d.\n", Nreq, APC_MAXEQ);
		printf("Change APC_MAXEQ in aph_std_lin_prism.h and recompile the code.\n");
		exit(-1);
	}
#endif

	/* get the values of shape functions and their LOCAL derivatives */
	if (Base_phi != NULL) {
		if (Control == APC_ELEM_CALC_SHAPE_FUNC_N_VALUES) {
			num_shap = apr_shape_fun_3D_std_lin_tetra(Eta, Base_phi,
							NULL, NULL, NULL);
		} else
			num_shap = apr_shape_fun_3D_std_lin_tetra(Eta, Base_phi,
							Base_dphix, Base_dphiy, Base_dphiz);

		/* check parameter num_shap */
#ifdef DEBUG_APM
		if (num_shap != 4) {
			printf("Wrong num_shap from apr_shape_fun_3D in linear approximation %d\n"
					, num_shap);
			exit(-1);
		}
#endif
	}

	num_shap = 4;
	//geo_phi[0] = (1.0- Eta[0] - Eta[1] - Eta[2])*one_sixth;
	//geo_phi[1] = Eta[0]*one_sixth;
	//geo_phi[2] = Eta[1]*one_sixth;
	//geo_phi[3] = Eta[2]*one_sixth;

	//assert(geo_phi[0] >= -SMALL && geo_phi[0]<=1.0+SMALL);
	//assert(geo_phi[1] >= -SMALL && geo_phi[1]<=1.0+SMALL);
	//assert(geo_phi[2] >= -SMALL && geo_phi[2]<=1.0+SMALL);
	//assert(geo_phi[3] >= -SMALL && geo_phi[3]<=1.0+SMALL);


#ifdef DEBUG_APM
	if (Node_coor == NULL) {
		printf("Error 23769483 in shape_fun. Exiting!\n");
		exit(-1);
	}
#endif


	dxdeta[0] = Node_coor[3 + 0] - Node_coor[0 + 0];
	dxdeta[3] = Node_coor[3 + 1] - Node_coor[0 + 1];
	dxdeta[6] = Node_coor[3 + 2] - Node_coor[0 + 2];

	dxdeta[1] = Node_coor[6 + 0] - Node_coor[0 + 0];
	dxdeta[4] = Node_coor[6 + 1] - Node_coor[0 + 1];
	dxdeta[7] = Node_coor[6 + 2] - Node_coor[0 + 2];

	dxdeta[2] = Node_coor[9 + 0] - Node_coor[0 + 0];
	dxdeta[5] = Node_coor[9 + 1] - Node_coor[0 + 1];
	dxdeta[8] = Node_coor[9 + 2] - Node_coor[0 + 2];

#ifdef DEBUG_APM
	{
		double v1[3] = {dxdeta[0], dxdeta[3], dxdeta[6]};
		double v2[3] = {dxdeta[1], dxdeta[4], dxdeta[7]};
		double v3[3] = {dxdeta[2], dxdeta[5], dxdeta[8]};
		double hsize = pow(fabs(utr_vec3_mxpr(v1, v2, v3)), 1.0 / 3.0);
		if (hsize < 0.0001) {
			printf("degenerate element in elem_calc: hsize %.12lf\n", hsize);
			printf("v1: %lf, %lf, %lf\n", v1[0], v1[1], v1[2]);
			printf("v2: %lf, %lf, %lf\n", v2[0], v2[1], v2[2]);
			printf("v3: %lf, %lf, %lf\n", v3[0], v3[1], v3[2]);
			exit(-1);
		}
	}
#endif

	/*physical coordinates*/
	if (Xcoor != NULL) {

		Xcoor[0] = Node_coor[0] + Eta[0] * dxdeta[0] + Eta[1] * dxdeta[1] + Eta[2] * dxdeta[2];
		Xcoor[1] = Node_coor[1] + Eta[0] * dxdeta[3] + Eta[1] * dxdeta[4] + Eta[2] * dxdeta[5];
		Xcoor[2] = Node_coor[2] + Eta[0] * dxdeta[6] + Eta[1] * dxdeta[7] + Eta[2] * dxdeta[8];

		//Xcoor[0] = 0;
		//Xcoor[1] = 0;
		//Xcoor[2] = 0;
		//for (i=0;i<nrgeo;i++) {
		//    Xcoor[0] += Node_coor[3*i]*geo_phi[i];
		//    Xcoor[1] += Node_coor[3*i+1]*geo_phi[i];
		//    Xcoor[2] += Node_coor[3*i+ 2]*geo_phi[i];
		//}
	}

	/* function value */
	if (Sol != NULL) {

#ifdef DEBUG_APM
		if (Sol_dofs == NULL) {
			printf("Error 23945483 in shape_fun. Exiting!\n");
			exit(-1);
		}
#endif

		for (ieq = 0; ieq < Nreq; ieq++) Sol[ieq] = 0.0;
		for (i = 0; i < num_shap; i++) {
			for (ieq = 0; ieq < Nreq; ieq++) {
				assert(Sol_dofs[i * Nreq + ieq] > -10e200 && Sol_dofs[i * Nreq + ieq] < 10e200);
				Sol[ieq] += Sol_dofs[i * Nreq + ieq] * Base_phi[i];
			}
		}

		// check if output is reasonable
		assert(Sol[0] > -10e200 && Sol[0] < 10e200);
		assert(Sol[1] > -10e200 && Sol[1] < 10e200);


	}

	/*if no computations involving derivatives of shape functions*/
	if (Control == APC_ELEM_CALC_SHAPE_FUNC_N_VALUES) return (0.0);

	// det = utr_vec3_mxpr(v1,v2,v3)*one_sixth;

	/* local derivatives of geometrical shape functions*/

	/*
	  geo_dphix[0] = -1.0;	// -1.0/6.0;
	  geo_dphix[1] =  1.0;	//1.0/6.0;
	  geo_dphix[2] =  0.0;	//0.0/6.0;
	  geo_dphix[3] =  0.0;	//0.0/6.0;

	  geo_dphiy[0] =  -1.0;	//-1.0/6.0 ;
	  geo_dphiy[1] =  0.0;	//0.0/6.0;
	  geo_dphiy[2] =  1.0;	//1.0/6.0;
	  geo_dphiy[3] =  0.0;	//0.0/6.0;


	  geo_dphiz[0] = -1.0;	//-1.0/6.0;
	  geo_dphiz[1] =  0.0;	//0.0/6.0;
	  geo_dphiz[2] =  0.0;	//0.0/6.0;
	  geo_dphiz[3] =  1.0;	//1.0/6.0;

	*/

	/* Jacobian matrix J */
	/*
		dxdeta[0] = 0.0; dxdeta[1] = 0.0; dxdeta[2] = 0.0;
		dxdeta[3] = 0.0; dxdeta[4] = 0.0; dxdeta[5] = 0.0;
		dxdeta[6] = 0.0; dxdeta[7] = 0.0; dxdeta[8] = 0.0;
		for (i=0;i<nrgeo;i++) {
			dxdeta[0] += Node_coor[3*i]  *geo_dphix[i];
			dxdeta[1] += Node_coor[3*i]  *geo_dphiy[i];
			dxdeta[2] += Node_coor[3*i]  *geo_dphiz[i];
			dxdeta[3] += Node_coor[3*i+1]*geo_dphix[i];
			dxdeta[4] += Node_coor[3*i+1]*geo_dphiy[i];
			dxdeta[5] += Node_coor[3*i+1]*geo_dphiz[i];
			dxdeta[6] += Node_coor[3*i+2]*geo_dphix[i];
			dxdeta[7] += Node_coor[3*i+2]*geo_dphiy[i];
			dxdeta[8] += Node_coor[3*i+2]*geo_dphiz[i];
		}
	*/

	/* Jacobian |J| and inverse of the Jacobian matrix*/
	det = ut_mat3_inv(dxdeta, detadx) * one_sixth; // *1/6 fot tetrahedrons

	// global derivatives of geometrical shape functions
	/* for (i=0;i<nrgeo;i++) { */
	/*     daux = geo_dphix[i]*detadx[0] */
	/*            + geo_dphiy[i]*detadx[3] */
	/*            + geo_dphiz[i]*detadx[6]; */
	/*     faux = geo_dphix[i]*detadx[1] */
	/*            + geo_dphiy[i]*detadx[4] */
	/*            + geo_dphiz[i]*detadx[7]; */
	/*     geo_dphiz[i] = geo_dphix[i]*detadx[2] */
	/*                    + geo_dphiy[i]*detadx[5] */
	/*                    + geo_dphiz[i]*detadx[8]; */
	/*     geo_dphix[i] = daux; */
	/*     geo_dphiy[i] = faux; */
	/* } */

	/* global derivatives of solution shape functions */
	if (Base_dphix != NULL) {

#ifdef DEBUG_APM
		if (Base_dphiy == NULL && Base_dphiz == NULL) {
			printf("Error 23983483 in shape_fun. Exiting!\n");
			exit(-1);
		}
#endif

		for (i = 0; i < num_shap; i++) {
			daux = Base_dphix[i] * detadx[0]
					+   Base_dphiy[i] * detadx[3]
					+   Base_dphiz[i] * detadx[6];
			faux = Base_dphix[i] * detadx[1]
					+   Base_dphiy[i] * detadx[4]
					+   Base_dphiz[i] * detadx[7];
			Base_dphiz[i] =   Base_dphix[i] * detadx[2]
					+   Base_dphiy[i] * detadx[5]
					+   Base_dphiz[i] * detadx[8];
			Base_dphix[i] = daux;
			Base_dphiy[i] = faux;
		}
	}

	/* global derivatives of solution */
	if (Dsolx != NULL) {

#ifdef DEBUG_APM
		if (Dsoly == NULL && Dsolz == NULL) {
			printf("Error 239828483 in shape_fun. Exiting!\n");
			exit(-1);
		}
#endif

		assert(Dsoly != NULL);
		assert(Dsolz != NULL);

		for (ieq = 0; ieq < Nreq; ieq++) {
			Dsolx[ieq] = 0.0;
			Dsoly[ieq] = 0.0;
			Dsolz[ieq] = 0.0;
		}
		for (i = 0; i < num_shap; i++) {
			for (ieq = 0; ieq < Nreq; ieq++) {
				Dsolx[ieq] += Sol_dofs[i * Nreq + ieq] * Base_dphix[i];
				Dsoly[ieq] += Sol_dofs[i * Nreq + ieq] * Base_dphiy[i];
				Dsolz[ieq] += Sol_dofs[i * Nreq + ieq] * Base_dphiz[i];
			}
		}

		assert(Dsolx[0] > -10e200 && Dsolx[0] < 10e200);
		assert(Dsoly[0] > -10e200 && Dsoly[0] < 10e200);
		assert(Dsolz[0] > -10e200 && Dsolz[0] < 10e200);

	}





	switch (Control) {
		case APC_ELEM_CALC_DERIVATIVES_AND_JACOBIAN:
			return (fabs(det));
			break;
		/* area element dS = vector normal, basing at Minors of J matrix */
		case APC_ELEM_CALC_ON_FACE(0): { // Face 0 - base (pointing inside elem, so swapping v1 and v2
			n0 = MMC_FACE_NODES_FOR_TETRA[0][0] * 3;
			n1 = MMC_FACE_NODES_FOR_TETRA[0][1] * 3;
			n2 = MMC_FACE_NODES_FOR_TETRA[0][2] * 3;
		}
		break;
		case APC_ELEM_CALC_ON_FACE(1): { // Face 1 - at 0-1 edge - pointing outside
			n0 = MMC_FACE_NODES_FOR_TETRA[1][0] * 3;
			n1 = MMC_FACE_NODES_FOR_TETRA[1][1] * 3;
			n2 = MMC_FACE_NODES_FOR_TETRA[1][2] * 3;
		}
		break;
		case APC_ELEM_CALC_ON_FACE(2): { // Face 2 - at 0-2 edge - pointing inside, so swap
			n0 = MMC_FACE_NODES_FOR_TETRA[2][0] * 3;
			n1 = MMC_FACE_NODES_FOR_TETRA[2][1] * 3;
			n2 = MMC_FACE_NODES_FOR_TETRA[2][2] * 3;
		}
		break;
		case APC_ELEM_CALC_ON_FACE(3): { // Face 3 - at 1-2 edge - pointing outside
			n0 = MMC_FACE_NODES_FOR_TETRA[3][0] * 3;
			n1 = MMC_FACE_NODES_FOR_TETRA[3][1] * 3;
			n2 = MMC_FACE_NODES_FOR_TETRA[3][2] * 3;
		}
		break;
		default: {
			printf("Error: Apr_elem_calc_3D_tetra: bad control value!");
			exit(-1);
		}
		break;
	} //!switch

	v1[0] = Node_coor[n1 + 0] - Node_coor[n0 + 0];
	v1[1] = Node_coor[n1 + 1] - Node_coor[n0 + 1];
	v1[2] = Node_coor[n1 + 2] - Node_coor[n0 + 2];

	v2[0] = Node_coor[n2 + 0] - Node_coor[n1 + 0];
	v2[1] = Node_coor[n2 + 1] - Node_coor[n1 + 1];
	v2[2] = Node_coor[n2 + 2] - Node_coor[n1 + 2];

	utr_vec3_prod(v1, v2, ds);

	det = ut_vec3_length(ds);

	if (Vec_nor != NULL) {

		/* normalize vector normal */
		Vec_nor[0] = ds[0] / det;
		Vec_nor[1] = ds[1] / det;
		Vec_nor[2] = ds[2] / det;

	}


	return (det);

}


/*---------------------------------------------------------
apr_shape_fun_3D - to compute values of shape functions and their
	local derivatives at a point within the master 3D element
----------------------------------------------------------*/
int apr_shape_fun_3D( /* returns: the number of shape functions (<=0 - failure) */
		int Base_type,	   /* in: type of basis functions: */
		/* REMARK: type of basis functions differentiates element types as well */
		/*  #define APC_BASE_TENSOR_DG 1 - tensor product for prismatic elements and DG*/
		/*  #define APC_BASE_COMPLETE_DG 2 - complete polynomials for prismatic elements and DG*/
		/*  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements */
		/*  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements  */
		int Pdeg, 	   /* in: degree of polynomial - can be either */
		/*	a single number, for isotropic p, */
		/*	or a combination pdegy*10+pdegx */
		double * Eta,	  /* in: local coord of the considered point */
		double * Base_phi, /* out: basis functions */
		double * Base_dphix, /* out: x derivative of basis functions */
		double * Base_dphiy, /* out: y derivative of basis functions */
		double * Base_dphiz /* out: z derivative of basis functions */
)
{

	int num_shap;

	if (Base_type == APC_BASE_TETRA_STD) {

		if (Base_phi != NULL) {
			if (Base_dphix == NULL) {
				num_shap = apr_shape_fun_3D_std_lin_tetra(Eta, Base_phi,
								NULL, NULL, NULL);
			} else {
				num_shap = apr_shape_fun_3D_std_lin_tetra(Eta, Base_phi,
								Base_dphix, Base_dphiy, Base_dphiz);
			}
		}

	} else {

		if (Base_phi != NULL) {
			if (Base_dphix == NULL) {
				num_shap = apr_shape_fun_3D_std_lin_prism(Eta, Base_phi, NULL, NULL, NULL);
			}
			else {
				num_shap = apr_shape_fun_3D_std_lin_prism(Eta,
								Base_phi, Base_dphix, Base_dphiy, Base_dphiz);
			}
		}

	}

	return num_shap;
}

/*------------------------------------------------------------------
apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
-------------------------------------------------------------------*/
double apr_elem_calc_3D(
		/* returns: Jacobian determinant at a point, either for */
		/* 	volume integration if Vec_norm==NULL,  */
		/* 	or for surface integration otherwise */
		int Control,	    /* in: control parameter (what to compute): */
		/*	1  - shape functions and values */
		/*	2  - derivatives and jacobian */
		/* 	>2 - computations on the (Control-2)-th */
		/*	     element's face */
		int Nreq,	    /* in: number of equations */
		int * Pdeg_vec, /* in: element degree of polynomial */
		int Base_type,	    /* in: type of basis functions: */
		/* REMARK: type of basis functions differentiates element types as well
			examples for standard linear approximation (from include/aph_intf.h):
			  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements
			  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements  */
		double * Eta,	   /* in: local coordinates of the input point */
		double * Node_coor, /* in: array of coordinates of vertices of element */
		double * Sol_dofs,  /* in: array of element' dofs */
		double * Base_phi,  /* out: basis functions */
		double * Base_dphix, /* out: x-derivatives of basis functions */
		double * Base_dphiy, /* out: y-derivatives of basis functions */
		double * Base_dphiz, /* out: z-derivatives of basis functions */
		double * Xcoor,	   /* out: global coordinates of the point*/
		double * Sol,       /* out: solution at the point */
		double * Dsolx,     /* out: derivatives of solution at the point */
		double * Dsoly,     /* out: derivatives of solution at the point */
		double * Dsolz,     /* out: derivatives of solution at the point */
		double * Vec_nor   /* out: outward unit vector normal to the face */
) {
	double determ = 0.0;


	if (Base_type == APC_BASE_TETRA_STD) {
		determ = apr_elem_calc_3D_tetra(Control, Nreq, Pdeg_vec, Base_type,
						Eta, Node_coor, Sol_dofs,
						Base_phi, Base_dphix, Base_dphiy, Base_dphiz,
						Xcoor, Sol, Dsolx, Dsoly, Dsolz, Vec_nor);
	} else {
		determ = apr_elem_calc_3D_prism(Control, Nreq, Pdeg_vec, Base_type,
						Eta, Node_coor, Sol_dofs,
						Base_phi, Base_dphix, Base_dphiy, Base_dphiz,
						Xcoor, Sol, Dsolx, Dsoly, Dsolz, Vec_nor);
	}

	return determ;
}


/*---------------------------------------------------------
apr_set_quadr_3D - to prepare quadrature data for a given element
---------------------------------------------------------*/
int apr_set_quadr_3D(
		int Base_type,	   /* in: type of basis functions: */
		/*
		 REMARK: type of basis functions differentiates element types as well
			examples for standard linear approximation (from include/aph_intf.h):
			  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements
			  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements
		*/
		int * Pdeg_vec,	/* in: element degree of polynomial */
		int * Ngauss,	/* out: number of gaussian points */
		double * Xg,	/* out: coordinates of gaussian points */
		double * Wg	/* out: weights associated with points */
)
{

	/* local variables */
	int pdeg;
	int orderx, orderz; 	/* orders of approximation in x,y,z */
	int ngauss2, ngaussz;  /* numbers of gaussian points */
	double * xg2, *xgz; /* gauss points in 2D and 1D*/
	double * wg2, *wgz; /* gauss weights in 2D and 1D*/
	int ki, kj; //iaux,

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* currently only uniform p elements are considered */
	pdeg = *Pdeg_vec;

	/* currently second order integration is assumed for linear elements */
	/* it is possible to change it by modifying this routine */


	if (Base_type == APC_BASE_TETRA_STD) {
		if (pdeg == 0) {
			orderz = 0;
		} else if (pdeg == 101) {
			orderz = 2;
		} else {
			orderz = 2;
		}

		apr_gauss_select(3, orderz, &ngaussz, &xgz, &wgz);
		*Ngauss = ngaussz;
		/* fill the arrays of Gauss points and Gauss coefficients */
		for (ki = 0; ki < ngaussz; ki++) {
			Xg[3 * ki]   = xgz[3 * ki];
			Xg[3 * ki + 1] = xgz[3 * ki + 1];
			Xg[3 * ki + 2] = xgz[3 * ki + 2];
			Wg[ki]     = wgz[ki];
		}
	} else { // prism elements
		if (pdeg == 0) {
			orderz = 0;
			orderx = 0;
		} else if (pdeg == 101) {
			orderz = 2;
			orderx = 2;
		} else {
			orderz = 2;
			orderx = 2;
		}

		//iaux=1; /* 1D integration for vertical direction */
		apr_gauss_select(1, orderz, &ngaussz, &xgz, &wgz);
		//iaux=2; /* 2D integration for triangular bases */
		apr_gauss_select(2, orderx, &ngauss2, &xg2, &wg2);

		/*kbw
		  printf("selecting gauss for orderz %d and orderx %d\n",
			 orderz, orderx);
		  for(ki=0;ki<ngaussz;ki++){
			printf("point %d, weight %lf, coor %lf\n",
			   ki, wgz[ki], xgz[ki]);
		  }
		  for(ki=0;ki<ngauss2;ki++){
			printf("point %d, weight %lf, coor %lf %lf\n",
			   ki, wg2[ki], xg2[3*ki+1], xg2[3*ki+2]);
		  }

		/*kew*/


		/* set the total number of Gauss points */
		*Ngauss = ngaussz * ngauss2;

		/* fill the arrays of Gauss points and Gauss coefficients */
		for (ki = 0; ki < ngaussz; ki++) {
			for (kj = 0; kj < ngauss2; kj++) {

				Xg[3 * (ki * ngauss2 + kj)] = xg2[3 * kj + 1];
				Xg[3 * (ki * ngauss2 + kj) + 1] = xg2[3 * kj + 2];
				Xg[3 * (ki * ngauss2 + kj) + 2] = xgz[ki];
				/* we correct weights for 2D integration on triangles */
				Wg[ki * ngauss2 + kj] = 0.5 * wg2[kj] * wgz[ki];

			}
		}
	}
	return (1);
}

/*---------------------------------------------------------
apr_set_quadr_2D - to prepare quadrature data for a given face
---------------------------------------------------------*/
int apr_set_quadr_2D(
		int Fa_type,	/* in: type of a face */
		int Base_type,	   /* in: type of basis functions: */
		/*   (APC_TENSOR) - tensor product */
		/* 	(APC_COMPLETE) - complete polynomials */
		int * Pdeg_vec,	/* in: element degree of polynomial */
		int * Ngauss,	/* out: number of gaussian points */
		double * Xg,	/* out: coordinates of gaussian points */
		double * Wg	/* out: weights associated with points */
)
{

	/* local variables */
	int pdeg;
	int orderx, orderz; 	/* orders of approximation in x,y,z */
	int ngauss2, ngaussx, ngaussy; /* numbers of gaussian points */
	double * xg2, *xgx, *xgy; /* gauss points in 2D and 1D*/
	double * wg2, *wgx, *wgy; /* gauss weights in 2D and 1D*/
	int iaux, ki, kj;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* currently only uniform p elements are considered */
	pdeg = *Pdeg_vec;

	/* currently second order integration is assumed for linear elements */
	/* it is possible to change it by modifying this routine */
	if (pdeg == 0) {
		orderz = 0;
		orderx = 0;
	} else if (pdeg == 101) {
		orderz = 2;
		orderx = 2;
	} else {
		orderz = 2;
		orderx = 2;
	}

	/* fill the arrays of Gauss points and Gauss coefficients */

	if (Fa_type == MMC_QUAD) {

		iaux = 1; /* 1D line integration */
		apr_gauss_select(iaux, orderx, &ngaussx, &xgx, &wgx);
		apr_gauss_select(iaux, orderz, &ngaussy, &xgy, &wgy);
		*Ngauss = ngaussx * ngaussy;

		if (Xg != NULL) {

			for (ki = 0; ki < ngaussx; ki++) {
				for (kj = 0; kj < ngaussy; kj++) {

					Xg[2 * (ki * ngaussy + kj)] = xgx[ki];
					Xg[2 * (ki * ngaussy + kj) + 1] = xgy[kj];
					/* we correct weights because 2D master element side is from 0 to 1 */
					/* but Xg is unchanged so mmr_fa_elem_coor must take care of this */
					Wg[ki * ngaussy + kj] = 0.5 * wgx[ki] * wgy[kj];

				}
			}

		}
	}
	else if (Fa_type == MMC_TRIA) {

		iaux = 2; /* 2D integration for triangles */
		apr_gauss_select(iaux, orderx, &ngauss2, &xg2, &wg2);
		*Ngauss = ngauss2;

		if (Xg != NULL) {

			for (ki = 0; ki < ngauss2; ki++) {
				Xg[2 * ki] = xg2[3 * ki + 1];
				Xg[2 * ki + 1] = xg2[3 * ki + 2];
				/* we correct weights for 2D integration on triangles */
				Wg[ki] = 0.5 * wg2[ki];
			}

		}
	}

	return (1);
}

/*---------------------------------------------------------
apr_set_quadr_2D_penalty - to prepare penalty quadrature data for a given face
			  in general should be Gauss-Lobatto integration
			  for linear elements just vertex based integration
---------------------------------------------------------*/
int apr_set_quadr_2D_penalty(
		int Fa_type,	/* in: type of a face */
		int Base_type,	   /* in: type of basis functions: */
		/*   (APC_TENSOR) - tensor product */
		/* 	(APC_COMPLETE) - complete polynomials */
		int * Pdeg_vec,	/* in: element degree of polynomial */
		int * Ngauss,	/* out: number of gaussian points */
		double * Xg,	/* out: coordinates of gaussian points */
		double * Wg	/* out: weights associated with points */
)
{

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* simple vertex based integration */
	if (Fa_type == MMC_QUAD) {

		*Ngauss = 4;
		if (Xg != NULL) {
			Xg[0] = -1.0;
			Xg[1] = -1.0;
			/* we correct weights because 2D master element side is from 0 to 1 */
			/* but Xg is unchanged so mmr_fa_elem_coor must take care of this */
			Wg[0] = 0.5;
			Xg[2] = -1.0;
			Xg[3] = 1.0;
			Wg[1] = 0.5;
			Xg[4] = 1.0;
			Xg[5] = 1.0;
			Wg[2] = 0.5;
			Xg[6] = 1.0;
			Xg[7] = -1.0;
			Wg[3] = 0.5;

		}
	}
	else if (Fa_type == MMC_TRIA) {

		*Ngauss = 4;

		if (Xg != NULL) {
			Xg[0] = 0.0;
			Xg[1] = 0.0;
			Wg[0] = 0.166666666666667;
			/* we correct weights for 2D integration on triangles */
			Xg[2] = 0.0;
			Xg[3] = 1.0;
			Wg[1] = 0.166666666666667;
			Xg[4] = 1.0;
			Xg[5] = 0.0;
			Wg[2] = 0.166666666666667;

		}
	}


	return (1);
}


/*---------------------------------------------------------
apr_L2_proj - to L2 project a function onto an element
WARNING: projections from different elements do not form a continuous function!!!
---------------------------------------------------------*/
int apr_L2_proj(	/* returns: >0 - success, <=0 - failure */
		int Field_id,   /* in: field ID */
		int Mode,	/* in: mode of operation */
		/*    <-1 - projection from ancestor to father */
		/*          the value is the number of ancestors */
		/*     -1 - projection from father to son */
		/*     >0 - projection of function, the routine */
		/*          returning value at point is specified */
		/*          by the pointer in the last argument */
		int El,		/* in: element number */
		int * Pdeg_vec,	/* in: element degree of approximation */
		double * Dofs,	/* out: workspace for degress of freedom of El */
		/* 	NULL - write to  data structure */
		int * El_from,	/* in: list of elements to provide function */
		int * Pdeg_vec_from,	/* in: degree of polynomial for each El_from */
		double * Dofs_from, /* in: Dofs of El_from or...*/
		double (*Fun_p)(double *, double *, double *, double *) /* in: pointer to */
		/* function with field values and its derivatives */
)
{

	/* local variables */
	int mesh_id;
	int pdeg;
//  int *pdeg_from;
	int base_q;		/* type of basis functions for quadrilaterals */
	int num_shap;         /* number of element shape functions */
	int ndofs;            /* local dimension of the problem */
	int ngauss;           /* number of gauss points */
	double xg[3000];   	/* coordinates of gauss points in 3D */
	double wg[1000];      /* gauss weights */
	double determ;        /* determinant of jacobi matrix */
	double vol;           /* volume for integration rule */
	double xcoor[3];      /* global coord of gauss point */
	double xloc[3];       /* local point coordinates */
	double base_phi[APC_MAXELVD];    /* basis functions */
	double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
	double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
	double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
	double base_aux[APC_MAXELVD];    /* basis functions */
	int el_nodes[MMC_MAXELVNO + 1];    /* list of nodes of El */
	double node_coor[3 * MMC_MAXELVNO]; /* coord of nodes of El */
	double stiff_loc[APC_MAXELSD * APC_MAXELSD];
	/* stiffness matrix for local problem */
	double f_loc[APC_MAXELSD]; /* rhs vector for local problem */
	int ipiv[APC_MAXELSD];
	double value[APC_MAXEQ]; /* solution vector to be projected */
	double u_val[APC_MAXEQ]; /* projected solution vector */

	/* auxiliary variables */
	int i, ki, kk, iaux, iel, iel_from, ifound, nreq;
	int idofs, jdofs, ieq;
//  double daux, time;

	/* constatnts */
	int ione = 1;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(Field_id);

	/* get formulation parameters */
	base_q = apr_get_base_type(Field_id, El);
	nreq = apr_get_nreq(Field_id);

	/* cuurently only uniform p elements are considered */
	pdeg = Pdeg_vec[0];

	/* find number of element shape functions and scalar dofs */
	num_shap = apr_get_el_pdeg_numshap(Field_id, El, &pdeg);
	ndofs = num_shap * nreq;

	/* initialize the matrices to zero */
	for (i = 0; i < ndofs * ndofs; i++) stiff_loc[i] = 0.0;
	for (i = 0; i < ndofs; i++) f_loc[i] = 0.0;

	/* get the coordinates of the nodes of El in the right order */
	mmr_el_node_coor(mesh_id, El, el_nodes, node_coor);

	/* prepare data for gaussian integration */
	apr_set_quadr_3D(base_q, &pdeg, &ngauss, xg, wg);

	/*kbw
	if(Mode<-1){
	printf("In l2_proj for element %d\n",El);
	printf("pdeg %d, ngauss %d\n",pdeg,ngauss);
	printf("nreq %d, ndof %d, local_dim %d\n",nreq,num_shap,ndofs);
	printf("%d nodes with coordinates:\n",el_nodes[0]);
	for(i=0;i<el_nodes[0];i++){
	  printf("node %d (global - %d): x - %f, y - %f, y - %f\n", i, el_nodes[i+1],
		node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
	}
	getchar();
	}
	/*kew*/

	for (ki = 0; ki < ngauss; ki++) {

		/* at the gauss point, compute basis functions, determinant etc*/
		iaux = 2; /* calculations with jacobian but not on the boundary */
		determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
						&xg[3 * ki], node_coor, NULL,
						base_phi, base_dphix, base_dphiy, base_dphiz,
						xcoor, NULL, NULL, NULL, NULL, NULL);

		vol = determ * wg[ki];

		/*kbw
		if(Mode<-1){
		printf("at gauss point %d, local coor %lf, %lf, %lf\n",
		ki,xg[3*ki],xg[3*ki+1],xg[3*ki+2]);
		printf("global coor %lf %lf %lf\n",xcoor[0],xcoor[1],xcoor[2]);
		printf("weight %lf, determ %lf, coeff %lf\n",
			wg[ki],determ,vol);
		getchar();
		}
		/*kew*/

		/* get value of projected function at integration point */
		if (Mode > 0) {
			//	(*Fun_p)(xcoor, fun_val, fun_dx, fun_dy, fun_dz);
		}
		else if (Mode <= -1) {

			/* we assume all elements are of the same type and have the same
			basis functions */

			/* initiate counter for Dofs */
			idofs = 0;

			/* for each element on a list */
			for (iel = 0; iel < -Mode; iel++) {
				iel_from = El_from[iel];

				/* find local coordinates within the element */

				ifound = mmr_loc_loc(mesh_id, El, &xg[3 * ki], iel_from, xloc);

				/*kbw
				if(Mode<-1){
				if(ifound==1) printf("in element %d, pdeg %d, idofs %d, at local coor %lf, %lf, %lf\n",
				iel_from,pdeg_from[iel], idofs,xloc[0],xloc[1],xloc[2]);
				else printf("in element %d, pdeg %d, idofs %d - not found\n",
				iel_from,pdeg_from[iel], idofs);
				printf("Dofs_from:\n");
				for (i=0;i<ndofs;i++) {
				  printf("%20.15lf",Dofs_from[idofs+i]);
				}
				printf("\n");
				}
				/*kew*/


				if (ifound == 1) {

					/* find value of function at integration point */
					iaux = 1; /* calculations without jacobian */

					apr_elem_calc_3D(iaux, nreq, &Pdeg_vec_from[iel],
							base_q, xloc, NULL, &Dofs_from[idofs],
							base_aux, NULL, NULL, NULL,
							NULL, value, NULL, NULL, NULL, NULL);

					break;
				}

				/* update counter for dofs */
				idofs += nreq * apr_get_el_pdeg_numshap(Field_id, iel, &Pdeg_vec_from[iel]);

			}

		}
		else {
			printf("Unknown mode for L2 proection!\n");
			exit(-1);
		}

		/*kbw
		if(Mode<-1){
		  printf("found value(s) for projection: ");
		  for(ieq=0;ieq<nreq;ieq++){
			printf("%20.15lf", value[ieq]);
		  }
		  printf("\n");
		}
		/*kew*/

		for (ieq = 0; ieq < nreq; ieq++) {

			kk = (ieq * ndofs + ieq) * num_shap;

			for (jdofs = 0; jdofs < num_shap; jdofs++) {
				for (idofs = 0; idofs < num_shap; idofs++) {


					/* mass matrix */
					stiff_loc[kk + idofs] +=
							base_phi[jdofs] * base_phi[idofs] * vol;

				}/* idofs */
				kk += ndofs;

			} /* jdofs */
		}/* ieq */

		kk = 0;
		for (ieq = 0; ieq < nreq; ieq++) {
			for (idofs = 0; idofs < num_shap; idofs++) {

				/* right hand side vector */
				f_loc[kk] += value[ieq] * base_phi[idofs] * vol;
				kk++;

			}/* idofs */
		}/* ieq */

	} /* ki */

	/*kbw
	if(El>0){
	printf("Stiffness matrix:\n");
	for (idofs=0;idofs<ndofs;idofs++) {
	  for (jdofs=0;jdofs<ndofs;jdofs++) {
		printf("%20.15lf",stiff_loc[idofs+jdofs*ndofs]);
	  }
	  printf("\n");
	}
	printf("F_loc:\n");
	for (idofs=0;idofs<ndofs;idofs++) {
	  printf("%20.15lf",f_loc[idofs]);
	}
	printf("\n");
	getchar();
	}
	kew*/

	/* solve the local problem */
	dgetrf_(&ndofs, &ndofs, stiff_loc, &ndofs, ipiv, &iaux);
	mf_check(iaux == 0, "DGETRF exited with error %d", iaux);
	dgetrs_("N", &ndofs, &ione, stiff_loc, &ndofs, ipiv, f_loc, &ndofs, &iaux);

	/*kbw
	if(El>0){
	printf("In l2_proj for element %d, ",El);
	printf("Solution:\n");
	for (idofs=0;idofs<ndofs;idofs++) {
	  printf("%20.15lf",f_loc[idofs]);
	}
	printf("\n");
	getchar();
	}
	kew*/

	/* rewrite solution */
	if (Dofs == NULL) {
		/* cannot be stored to global data structure based on vertices (FEM nodes) */
	}
	else {
		/* to vector Dofs */
		for (i = 0; i < ndofs; i++) Dofs[i] = f_loc[i];
	}

	/* check for Mode = -1 */
#ifdef DEBUG_APM
	if (Mode == -1 && Dofs != NULL) {

		/* in a loop over integration points */
		for (ki = 0; ki < ngauss; ki++) {

			/* at a gauss point compute solution etc*/
			iaux = 1; /* calculations without jacobian */
			apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
					&xg[3 * ki], NULL, Dofs,
					base_phi, NULL, NULL, NULL,
					NULL, u_val, NULL, NULL, NULL, NULL);

			/* get solution in element_from */
			iel_from = *El_from;

			/* find local coordinates within the element */
			ifound = mmr_loc_loc(mesh_id, El, &xg[3 * ki], iel_from, xloc);

			if (ifound == 1) {

				/* find value of function at integration point */
				iaux = 1; /* calculations without jacobian */
				apr_elem_calc_3D(iaux, nreq, &Pdeg_vec_from[iel_from],
						base_q, xloc, NULL, Dofs_from,
						base_phi, NULL, NULL, NULL,
						NULL, value, NULL, NULL, NULL, NULL);

			}
			else {
				printf("Something wrong in checking L2 proj!\n");
				exit(-1);
			}

			/*kbw
			if(El>0){
			printf("Checking l2_proj from element %d (pdeg %d0 to elemnt %d, (pdeg %d)",
			*El_from,*Pdeg_from,El,Pdeg);
			printf("Dofs_from:\n");
			for (idofs=0;idofs<ndofs;idofs++) {
			  printf("%20.15lf",Dofs_from[idofs]);
			}
			printf("\n");
			printf("Dofs_to:\n");
			for (idofs=0;idofs<ndofs;idofs++) {
			  printf("%20.15lf",Dofs[idofs]);
			}
			printf("\n");
			printf("Data:\n");
			for(ieq=0;ieq<nreq;ieq++){
			  printf("%20.15lf",value[ieq]);
			}
			printf("\n");
			printf("Solution:\n");
			for(ieq=0;ieq<nreq;ieq++){
			  printf("%20.15lf",u_val[ieq]);
			}
			printf("\n");
			getchar();
			}
			kew*/

			for (ieq = 0; ieq < nreq; ieq++) {

				if (fabs(value[ieq] - u_val[ieq]) > APC_SMALL) {
					printf("Something wrong in checking L2 proj! %lf != %lf\n",
							value[ieq], u_val[ieq]);
					exit(-1);
				}
			}

		}
	}
#endif

	return (1);
}


/*inline*/ int apr_is_pt_in_t4elem(// returns 0 - false, 1- true
		const int Mesh_id,
		const int El_id,
		const double Xglob[],
		double Xloc[],
		const double node_coor[],
		double D[],
		const double Precision)
{
	double amat[16] = {0.0};
	int ret_val = 0;
	int	fOrient[4] = {0};
	/* find coefficients of the transformation by solving a system of
	   linear equations */
	int ndofs = 4; /* four unknown coefficients of the transformation */
	int i = 0, col, row;

	mmr_el_faces(Mesh_id, El_id, NULL, fOrient);

	for (col = 0; col < 3; ++col) {
		for (row = 0; row < 4; ++row, ++i) {
			amat[i] = node_coor[3 * row + col];
		}
	}
	// amat stored COLUMNWISE
	for (row = 0; row < 4; row++) {
		amat[i++] = 1.0;
	}

	utr_mat_det(amat, ndofs, 'c', &D[0]);

	for (i = 0; i < 3; ++i) {
		amat[4 * i] = Xglob[i];
	}
	utr_mat_det(amat, ndofs, 'c', &D[1]);
	//D1*=fOrient[0];
	assert(D[0] != 0.0);
	if ((D[0] > 0 && D[1] >= -Precision) || (D[0] < 0 && D[1] <= Precision)) {
		for (i = 0; i < 3; ++i) {
			amat[4 * i] = node_coor[i];
			amat[1 + 4 * i] = Xglob[i];
		}

		utr_mat_det(amat, ndofs, 'c', &D[2]);
		//D2*=fOrient[1];
		if ((D[0] > 0 && D[2] >= -Precision) || (D[0] < 0 && D[2] <= Precision)) {
			for (i = 0; i < 3; ++i) {
				amat[1 + 4 * i] = node_coor[3 + i];
				amat[2 + 4 * i] = Xglob[i];
			}
			utr_mat_det(amat, ndofs, 'c', &D[3]);
			//D3*=fOrient[2];
			if ((D[0] > 0 && D[3] >= -Precision) || (D[0] < 0 && D[3] <= Precision)) {
				for (i = 0; i < 3; ++i) {
					amat[2 + 4 * i] = node_coor[6 + i];
					amat[3 + 4 * i] = Xglob[i];
				}
				utr_mat_det(amat, ndofs, 'c', &D[4]);
				//D4*=fOrient[3];
				if ((D[0] > 0 && D[4] >= -Precision) || (D[0] < 0 && D[4] <= Precision)) {
					assert(fabs(D[0] - (D[1] + D[2] + D[3] + D[4])) < Precision);
					ret_val = 1;

				}
			}
		}
	}

	return ret_val;
}

//int apr_sol_xglob_el(
//        const int Field_id, /* in: field ID */
//        const int El,       //in: number of element to check
//        const double Xglob[3],  /* in: global coordinates of a point */
//        double Xloc[3], /* out: list of local coordinates within element */
//        double *Sol,	/* out: list of solutions at the point */
//        double *Dxsol,  /* out: list of derivatives wrt x of solution */
//        double *Dysol,  /* out: list of derivatives wrt y of solution */
//        double *Dzsol   /* out: list of derivatives wrt z of solution */
//        )
//{
//    /* if solution required */
//    if (Sol!=NULL) {
//        /* find degree of polynomial and number of scalar dofs */
//        apr_get_el_pdeg(Field_id, El, &pdeg);
//        /* get the coordinates of nodes of El in the right order */
//        mmr_el_node_coor(mesh_id,El,el_nodes,node_coor);
//        /* get the most recent solution degrees of freedom */
//        apr_get_el_dofs(Field_id, El, Nb_sol, dofs_loc2);
//        base_q=apr_get_base_type(Field_id, El);
//        num_shap = apr_get_el_pdeg_numshap(Field_id, El, &pdeg);
//        ndofs = nreq*num_shap;
//        if (Dxsol!=NULL) {
//            /* find value of solution at a point */
//            iaux = 2; /* calculations with jacobian */
//            jaux = (nrel-1)*nreq;
//            apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
//                             &Xloc[3*nrel_check],node_coor,dofs_loc2,
//                    base_phi,base_dphix,base_dphiy,base_dphiz,
//                    xcoor,&Sol[jaux],&Dxsol[jaux],
//                    &Dysol[jaux],&Dzsol[jaux],NULL);
//        } else {
//            /* find value of solution at a point */
//            iaux = 1; /* calculations without jacobian */
//            jaux = (nrel-1)*nreq;
//            apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
//                             &Xloc[3*nrel_check],node_coor,dofs_loc2,
//                    base_phi,NULL,NULL,NULL,
//                    xcoor,&Sol[jaux],NULL,NULL,NULL,NULL);

//        }
//        mf_check( fabs(Xglob[0]-xcoor[0]) < SMALL, "Global solution(sol_xglob) make no sense..." );
//        mf_check( fabs(Xglob[1]-xcoor[1]) < SMALL, "Global solution(sol_xglob) make no sense..." );
//        mf_check( fabs(Xglob[2]-xcoor[2]) < SMALL, "Global solution(sol_xglob) make no sense..." );
//        mf_check( Sol[jaux] < 1.0e200 && Sol[jaux] > -1.0e200, "Global solution(sol_xglob) make no sense..." );
//    }
//    return 1;
//}

/*---------------------------------------------------------
apr_sol_xglob - to return the solution at a point with global
	coordinates specified. The procedure finds,
	for a given global point, the respective local coordinates
	within the proper initial mesh element, than computes
	corresponding local coordinates within an active ancestor of
	the initial element and finally finds the value of solution.
	There may be several initial mesh elements, several ancestors
	and several values due to the discontinuity of approximate solution.
---------------------------------------------------------*/
int apr_sol_elem(//INPUT:
		const int Field_id, const int Elem_id, const double Local_point[3], const int Solution_number,
		// OUTPUT:
		double * Solution, double * Dxsol, double * Dysol, double * Dzsol,
		double Global_point_coords[3])
{
//    if (Solution!=NULL) {

//        /* find degree of polynomial and number of scalar dofs */
//        int pdeg=0;
//        apr_get_el_pdeg(Field_id, Elem_id, &pdeg);

//        /* get the coordinates of nodes of El in the right order */
//        const int mesh_id = apr_get_mesh_id(Field_id);
//        int el_nodes[MMC_MAXELVNO]={0}, node_coor[MMC_MAXELVNO*3]={0};
//        mmr_el_node_coor(mesh_id,Elem_id,el_nodes,node_coor);

//        /* get the most recent solution degrees of freedom */
//        double dofs_loc2[APC_MAXELSD]={0.0};
//        apr_get_el_dofs(Field_id, Elem_id, Solution_number, dofs_loc2);
//        base_q=apr_get_base_type(Field_id, Elem_id);
//        num_shap = apr_get_el_pdeg_numshap(Field_id, Elem_id, &pdeg);
//        ndofs = nreq*num_shap;
//        if (Dxsol!=NULL) {
//            /* find value of solution at a point */
//            jaux = (nrel-1)*nreq;
//            apr_elem_calc_3D(APC_ELEM_CALC_DERIVATIVES_AND_JACOBIAN, nreq, &pdeg, base_q,
//                             &Local_point[3*nrel_check],node_coor,dofs_loc2,
//                    base_phi,base_dphix,base_dphiy,base_dphiz,
//                    xcoor,&Solution[jaux],&Dxsol[jaux],
//                    &Dysol[jaux],&Dzsol[jaux],NULL);
//        } else {
//            /* find value of solution at a point */
//            jaux = (nrel-1)*nreq;
//            apr_elem_calc_3D(APC_ELEM_CALC_SHAPE_FUNC_N_VALUES, nreq, &pdeg, base_q,
//                             &Local_point[3*nrel_check],node_coor,dofs_loc2,
//                    base_phi,NULL,NULL,NULL,
//                    xcoor,&Solution[jaux],NULL,NULL,NULL,NULL);

//        }
//        mf_check( fabs(Global_point_coords[0]-xcoor[0]) < SMALL, "Global solution(sol_xglob) make no sense..." );
//        mf_check( fabs(Global_point_coords[1]-xcoor[1]) < SMALL, "Global solution(sol_xglob) make no sense..." );
//        mf_check( fabs(Global_point_coords[2]-xcoor[2]) < SMALL, "Global solution(sol_xglob) make no sense..." );
//        mf_check( Solution[jaux] < 1.0e200 && Solution[jaux] > -1.0e200, "Global solution(sol_xglob) make no sense..." );
//    }

	return 0;
}



int 	apr_sol_xglob( /* returns: >0 - success, <=0 - failure */
		int Field_id,   /* in: field ID */
		double * Xglob,	/* in: global coordinates of a point */
		int Nb_sol,    /* in: which solution to take: 1 - sol_1, 2 - sol_2 */
		int * El,	/* out: list of element numbers,  */
		/*      El[0] - number of elements on the list */
		double * Xloc,	/* out: list of local coordinates within elements */
		double * Sol,	/* out: list of solutions at the point */
		double * Dxsol, /* out: list of derivatives wrt x of solution */
		double * Dysol, /* out: list of derivatives wrt y of solution */
		double * Dzsol,  /* out: list of derivatives wrt z of solution */
		double Close_proximity_precision, //< default to APC_CLOSE; how far from element point is considered to be within an element
		int Sol_xglob_flags //< see ape_sol_xlgob
)
{
	/* pointer to field structure */
	apt_field * field_p;
	int mesh_id;

	int i, iel, iaux, jaux, row, col, isys, pdeg, nrel, el_aux; //j,
	int nmel, ndofs, base_q, ifound, nreq, num_shap;

	int ipiv[4] = {0}, el_son = 0, ison = 0, elsons[MMC_MAXELSONS + 1] = {0};
	int nrel_check = 0, el_check[10 * MMC_MAXELSONS + 1] = {0};
	double bvec[4] = {0.0}, atrans[9] = {0.0}, btrans[3] = {0.0};

	double dofs_loc[APC_MAXELSD] = {0.0};
	double dofs_loc2[APC_MAXELSD] = {0.0};
	double base_phi[APC_MAXELVD] = {0.0};
	double base_dphix[APC_MAXELVD] = {0.0}; /* x-derivatives of basis function */
	double base_dphiy[APC_MAXELVD] = {0.0}; /* y-derivatives of basis function */
	double base_dphiz[APC_MAXELVD] = {0.0}; /* y-derivatives of basis function */
	double node_coor[3 * MMC_MAXELVNO] = {0.0};
	double xloc[10 * MMC_MAXELSONS * 3] = {0.0}, xloc_aux[3] = {0.0}, xcoor[3] = {0.0};
	int el_nodes[MMC_MAXELVNO + 1];      /* list of nodes of El */
	const double ref_points_prism[6][3] =
			/* data for reference points for three dimensional prismatic element */
	{	{ 0.0, 0.0, -1.0 }, { 1.0, 0.0, -1.0 }, { 0.0, 1.0, -1.0 },
		{ 0.0, 0.0,  1.0 }, { 1.0, 0.0,  1.0 }, { 0.0, 1.0,  1.0 }
	};

	int ione = 1, keep_searching = 1;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (Close_proximity_precision < 1e-15) Close_proximity_precision = 1e-15;

	/* select the proper approximation field */
	field_p = apr_select_field(Field_id);

	/* get the corresponding mesh ID */
	mesh_id = field_p->mesh_id;

	/* get necessary mesh parameters */
	nmel = mmr_get_max_elem_id(mesh_id);

	/* get formulation parameters */
	nreq = apr_get_nreq(Field_id);

	/* prepare the matrix for finding coefficients of linear transformation */

	nrel = 0;
	int first_id = 1;
	if (Sol_xglob_flags & APE_SOL_XGLOB_CHECK_ONLY_GIVEN_ELEMENT) {
		mf_check(El[0] > 0, "Given elem to check have invalid id=%d!", El[0]);
		first_id = El[0];
		nmel = El[0];

	}

	/* loop over all initial mesh elements */
	for (iel = first_id; iel <= nmel && keep_searching; ++iel) {
		if (mmr_el_status(mesh_id, iel) != MMC_FREE && mmr_el_gen(mesh_id, iel) == 0) {
			int nodes[10] = {0};

			/* get the coordinates of element nodes in the right order */
			mmr_el_node_coor(mesh_id, iel, nodes, node_coor);

			/*kbw
			if(iel>0){
						printf("In sol_xglob for element %d\n",iel);
			printf("nodes with coordinates:\n");
						for(i=0;i<4;i++){
			  printf("node %d : x - %f, y - %f, y - %f\n", i,
				node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
			}
			getchar();
			}
			/*kew*/


			/*kb!!!
			we assume that the geometrical transformation from the master element
			is linear
			kb!!!*/

			if (mmr_el_type(mesh_id, iel) == MMC_TETRA) {
				double D[5] = {0.0};
				if (apr_is_pt_in_t4elem(mesh_id, iel, Xglob, Xloc, node_coor, D, Close_proximity_precision) == 1) {
					if ( (mmr_el_type_ref(mesh_id, iel) == MMC_NOT_REF)
							|| !(Sol_xglob_flags & APE_SOL_XGLOB_MATCH_WITH_ADAPTATION) ) {
						++nrel;
						El[nrel] = iel;
						el_aux = iel;
					}
					else {
						int * el_sons = (int *) calloc(MMC_MAXELSONS * mmr_get_max_gen(mesh_id), sizeof(int));
						mmr_el_fam_all(mesh_id, iel, el_sons);
						el_aux = 0; // this means that target element not yet found
						for (i = 1; (i < el_sons[0]) && (el_aux == 0); ++i) {
							ison = el_sons[i];
							if (mmr_el_type_ref(mesh_id, ison) == MMC_NOT_REF) {
								mmr_el_node_coor(mesh_id, iel, nodes, node_coor);
								if (apr_is_pt_in_t4elem(mesh_id, ison, Xglob, Xloc, node_coor, D, Close_proximity_precision) == 1) {
									// we found real active element including given Xglob point
									++nrel; // we found active element
									El[nrel] = ison;
									el_aux = ison;
								}
							}
						}

						free(el_sons);
					}

					Xloc[3 * nrel_check  ] = D[2] / D[0];
					Xloc[3 * nrel_check + 1] = D[3] / D[0];
					Xloc[3 * nrel_check + 2] = D[4] / D[0];
//              assert(Xloc[3*nrel_check  ] > -Close_proximity_precision);
//              assert(Xloc[3*nrel_check+1] > -Close_proximity_precision);
//              assert(Xloc[3*nrel_check+2] > -Close_proximity_precision);
//              assert(Xloc[3*nrel_check  ]-1 < Close_proximity_precision);
//              assert(Xloc[3*nrel_check+1]-1 < Close_proximity_precision);
//              assert(Xloc[3*nrel_check+2]-1 < Close_proximity_precision);

					if (      Xloc[3 * nrel_check  ] < 0.0) Xloc[3 * nrel_check  ] = 0.0;
					else if ( Xloc[3 * nrel_check  ] > 1.0) Xloc[3 * nrel_check  ] = 1.0;
					if (      Xloc[3 * nrel_check + 1] < 0.0) Xloc[3 * nrel_check + 1] = 0.0;
					else if ( Xloc[3 * nrel_check + 1] > 1.0) Xloc[3 * nrel_check + 1] = 1.0;
					if (      Xloc[3 * nrel_check + 2] < 0.0) Xloc[3 * nrel_check + 2] = 0.0;
					else if ( Xloc[3 * nrel_check + 2] > 1.0) Xloc[3 * nrel_check + 2] = 1.0;

					/* if solution required */
					if (Sol != NULL) {
						/* find degree of polynomial and number of scalar dofs */
						apr_get_el_pdeg(Field_id, el_aux, &pdeg);
						/* get the coordinates of nodes of El in the right order */
						mmr_el_node_coor(mesh_id, el_aux, el_nodes, node_coor);
						/* get the most recent solution degrees of freedom */
						apr_get_el_dofs(Field_id, el_aux, Nb_sol, dofs_loc2);
						base_q = apr_get_base_type(Field_id, el_aux);
						num_shap = apr_get_el_pdeg_numshap(Field_id, el_aux, &pdeg);
						ndofs = nreq * num_shap;
						if (Dxsol != NULL) {
							/* find value of solution at a point */
							iaux = 2; /* calculations with jacobian */
							jaux = (nrel - 1) * nreq;
							apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
									&Xloc[3 * nrel_check], node_coor, dofs_loc2,
									base_phi, base_dphix, base_dphiy, base_dphiz,
									xcoor, &Sol[jaux], &Dxsol[jaux],
									&Dysol[jaux], &Dzsol[jaux], NULL);
						} else {
							/* find value of solution at a point */
							iaux = 1; /* calculations without jacobian */
							jaux = (nrel - 1) * nreq;
							apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
									&Xloc[3 * nrel_check], node_coor, dofs_loc2,
									base_phi, NULL, NULL, NULL,
									xcoor, &Sol[jaux], NULL, NULL, NULL, NULL);

						}
						mf_check_debug( fabs(Xglob[0] - xcoor[0]) < Close_proximity_precision, "Global solution(sol_xglob) make no sense..." );
						mf_check_debug( fabs(Xglob[1] - xcoor[1]) < Close_proximity_precision, "Global solution(sol_xglob) make no sense..." );
						mf_check_debug( fabs(Xglob[2] - xcoor[2]) < Close_proximity_precision, "Global solution(sol_xglob) make no sense..." );
						mf_check_debug( Sol[jaux] < 1.0e200 && Sol[jaux] > -1.0e200, "Global solution(sol_xglob) make no sense..." );
					}
					/* end rewriting solution */
					++nrel_check;
				}

			}
			else { // Assume mmr_el_type(mesh_id,iel)==MMC_PRISM.

				/* find coefficients of the transformation by solving a system of
				   linear equations */

				ndofs = 4; /* four unknown coefficients of the transformation */
				double amat[16] = {0.0};
				i = 0;
				for (col = 0; col < 3; col++) {
					for (row = 0; row < 4; row++) {
						amat[i++] = node_coor[3 * row + col];
					}
				}

				for (row = 0; row < 4; row++) {
					amat[i++] = 1.0;
				}

				/*kbw
				printf("before decomposition amat: \n");
				for(i=0;i<(ndofs)*(ndofs);i++) printf("%20.15lf",amat[i]) ;
				printf("\n");
				/*kew*/

				dgetrf_(&ndofs, &ndofs, amat, &ndofs, ipiv, &iaux);
				mf_check(iaux == 0, "DGETRF exited with error %d", iaux);

				/*kbw
				printf("after decomposition amat: \n");
				for(i=0;i<(ndofs)*(ndofs);i++) printf("%20.15lf",amat[i]) ;
				printf("\n");
				/*kew*/

				for (isys = 0; isys < 3; isys++) {

					for (row = 0; row < 4; row++) {
						bvec[row] = ref_points_prism[row][isys];
					}

					/*kbw
					printf("bvec: \n");
					for(i=0;i<4;i++) printf("%20.15lf",bvec[i]) ;
					printf("\n");
					/*kew*/

					/* solve the local problem */
					dgetrs_("N", &ndofs, &ione, amat, &ndofs, ipiv, bvec, &ndofs, &iaux);
					for (i = 0; i < 3; i++) atrans[3 * isys + i] = bvec[i] ;
					btrans[isys] = bvec[3];

					/*kbw
					printf("solution: \n");
					for(i=0;i<ndofs;i++) printf("%20.15lf",bvec[i]) ;
					printf("\n");
					/*kew*/


					/*kbw
							  for(row=0;row<3;row++){

								daux=0;
								for(col=0;col<3;col++){

								  daux += atrans[3*isys+col]*node_coor[3*row+col];

								}

								daux += btrans[isys];

								printf("Checking: %20.15lf ?= %20.15lf\n",
							daux,ref_points_prism[row][isys]);

							  }
					/*kew*/

				} /* end of loop over systems */

				/* checking remaining element vertices to confirm the transformation
				is linear */
				/*kbw
				printf("Vertex 5: \n");
						for(row=0;row<3;row++){

						  daux=0;
						  for(col=0;col<3;col++){

							daux += atrans[3*row+col]*node_coor[12+col];

						  }

						  daux += btrans[row];

						  printf("Checking coordinate %d : %20.15lf ?= %20.15lf\n",
						row, daux,ref_points_prism[4][row]);

						}
				printf("Vertex 6: \n");
						for(row=0;row<3;row++){

						  daux=0;
						  for(col=0;col<3;col++){

							daux += atrans[3*row+col]*node_coor[15+col];

						  }

						  daux += btrans[row];

						  printf("Checking coordinate %d : %20.15lf ?= %20.15lf\n",
						row, daux,ref_points_prism[5][row]);

						}
				/*kew*/


				ut_mat3vec(atrans, Xglob, xloc);
				for (i = 0; i < 3; i++) xloc[i] += btrans[i];

				/*checking*/
#ifdef DEBUG_APM
				iaux = 1; /* calculations without jacobian */
				pdeg = 101;
				apr_elem_calc_3D(iaux, nreq, &iaux, iaux,
						xloc, node_coor, NULL,
						NULL, NULL, NULL, NULL,
						xcoor, NULL, NULL, NULL, NULL, NULL);
				for (i = 0; i < 3; i++)
					if (fabs(Xglob[i] - xcoor[i]) > APC_SMALL) {
						printf("2Error in finding local coordinates for global point\n");
						return (-1);
					}

				/*kbw
				printf("Checking :\n");
				for(i=0;i<3;i++)
					printf("coord %d, local %15.12lf, global %15.12lf ?= %15.12lf\n",
						i,xloc[i],Xglob[i],xcoor[i]);
				/*kew*/
#endif

				/* if global point within initial mesh element */
				if ( xloc[0] > -Close_proximity_precision && xloc[1] > -Close_proximity_precision &&
						xloc[0] + xloc[1] < 1.0 + Close_proximity_precision &&
						xloc[2] > -1.0 - Close_proximity_precision && xloc[2] < 1.0 + Close_proximity_precision ) {

					if (xloc[0] < 0.0) xloc[0] = 0.0;
					if (xloc[1] < 0.0) xloc[1] = 0.0;
					if (xloc[0] + xloc[1] > 1.0) xloc[0] = 1.0 - xloc[1];
					if (xloc[2] < -1.0) xloc[2] = -1.0;
					if (xloc[2] > 1.0) xloc[2] = 1.0;

					/*kbw
					printf("Found respective point in initial element %d\n local coordinates:",iel);
					for(i=0;i<3;i++) printf("%15.12lf",xloc[i]);printf("\n");
					/*kew*/
					nrel_check = 1;
					el_check[0] = iel;

					while (nrel_check > 0) {

						nrel_check--;
						el_aux = el_check[nrel_check];

						/*kbw
						printf("checking element %d, waiting list:",el_aux);
						for(i=0;i<nrel_check;i++) printf(" %d",el_check[i]);
						printf("\n");
						/*kew*/

						while ( mmr_el_status(mesh_id, el_aux) <= 0) {

							for (i = 0; i < 3; i++) xloc_aux[i] = xloc[3 * nrel_check + i];

							mmr_el_fam(mesh_id, el_aux, elsons, NULL);

							/*find proper son and its local coordinates*/
							ifound = 0;
							if (xloc_aux[2] < Close_proximity_precision) {
								if (xloc_aux[0] > 0.5 - Close_proximity_precision) {
									el_check[nrel_check] = elsons[2];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
								if (xloc_aux[1] > 0.5 - Close_proximity_precision) {
									el_check[nrel_check] = elsons[3];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
								if (xloc_aux[0] + xloc_aux[1] < 0.5 + Close_proximity_precision)  {
									el_check[nrel_check] = elsons[1];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
								if (xloc_aux[0] + xloc_aux[1] > 0.5 - Close_proximity_precision &&
										xloc_aux[0] < 0.5 + Close_proximity_precision &&
										xloc_aux[1] < 0.5 + Close_proximity_precision )  {
									el_check[nrel_check] = elsons[4];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
							}
							if (xloc_aux[2] > -Close_proximity_precision) {
								if (xloc_aux[0] > 0.5 - Close_proximity_precision) {
									el_check[nrel_check] = elsons[6];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
								if (xloc_aux[1] > 0.5 - Close_proximity_precision) {
									el_check[nrel_check] = elsons[7];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
								if (xloc_aux[0] + xloc_aux[1] < 0.5 + Close_proximity_precision)  {
									el_check[nrel_check] = elsons[5];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
								if (xloc_aux[0] + xloc_aux[1] > 0.5 - Close_proximity_precision &&
										xloc_aux[0] < 0.5 + Close_proximity_precision &&
										xloc_aux[1] < 0.5 + Close_proximity_precision )  {
									el_check[nrel_check] = elsons[8];
									ifound += mmr_loc_loc(mesh_id, el_aux, xloc_aux,
													el_check[nrel_check], &xloc[3 * nrel_check]);
									nrel_check++;
								}
							}

							if (ifound <= 0) {
								printf("Local coordinates not found within a family\n");
								return (-1);
							}

							/*kbw
							printf("after checking family for element %d, waiting list:\n",el_aux);
							for(i=0;i<nrel_check;i++) {
							  double ncoor[20];int ii;
							  printf("element %d, xcoor %lf %lf %lf\n",
							  el_check[i],xloc[3*i],xloc[3*i+1],xloc[3*i+2]);
							  mmr_el_node_coor(mesh_id, el_check[i], NULL, ncoor);
							  printf("nodes with coordinates:\n");
							  for(ii=0;ii<6;ii++){
								printf("node %d : x - %f, y - %f, y - %f\n", ii,
									 ncoor[3*ii],ncoor[3*ii+1],ncoor[3*ii+2]);
							  }
							}
							getchar();
							/*kew*/

							--nrel_check;
							el_aux = el_check[nrel_check];

						} /* end while elements inactive */

						/* found active element with global coordinates Xglob */
						++nrel;
						El[nrel] = el_aux;
						for (i = 0; i < 3; i++) Xloc[(nrel - 1) * 3 + i] = xloc[3 * nrel_check + i];

						/*kbw
						printf("found next (%d) element %d, coor %lf %lf %lf\n",
							   nrel,El[nrel],Xloc[(nrel-1)*3],Xloc[(nrel-1)*3+1],Xloc[(nrel-1)*3+2]);

						getchar();
						/*kew*/

						/* if solution required */
						if (Sol != NULL) {

							/* find degree of polynomial and number of element scalar dofs */
							apr_get_el_pdeg(Field_id, el_aux, &pdeg);

							/* get the coordinates of the nodes of El in the right order */
							mmr_el_node_coor(mesh_id, el_aux, el_nodes, node_coor);

							/* get the most recent solution degrees of freedom */
							apr_get_el_dofs(Field_id, el_aux, Nb_sol, dofs_loc2);

							/*kbw
							printf("found next (%d) element %d, coor %lf %lf %lf\n",
								   nrel,El[nrel],Xloc[(nrel-1)*3],Xloc[(nrel-1)*3+1],Xloc[(nrel-1)*3+2]);
							printf("dofs:\n");
							for(j=0;j<6;j++){
								for(i=0;i<3;i++){
									printf("%lf ", dofs_loc2[(j)*3+i]);
								}
								printf("\n");
							}
							getchar();
							/*kew*/

							base_q = apr_get_base_type(Field_id, el_aux);
							num_shap = apr_get_el_pdeg_numshap(Field_id, el_aux, &pdeg);
							ndofs = nreq * num_shap;

							if (Dxsol != NULL) {

								/* find value of solution at a point */
								iaux = 2; /* calculations with jacobian */
								jaux = (nrel - 1) * nreq;
								apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
										&xloc[3 * nrel_check], node_coor, dofs_loc2,
										base_phi, base_dphix, base_dphiy, base_dphiz,
										xcoor, &Sol[jaux],
										&Dxsol[jaux], &Dysol[jaux], &Dzsol[jaux], NULL);

							}
							else {

								/* find value of solution at a point */
								iaux = 1; /* calculations without jacobian */
								jaux = (nrel - 1) * nreq;
								apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
										&xloc[3 * nrel_check], node_coor, dofs_loc2,
										base_phi, NULL, NULL, NULL,
										xcoor, &Sol[jaux], NULL, NULL, NULL, NULL);

							}

							/*kbw
							{
							double f, fx,fy,fz, lalp;

							printf("\nChecking element %d:\n", el_aux);
							for(i=0;i<3;i++)
								printf("coord %d, local %15.12lf, global %15.12lf ?= %15.12lf\n",
									i,xloc[3*nrel_check+i],Xglob[i],xcoor[i]);

							//dgpr_exact_sol(name,xcoor[0],xcoor[1],xcoor[2],&f,&fx,&fy,&fz,&lalp);
							for(i=0;i<nreq;i++){
							printf("Solution (%d) : exact %20.15lf, approximate %20.15lf\n",
								jaux+i,f,Sol[jaux+i]);
							if(Dxsol!=NULL){
							printf("x derivative  : exact %20.15lf, approximate %20.15lf\n",
								fx,Dxsol[jaux+i]);
							printf("y derivative  : exact %20.15lf, approximate %20.15lf\n",
								fy,Dysol[jaux+i]);
							printf("z derivative  : exact %20.15lf, approximate %20.15lf\n",
								fz,Dzsol[jaux+i]);
							}
							}
							getchar();

							}

							/*kew*/

							/*checking*/
#ifdef DEBUG_APM
							for (i = 0; i < 3; i++) {
								if (fabs(Xglob[i] - xcoor[i]) > APC_SMALL) {
									printf("3Error in finding local coordinates for global point\n");
									return (-1);
								}
							}
#endif

						} /* end rewriting solution */

					} /* end different branches from initial mesh element */

				} /* end if point found within initial mesh lement */

			} // end else if type==MMC_PRISM

		} /* end if element in initial mesh */

		// Check if more elements are required, or just one.
		if (! (Sol_xglob_flags & APE_SOL_XGLOB_MATCH_ALL_ELEMENTS) ) {
			if (nrel > 0) {
				keep_searching = 0; // break searching if found at least one element
			}
		}

	} /* end of loop over elements */

	El[0] = nrel;

	/*kbw
	printf("Solution:");
	for(i=0;i<nreq*nrel;i++){
	printf("%20.15lf",Sol[i]);
	}
	printf("\n");
	/*kew*/

	return (1);
}


int apr_get_profile( // returns: >0 - success, <=0 - failure
		FILE * filePtr, //in: pointer to file, where profile will be printed
		int fieldId,		//in: field ID
		int solNr,			//in: solution number
		int nSol,			//in: solution length (no. of components)
		double * pt1,		//in: begin point coords
		double * pt2,		//in: end point coords
		int nPoints			//in: number of points (including end points)
)
{
	if (filePtr != NULL
			&& fieldId > 0
			&& solNr > 0
			&& nSol > 0
			&& pt1 != NULL
			&& pt2 != NULL
			&& nPoints > 1)
	{
		double nDp = nPoints - 1;
		int p = 0, s = 0;
		int el[30] = {0};
		double dp[3] = {(pt2[0] - pt1[0]) / nDp, (pt2[1] - pt1[1]) / nDp, (pt2[2] - pt1[2]) / nDp};
		double coord[3] = {pt1[0], pt1[1], pt1[2]};
		double solution[100] = {0.0}, xloc[100] = {0.0};

		fprintf(filePtr, "x\ty\tz\t");
		for (s = 0; s < nSol; ++s)
		{
			fprintf(filePtr, "sol%d\t", s);
		}

		for (p = 0; p < nPoints; ++p)
		{
			/*kbw
			  printf("\nfield_id %d, coord %lf, %lf, %lf, solNr %d\n",
					   fieldId,coord[0],coord[1],coord[2],solNr  );
			/*kew*/

			int iaux = apr_sol_xglob(fieldId, coord, solNr, el,
							xloc, solution, NULL, NULL, NULL, APC_CLOSE,
							APE_SOL_XGLOB_DEFAULT
							| APE_SOL_XGLOB_MATCH_WITH_ADAPTATION);
			/*kbw
			  printf("\nfield_id %d, coord %lf, %lf, %lf, solNr %d solution[0] %lf\n",
					 fieldId,coord[0],coord[1],coord[2],solNr ,solution[0] );

			  for(s=0;s<nSol; ++s)
			  {
			  printf("ieq %d, %Lf\t",s, solution[s]);	// val print
			  }
			/*kew*/

			if (iaux == 1) { // if point found
				fprintf(filePtr, "\n%lf\t%lf\t%lf\t", coord[0], coord[1], coord[2]);
				// coords print
				for (s = 0; s < nSol; ++s)
				{
					fprintf(filePtr, "%lf\t", solution[s]);	// val print
					solution[s] = 0.0; // reset to avoid printing incorrenct values
				}
			}
			coord[0] += dp[0];
			coord[1] += dp[1];
			coord[2] += dp[2];
		}
	}
	return 1;
}


/* Internal utilities */

/*---------------------------------------------------------
ut_mat3_inv - to invert a 3x3 matrix (stored as a vector!)
---------------------------------------------------------*/
double ut_mat3_inv(	/* returns: determinant of matrix to invert */
		double * mat,	/* matrix to invert */
		double * mat_inv	/* inverted matrix */
)
{
	double s0, s1, s2, rjac, rjac_inv;

	/* rjac = mat[0]*mat[4]*mat[8] + mat[3]*mat[7]*mat[2] */
	/*      + mat[6]*mat[1]*mat[5] - mat[6]*mat[4]*mat[2] */
	/*      - mat[0]*mat[7]*mat[5] - mat[3]*mat[1]*mat[8]; */

	rjac = mat[0] * (mat[4] * mat[8] - mat[7] * mat[5]) +
			mat[3] * (mat[7] * mat[2] - mat[1] * mat[8]) +
			mat[6] * (mat[1] * mat[5] - mat[4] * mat[2]);

#ifdef DEBUG_APM
	if (fabs(rjac) < 1.0e-15) {
		printf("Small jacobian of matrix being inverted %.15lf\n", rjac);
		exit(-1);
	}
#endif

	rjac_inv = 1.0 / rjac;

	s0 =   mat[4] * mat[8] - mat[7] * mat[5];
	s1 = - mat[3] * mat[8] + mat[6] * mat[5];
	s2 =   mat[3] * mat[7] - mat[6] * mat[4];

	mat_inv[0] = s0 * rjac_inv;
	mat_inv[3] = s1 * rjac_inv;
	mat_inv[6] = s2 * rjac_inv;

	s0 =   mat[7] * mat[2] - mat[1] * mat[8];
	s1 =   mat[0] * mat[8] - mat[6] * mat[2];
	s2 = - mat[0] * mat[7] + mat[6] * mat[1];

	mat_inv[1] = s0 * rjac_inv;
	mat_inv[4] = s1 * rjac_inv;
	mat_inv[7] = s2 * rjac_inv;

	s0 =   mat[1] * mat[5] - mat[4] * mat[2];
	s1 = - mat[0] * mat[5] + mat[3] * mat[2];
	s2 =   mat[0] * mat[4] - mat[3] * mat[1];

	mat_inv[2] = s0 * rjac_inv;
	mat_inv[5] = s1 * rjac_inv;
	mat_inv[8] = s2 * rjac_inv;

	return (rjac);
}

/*---------------------------------------------------------
ut_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
void ut_vec3_prod(
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
ut_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
double ut_vec3_mxpr( /* returns: mixed product [a,b,c] */
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
ut_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
double ut_vec3_length(	/* returns: vector length */
		double * vec	/* in: vector */
)
{
	return (sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
}

/*---------------------------------------------------------
ut_mat3vec - to compute matrix vector product in 3D space
---------------------------------------------------------*/
void ut_mat3vec(
		double * m1, 	/* in: matrix (stored by rows as a vector!) */
		double * v1, 	/* in: vector */
		double * v2	/* out: resulting vector */
)
{

	v2[0] = m1[0] * v1[0] + m1[1] * v1[1] + m1[2] * v1[2] ;
	v2[1] = m1[3] * v1[0] + m1[4] * v1[1] + m1[5] * v1[2] ;
	v2[2] = m1[6] * v1[0] + m1[7] * v1[1] + m1[8] * v1[2] ;
	return;
}

/*---------------------------------------------------------
ut_mat3mat - to compute matrix matrix product in 3D space
	(all matrices are stored by rows as vectors!)
---------------------------------------------------------*/
void ut_mat3mat(
		double * m1,	/* in: matrix */
		double * m2,	/* in: matrix */
		double * m3	/* out: matrix m1*m2 */
)
{

	m3[0] = m1[0] * m2[0] + m1[1] * m2[3] + m1[2] * m2[6] ;
	m3[1] = m1[0] * m2[1] + m1[1] * m2[4] + m1[2] * m2[7] ;
	m3[2] = m1[0] * m2[2] + m1[1] * m2[5] + m1[2] * m2[8] ;

	m3[3] = m1[3] * m2[0] + m1[4] * m2[3] + m1[5] * m2[6] ;
	m3[4] = m1[3] * m2[1] + m1[4] * m2[4] + m1[5] * m2[7] ;
	m3[5] = m1[3] * m2[2] + m1[4] * m2[5] + m1[5] * m2[8] ;

	m3[6] = m1[6] * m2[0] + m1[7] * m2[3] + m1[8] * m2[6] ;
	m3[7] = m1[6] * m2[1] + m1[7] * m2[4] + m1[8] * m2[7] ;
	m3[8] = m1[6] * m2[2] + m1[7] * m2[5] + m1[8] * m2[8] ;

}

/*---------------------------------------------------------
ut_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int ut_chk_list(	/* returns: */
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


/** @}*/
