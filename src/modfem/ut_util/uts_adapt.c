/************************************************************************
File uts_adapt - utility routines for adaptivity (common to all problem
		 modules, possibly used also by other modules)

  utr_create_patches - to create patches of elements containing a given node
			   (uses only mesh data - can be the same for all fields)
  utr_create_patches_small - the same as above but uses lists of elements nodes
							 instead of lists of elements dof entities
  utr_recover_derivatives - to recover deirvatives using patches of elements
			  (can be called separately for each field on the same mesh)
  utr_recover_derivatives_small-to recover derivatives using patches of elements
	(does not use apr_get_stiff_mat_data for getting lists of element nodes)
  utr_adapt - to enforce default adaptation strategy for SINGLE problem/field
  utr_test_refinements - to perform sequence of random refinements/unrefinements
  utr_manual_refinement - to perform manual refinement/unrefinement

------------------------------
History:
	08.2008 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>

/* interface for all approximation modules */
#include <modfem/aph_intf.h>

#include <modfem/uth_system.h>

#ifdef PARALLEL
/* interface for parallel mesh manipulation modules */
#include "mmph_intf.h"

/* interface with parallel communication library */
#include "pch_intf.h"
#endif

/* interface for general purpose utilities - for all problem dependent modules*/
#include <modfem/uth_intf.h>

/* interface for linear algebra packages */
#include <modfem/lin_alg_intf.h>

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>

#include <modfem/pdh_intf.h>

#include <modfem/uth_log.h>

// maximal number of solution components
#define UTC_MAXEQ PDC_MAXEQ


void utr_print_adapt_info(const int Field_id,
		FILE * Interactive_output) // in)
{
	const int mesh_id = apr_get_mesh_id(Field_id);
	const int nrno = mmr_get_nr_node(mesh_id);
	const int nmno = mmr_get_max_node_id(mesh_id);
	const int nred = mmr_get_nr_edge(mesh_id);
	const int nmed = mmr_get_max_edge_id(mesh_id);
	const int nrfa = mmr_get_nr_face(mesh_id);
	const int nmfa = mmr_get_max_face_id(mesh_id);
	const int nrel = mmr_get_nr_elem(mesh_id);
	const int nmel = mmr_get_max_elem_id(mesh_id);

	fprintf(Interactive_output, "\nBefore/after adaptation. Mesh info:\n");
	fprintf(Interactive_output, "Parameters (number of active, maximal index):\n");
	fprintf(Interactive_output, "Elements: nrel %d, nmel %d\n", nrel, nmel);
	fprintf(Interactive_output, "Faces:    nrfa %d, nmfa %d\n", nrfa, nmfa);
	fprintf(Interactive_output, "Edges:    nred %d, nmed %d\n", nred, nmed);
	fprintf(Interactive_output, "Nodes:    nrno %d, nmno %d\n", nrno, nmno);

	fprintf(Interactive_output, "Max_gen %d, Maxgendiff %d\n",
			mmr_get_max_gen(mesh_id), mmr_get_max_gen_diff(mesh_id));

	char approx_name[255];
	apr_module_introduce(approx_name);
	int nr_dof_ents = 0;
	if (strncmp(approx_name, "STANDARD_LINEAR", 15) == 0) {
		nr_dof_ents = nrno;
	}
	else if (strncmp(approx_name, "STANDARD_QUADRATIC", 18) == 0) {
		int pdeg = apr_get_el_pdeg(Field_id, 0, NULL);
		if (pdeg == 101) {
			nr_dof_ents = nrno;
		}
		if (pdeg == 202) {
			nr_dof_ents = nrno + nred;
		}
	}
	else if (strncmp(approx_name, "DG_SCALAR_PRISM ", 15) == 0) {
		// for DG approximation
		nr_dof_ents = nrel;
	}
	else {
		printf("Unknown approximation module in utr_print_adapt_info! Exiting.\n");
		exit(-1);
	}

	fprintf(Interactive_output, "Number of DOF entities:  %d\n", nr_dof_ents );
}

/*---------------------------------------------------------
  utr_create_patches - to create patches of elements containing a given node
---------------------------------------------------------*/
int utr_create_patches( /* returns: >0 - Nr_patches, <=0 - failure */
		int Field_id,       /* in: field data structure to be used */
		utt_patches ** patches_p /* in/out - array of patches for real nodes */
)
{
	int mesh_id, nreq, nr_patches;
	int i, ino, nno, idofent, dof_ent_id, nel, ipos; //j,
	int l_dof_ent_types[300], l_dof_ent_ids[300], l_dof_ent_nrdofs[300];
	utt_patches * patches;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	char field_module_name[100];
	apr_module_introduce(field_module_name);
	if ( strncmp(field_module_name, "STANDARD_LINEAR", 15) != 0 && strncmp(field_module_name, "STANDARD_QUADRATIC", 18) ) {
		printf("Patch creation algorithm works only with continuous approximation. Exiting!\n");
		exit(0);
	}

	mesh_id = apr_get_mesh_id(Field_id);

	nno = mmr_get_max_node_id(mesh_id);
	patches = (utt_patches *)malloc((nno + 2) * sizeof(utt_patches));
	nr_patches = nno + 1;
	for (ino = 0; ino <= nr_patches; ino++) patches[ino].deriv = NULL;

	/* initialize patches - central node is always the first !!! */
	ino = 0;
	while ((ino = mmr_get_next_node_all(mesh_id, ino)) != 0) {
		patches[ino].nr_elems = 0;
		for (i = 0; i < UTC_MAXEL_PATCH; i++)
			patches[ino].elems[i] = UTC_LIST_END_MARK;
		for (i = 0; i < UTC_MAXNO_PATCH; i++)
			patches[ino].nodes[i] = UTC_LIST_END_MARK;
		patches[ino].nr_nodes = 1;
		patches[ino].nodes[0] = ino;
	}

	nel = 0;
	while ((nel = mmr_get_next_act_elem(mesh_id, nel)) != 0) {

		int nr_dof_ent_loc = 300; // each time must be initialized!!!
		pdr_comp_stiff_mat(Field_id, PDC_ELEMENT, nel, PDC_NO_COMP, NULL,
				&nr_dof_ent_loc, l_dof_ent_types, l_dof_ent_ids,
				l_dof_ent_nrdofs, NULL, NULL, NULL, NULL);

		for (idofent = 0; idofent < nr_dof_ent_loc; idofent++) {

			if (l_dof_ent_types[idofent] == APC_VERTEX || l_dof_ent_types[idofent] == APC_MIXED_VERTEX)
			{

				dof_ent_id = l_dof_ent_ids[idofent];
				ipos = utr_put_list(nel, patches[dof_ent_id].elems, UTC_MAXEL_PATCH);

#ifdef DEBUG
				if ( ipos == 0 || ipos == -UTC_MAXEL_PATCH) {
					printf("Number of elements %d in patch %d > maximum allowed %d\n",
							ipos, dof_ent_id, UTC_MAXEL_PATCH);
					printf("Recompile code with greater UTC_MAXEL_PATCH. Exiting!\n");
					exit(0);
				}
#endif

				if (ipos < 0) patches[dof_ent_id].nr_elems++;

				for (ino = 0; ino < nr_dof_ent_loc; ino++) {

					nno =  l_dof_ent_ids[ino];
					if (nno != dof_ent_id) {

						ipos = utr_put_list(nno, patches[dof_ent_id].nodes, UTC_MAXNO_PATCH);

#ifdef DEBUG
						if ( ipos == 0 || ipos == -UTC_MAXNO_PATCH) {
							printf("Number of nodes %d in patch %d > maximum allowed %d\n",
									ipos, dof_ent_id, UTC_MAXNO_PATCH);
							printf("Recompile code with greater UTC_MAXNO_PATCH. Exiting!\n");
							exit(0);
						}
#endif
						if (ipos < 0) patches[dof_ent_id].nr_nodes++;

					}
				}
			}
		} // end loop over dof_ents of a given int_ent
	} //end loop over int_ent

#ifdef DEBUG
	nno = 0;
	while ((nno = mmr_get_next_node_all(mesh_id, nno)) != 0) {

		int ino, iel, jno;
		int el_nodes[MMC_MAXELVNO + 1];      /* list of nodes of El */
		int nr_alien = 0;

		for (ino = 0; ino < patches[nno].nr_nodes; ino++) {

			int ifound = 0;

			for (iel = 0; iel < patches[nno].nr_elems; iel++) {

				mmr_el_node_coor(mesh_id, patches[nno].elems[iel], el_nodes, NULL);
				for (jno = 1; jno <= el_nodes[0]; jno++) {

					if (patches[nno].nodes[ino] == el_nodes[jno]) {

						ifound = 1;
						break;

					}

				}

			}

			if (patches[nno].nr_elems > 0 && ifound == 0) {

				printf("Found alien node !!!\n");
				printf("Central node %d, alien node %d\n", nno, patches[nno].nodes[ino]);
				nr_alien++;

			}

		}

		if (patches[nno].nr_elems > 0 && nr_alien > 0) {

			printf("Patch with alien nodes !!!\n");
			//exit(0);

			/*kbw
				  printf("Patch for node %d, nr_elems %d, nr_nodes %d\n",
					 nno, patches[nno].nr_elems,patches[nno].nr_nodes);
				  for(i=0;i<patches[nno].nr_elems;i++)
				printf("%d\t", patches[nno].elems[i]);
				  printf("\n");

				  for(i=0;i<patches[nno].nr_nodes;i++)
				printf("%d\t", patches[nno].nodes[i]);
				  printf("\n");

				  for(iel=0;iel<patches[nno].nr_elems;iel++){
				mmr_el_node_coor(mesh_id,patches[nno].elems[iel],el_nodes,NULL);
				printf("Elem %d nodes: %d %d %d %d\n", patches[nno].elems[iel],
					   el_nodes[1], el_nodes[2], el_nodes[3], el_nodes[4]);
				  }


				  getchar();
			//kew*/


		}

	} // end checking

#endif

	/* nno=0; */
	/* while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){ */

	/* 	//utr_insert_sort(patches[nno].nodes, 0, patches[nno].nr_nodes-1); */

	/* 	/\*kbw */
	/* 	printf("Patch for node %d, nr_elems %d, nr_nodes %d\n",  */
	/* 	   nno, patches[nno].nr_elems,patches[nno].nr_nodes); */
	/* 	for(i=0;i<patches[nno].nr_elems;i++) */
	/* 	  printf("%d\t", patches[nno].elems[i]); */
	/* 	printf("\n"); */

	/* 	for(i=0;i<patches[nno].nr_nodes;i++) */
	/* 	  printf("%d\t", patches[nno].nodes[i]); */
	/* 	printf("\n"); */
	/* 	//kew*\/ */


	/* } */

	*patches_p = patches;
	return (nr_patches);

}



/*---------------------------------------------------------
  utr_recover_derivatives - to recover derivatives using patches of elements
  ---------------------------------------------------------*/
int utr_recover_derivatives( /* returns: >0 - success, <=0 - failure */
		int Field_id,         /* in: field data structure to be used */
		int Sol_vec_id,       /*in: which solution vector to take into account */
		int Nr_patches,
		utt_patches * patches /* in - array of patches for real nodes */
		/* out: array of patches with computed derivatives for ALL nodes */
)
{

	int k, l, nno, iel, nel, ider, mesh_id; //i,j,ino,
	int base;  /* type of basis functions  */
	int old_pdeg, pdeg, geo_pdeg, mixed_pdeg, nreq, mixed_nreq, sol_vec_id, ki, iaux;
	int el_nodes[APC_MAX_GEO_DOFS + 1] = {0};    /* list of nodes of El */
	int el_nodes_type[APC_MAX_GEO_DOFS + 1] = {0};    /* list of nodes of El */
	int el_vertices = 0;        /* number of verices nodes of El */
	double node_coor[3 * APC_MAX_GEO_DOFS] = {0.0}; /* coord of nodes of El */
	double dofs_loc[APC_MAXELSD] = {0.0}; /* element solution dofs */
	double xcoor[3];      /* global coord of gauss point */
	double u_val[UTC_MAXEQ] = {0.0}; /* computed solution */
	double u_x[UTC_MAXEQ] = {0.0}; /* gradient of computed solution */
	double u_y[UTC_MAXEQ] = {0.0}; /* gradient of computed solution */
	double u_z[UTC_MAXEQ] = {0.0}; /* gradient of computed solution */
	double base_phi[APC_MAXELVD] = {0.0};  /* basis functions */
	double base_dphix[APC_MAXELVD] = {0.0}; /* x-derivatives of basis function */
	double base_dphiy[APC_MAXELVD] = {0.0}; /* y-derivatives of basis function */
	double base_dphiz[APC_MAXELVD] = {0.0}; /* y-derivatives of basis function */
	double determ;
	int ngauss;           /* number of gauss points */
	double xg[3000] = {0.0};   	/* coordinates of gauss points in 3D */
	double wg[1000] = {0.0};    /* gauss weights */
	double vol;           /* volume for integration rule */
	double stiff_loc[APC_MAXELVD * APC_MAXELVD] = {0.0};
	/* element stiffness matrix for local problems - single value projection */
	double f_loc[3 * UTC_MAXEQ][APC_MAXELVD] = {{0.0}}; /* rhs vector for local problems */
	/* for each derivative and each component - one vector */

	double * stiff_patch; /* patch stiffness matrix for local problems */
	double ** f_patch;   /* patch rhs vectors for local problems */

	int kk, ieq, idofs, jdofs, num_shap, ndofs_patch, ndofs_el, nr_deriv, ino, i;
	int comp_sm;
	int * list_dof_ent_type;
	/* list of types of dofs for element (here only nodes-vertices) */
	int * list_dof_ent_id;  /* list of nodes' ids for element */
	int * list_dof_ent_nrdofs; /* list of no of dofs for nodes - here ieq or 1 */
	int * tabtrans; /* table to rewrite data for hanging nodes */
	int constr[5]; /* data  to rewrite data for hanging nodes */
	int ione = 1;

	int ipiv[APC_MAXELVD]; /* pivot information for gauss elimination */

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	mesh_id = apr_get_mesh_id(Field_id);

	nreq = apr_get_nreq(Field_id);
	mixed_nreq = apr_get_mixed_nreq(Field_id);

	old_pdeg = -1;

	// in our 3D code number of derivatives is equal 3 for each solution component
	nr_deriv = 3 * nreq;

	for (ino = 1; ino <= Nr_patches; ino++) {
		if (patches[ino].deriv != NULL) free(patches[ino].deriv);
		patches[ino].deriv = (double *)malloc(nr_deriv * sizeof(double));
		for (i = 0; i < nr_deriv; i++) patches[ino].deriv[i] = 0.0;
	}

	/* nno=0; */
	/* while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){ */

	/* 	/\* if(patches[nno].nr_elems != lista[nno].size)  *\/ */
	/* 	/\*   { *\/ */
	/* 	/\* 	printf("error 1 in new patches\n"); exit(-1); *\/ */
	/* 	/\*   } *\/ */
	/* 	/\* for(i=0;i<patches[nno].nr_elems;i++) *\/ */
	/* 	/\*   if(patches[nno].elems[i] != lista[nno].patch[i]) *\/ */
	/* 	/\* 	{ *\/ */
	/* 	/\* 	  printf("error 2 in new patches\n"); exit(-1); *\/ */
	/* 	/\* 	} *\/ */
	/* 	/\* if(patches[nno].nr_nodes != lista[nno].size_stiff_mat)  *\/ */
	/* 	/\*   { *\/ */
	/* 	/\* 	printf("error 3 in new patches\n"); exit(-1); *\/ */
	/* 	/\*   } *\/ */
	/* 	/\* for(i=0;i<patches[nno].nr_nodes;i++) *\/ */
	/* 	/\*   if(patches[nno].nodes[i] != lista[nno].posglob[i]) *\/ */
	/* 	/\* 	{ *\/ */
	/* 	/\* 	  printf("error 4 in new patches\n"); exit(-1); *\/ */
	/* 	/\* 	} *\/ */

	/* 	/\*kbw	  */
	/* 	if(nno==848 || nno==848){  */
	/* 	printf("Patch for node %d, nr_elems %d, nr_nodes %d\n",  */
	/* 	  nno, patches[nno].nr_elems,patches[nno].nr_nodes); */
	/* 	  for(i=0;i<patches[nno].nr_elems;i++) */
	/* 	  printf("%d\t", patches[nno].elems[i]); */
	/* 	  printf("\n"); */

	/* 	  for(i=0;i<patches[nno].nr_nodes;i++) */
	/* 	  printf("%d\t", patches[nno].nodes[i]); */
	/* 	  printf("\n"); */
	/* 	  getchar(); */
	/* 	} */
	/* 	  //kew*\/ */

	/* } */


	// PDC_MAX_DOF_PER_INT - maximal number of dof entities (here called nodes)
	//                       per single element
	// for hexahedra: 8 vertices + 12 edges + 6 faces + 1 interior = 27
	// this should work as well for constrained prisms and tetrahedra
	list_dof_ent_type = (int *)malloc(PDC_MAX_DOF_PER_INT * sizeof(int));
	list_dof_ent_id = (int *)malloc(PDC_MAX_DOF_PER_INT * sizeof(int));
	list_dof_ent_nrdofs = (int *)malloc(PDC_MAX_DOF_PER_INT * sizeof(int));

	nno = 0;
	while ((nno = mmr_get_next_node_all(mesh_id, nno)) != 0) {

		if (patches[nno].nr_elems > 0) {

			// allocate space for patch stiffness matrix - one dof per one node
			ndofs_patch = patches[nno].nr_nodes;
			stiff_patch = (double *)malloc(ndofs_patch * ndofs_patch * sizeof(double));

			for (k = 0; k < ndofs_patch * ndofs_patch; k++) stiff_patch[k] = 0.0;

			// allocate space for right hand sides - one RHS for each component
			// and each derivative

			f_patch = (double **)malloc(nr_deriv * sizeof(double *));
			for (k = 0; k < nr_deriv; k++)
				f_patch[k] = (double *)malloc(ndofs_patch * sizeof(double));

			for (l = 0; l < nr_deriv; l++)
				for (k = 0; k < ndofs_patch; k++)
					f_patch[l][k] = 0.0;

			// assemble patch stiffness matrix and load vectors
			for (iel = 0; iel < patches[nno].nr_elems; iel++) {

				nel = patches[nno].elems[iel];

				/* find degree of polynomial and number of element scalar dofs */
				apr_get_el_pdeg(Field_id, nel, &pdeg);
				apr_get_el_geo_pdeg(Field_id, nel, &geo_pdeg);

				/* get the coordinates of the nodes of El in the right order */
				if (pdeg == APC_QUADRATIC_HIERACHICAL_APPROXIMATION_PDEG) {
					mmr_el_node_coor(mesh_id, nel, el_nodes, node_coor);
					old_pdeg = pdeg;
					pdeg = APC_LINEAR_APPROXIMATION_PDEG;
				} else if (pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
					apr_get_el_geo_dofs(Field_id, nel, el_nodes, el_nodes_type, node_coor);
					el_vertices = 0;
					for (i = 1; i < el_nodes[0]; i++) {
						if (el_nodes_type[i] == APC_VERTEX || el_nodes_type[i] == APC_MIXED_VERTEX) {
							el_vertices++;
						}
					}
					old_pdeg = pdeg;
					pdeg = APC_LINEAR_APPROXIMATION_PDEG;
				} else if (pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
						pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
						pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG) {

					apr_get_mixed_pdeg(Field_id, &mixed_pdeg);
					mixed_nreq = apr_get_mixed_nreq(Field_id);
					if (mixed_pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
						apr_get_el_geo_dofs(Field_id, nel, el_nodes, el_nodes_type, node_coor);
						el_vertices = 0;
						for (i = 1; i < el_nodes[0]; i++) {
							if (el_nodes_type[i] == APC_VERTEX || el_nodes_type[i] == APC_MIXED_VERTEX) {
								el_vertices++;
							}
						}
						old_pdeg = mixed_pdeg;
						mixed_pdeg = APC_LINEAR_APPROXIMATION_PDEG;
					} else {
						apr_get_el_geo_dofs(Field_id, nel, el_nodes, el_nodes_type, node_coor);
					}

				} else {
					mmr_el_node_coor(mesh_id, nel, el_nodes, node_coor);
				}

				// get number of shape functions and element DOFs
				if (old_pdeg == APC_QUADRATIC_HIERACHICAL_APPROXIMATION_PDEG ||
						old_pdeg == APC_QUADRATIC_APPROXIMATION_PDEG ||
						old_pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
						old_pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
						old_pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG)
				{
					num_shap = apr_get_el_pdeg_numshap_related_with_dof_type(Field_id, nel, APC_VERTEX, &old_pdeg);
				}
				else
				{
					num_shap = apr_get_el_pdeg_numshap(Field_id, nel, &pdeg);
				}

				ndofs_el = num_shap;

				/* initialize the matrices to zero */
				for (k = 0; k < ndofs_el * ndofs_el; k++) stiff_loc[k] = 0.0;
				for (l = 0; l < nr_deriv; l++)
					for (k = 0; k < ndofs_el; k++) f_loc[l][k] = 0.0;

				base = apr_get_base_type(Field_id, nel);
				if (base == APC_BASE_PRISM_QUAD_STD) {
					base = APC_BASE_PRISM_STD;
				} else if (base == APC_BASE_TETRA_QUAD_STD) {
					base = APC_BASE_TETRA_STD;
				}

				/* prepare data for gaussian integration */
				if (pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
						pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
						pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG) {
					apr_set_quadr_3D(base, &mixed_pdeg, &ngauss, xg, wg);
				} else {
					apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
				}

				/* get solution degrees of freedom from Sol_vec_id vector*/
				apr_get_el_dofs(Field_id, nel, Sol_vec_id, dofs_loc);

				for (ki = 0; ki < ngauss; ki++) {

					/* at gauss point, compute basis functions, determinant etc*/
					iaux = 2; /* calculations with jacobian but not on the boundary */

					if (pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
							pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
							pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG) {

						determ = apr_elem_calc_3D(iaux, nreq, &mixed_pdeg, base,
										&xg[3 * ki], node_coor, dofs_loc,
										base_phi, base_dphix, base_dphiy, base_dphiz,
										xcoor, u_val, u_x, u_y, u_z, NULL);

					} else {

						determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base,
										&xg[3 * ki], node_coor, dofs_loc,
										base_phi, base_dphix, base_dphiy, base_dphiz,
										xcoor, u_val, u_x, u_y, u_z, NULL);
					}
					vol = determ * wg[ki];

					kk = 0;

					for (jdofs = 0; jdofs < num_shap; jdofs++) {
						for (idofs = 0; idofs < num_shap; idofs++) {

							/* mass matrix */
							stiff_loc[kk + idofs] += base_phi[jdofs] * base_phi[idofs] * vol;

						}/* idofs */
						kk += ndofs_el;

					} /* jdofs */

					for (ieq = 0; ieq < nreq; ieq++) {

						kk = 0;
						for (idofs = 0; idofs < num_shap; idofs++) {
							f_loc[3 * ieq][kk] +=  u_x[ieq] * base_phi[idofs] * vol;
							kk++;
						}
						kk = 0;
						for (idofs = 0; idofs < num_shap; idofs++) {
							f_loc[3 * ieq + 1][kk] += u_y[ieq] * base_phi[idofs] * vol;
							kk++;
						}
						kk = 0;
						for (idofs = 0; idofs < num_shap; idofs++) {
							f_loc[3 * ieq + 2][kk] += u_z[ieq] * base_phi[idofs] * vol;
							kk++;
						}

					}

				} // end loop over Gauss points

				/*kbw
					if(nno==848 || nno==848){
					  printf("\nEl %d, unconstrained nodes:", nel);
					  for(i=1; i<=el_nodes[0]; i++) printf("  %d,",el_nodes[i]);
					  printf("\n");
					  kk=0;
					  printf("\nStiff_mat for el %d:\n",nel);
					  for (jdofs=0;jdofs<ndofs_el;jdofs++)
						{
						  for (idofs=0;idofs<ndofs_el;idofs++)
						{
						  printf(" %15.10lf",stiff_loc[kk+idofs]);
						}
						  kk+=ndofs_el;
						  printf("\n");
						}

					  printf("\nRHS vectors:\n");
					  for(ider=0; ider < nr_deriv; ider++){
						for (idofs=0;idofs<ndofs_el;idofs++)
						  {
						printf(" %15.10lf",f_loc[ider][idofs]);
						  }
						printf("\n");
					  }
					  printf("\n\n");

					  getchar();
					}
				//kew*/

				/* obligatory procedure to fill Lists of dof_ents and rewite SM and RHSV */
				/* nrdofs_loc - allocated size for SM and RHS - checked in get_stiff_mat_data */
				int nrdofs_loc = APC_MAXELVD;
				/* nr_dof_ent - allocated size for list_dof_ent_.. lists - checked in get_stiff_mat_data */
				int nr_dof_ent = PDC_MAX_DOF_PER_INT;

				comp_sm = APC_REWR_BOTH;
				apr_get_stiff_mat_data(Field_id, nel, comp_sm,
						'N', 0, 1, &nr_dof_ent,
						list_dof_ent_type, list_dof_ent_id,
						list_dof_ent_nrdofs,
						&nrdofs_loc, stiff_loc, f_loc[0]);

				comp_sm = APC_REWR_RHS;
				for (ider = 1; ider < nr_deriv; ider++) {

					int nrdofs_loc = APC_MAXELVD; // each time must be initialized!!!
					int nr_dof_ent = PDC_MAX_DOF_PER_INT; // each time must be initialized!!!

					apr_get_stiff_mat_data(Field_id, nel, comp_sm,
							'N', 0, 1, &nr_dof_ent,
							list_dof_ent_type, list_dof_ent_id,
							list_dof_ent_nrdofs,
							&nrdofs_loc, NULL, f_loc[ider]);
				}

				// initialize table for rewriting values for constrained nodes
				tabtrans = (int *)malloc((nrdofs_loc + 1) * sizeof(int));

				for (k = 0; k < nrdofs_loc; k++) {
					int isuccess = 0;
					for (l = 0; l < ndofs_patch; l++) {
						if (list_dof_ent_id[k] == patches[nno].nodes[l]) {
							tabtrans[k] = l;
							isuccess = 1;
						}
					}
					if (isuccess == 0) {
						printf("Node in stiff_mat unconstrained, not present in patch!!! (utr_recover_derivatives). Exiting!!!\n");
						exit(0);
					}
				}

				//assembling
				for (idofs = 0; idofs < nrdofs_loc; idofs++) {
					for (jdofs = 0; jdofs < nrdofs_loc; jdofs++) {
						stiff_patch[tabtrans[idofs] + ndofs_patch * tabtrans[jdofs]]
						+= stiff_loc[idofs + jdofs * nrdofs_loc];
					}
				}

				for (ider = 0; ider < nr_deriv; ider++) {
					for (idofs = 0; idofs < nrdofs_loc; idofs++) {
						f_patch[ider][tabtrans[idofs]] += f_loc[ider][idofs];
					}
				}



				/*kbw
				if(nno==848 || nno==848){
				printf("\nPatch for node %d - %d nodes\n", nno, patches[nno].nr_nodes);

				for(k=0;k<ndofs_patch;k++)
				printf("%d\t", patches[nno].nodes[k]);
				printf("\n");

				printf("nrdofs_loc=%d\nlist_dof_ent_id:\n",nrdofs_loc);
				for (idofs=0;idofs<nrdofs_loc;idofs++) { //for each row!!!!
				printf("%d\t",list_dof_ent_id[idofs]);
				}
				printf("\n");

				printf("Transfer table\n:");
				for (idofs=0;idofs<nrdofs_loc;idofs++)
				printf("%d\t",tabtrans[idofs]);
				printf("\n");

				printf("\nElement %d: Modified stiffness matrix:\n",nel);
				for (idofs=0;idofs<nrdofs_loc;idofs++) { //for each row!!!!
				for (jdofs=0;jdofs<nrdofs_loc;jdofs++) { // for each element in row !!!
				 printf("%15.10lf",stiff_loc[idofs+jdofs*(nrdofs_loc)]);
				}
				printf("\n");
				}
				getchar();

				printf("Element %d: Rhs_vect:\n",nel);
				for(ider=0; ider < nr_deriv; ider++){
				for (idofs=0;idofs<nrdofs_loc;idofs++)
				 printf("%15.10lf",f_loc[ider][idofs]);
				printf("\n");
				}
				printf("\nPatch for node %d stiffness matrix:\n",nno);
				for (idofs=0;idofs<ndofs_patch;idofs++) { //for each row!!!!
				for (jdofs=0;jdofs<ndofs_patch;jdofs++) { // for each element in row !!!
				 printf("%12.8lf",stiff_patch[idofs+jdofs*ndofs_patch]);
				}
				printf("\n");
				}


				getchar();
				}
				//kew*/

				free(tabtrans);

			} //petla po el w patchu


			/*kbw
			if(nno==848 || nno==848){
			printf("\nPatch for node %d stiffness matrix:\n",nno);
			for (idofs=0;idofs<ndofs_patch;idofs++) { //for each row!!!!
			for (jdofs=0;jdofs<ndofs_patch;jdofs++) { // for each element in row !!!
			 printf("%12.8lf",stiff_patch[idofs+jdofs*ndofs_patch]);
			}
			printf("\n");
			}

			printf("Rhs_vect:\n");
			for(ider=0; ider < nr_deriv; ider++){
			for (idofs=0;idofs<ndofs_patch;idofs++)
			 printf("%12.8lf",f_patch[ider][idofs]);
			printf("\n");
			}
			getchar();
			}
			//kew*/

			dgetrf_(&ndofs_patch, &ndofs_patch, stiff_patch, &ndofs_patch, ipiv, &iaux);
			//mf_check(iaux == 0, "DGETRF exited with error %d",iaux);
			if (iaux != 0) {
				printf("\nPatch for node %d, ndofs_patch %d:\n", nno, ndofs_patch);
				printf("utr_recover_derivatives: DGETRF exited with error %d", iaux);
				exit(-1);
			}
			for (ider = 0; ider < nr_deriv; ider++) {
				dgetrs_("N", &ndofs_patch, &ione, stiff_patch, &ndofs_patch, ipiv,
						&f_patch[ider][0], &ndofs_patch, &iaux);
				//mf_check(iaux == 0, "DGETRF exited with error %d",iaux);
				if (iaux != 0) {
					printf("utr_recover_derivatives: derivative %d, DGETRS exited with error %d",
							ider, iaux);
					exit(-1);
				}
			}

			/*kbw
					printf("\nRecovered derivatives for node %d :\n",nno);
				  for(ider=0;ider<nr_deriv;ider++)
				{
				  for (jdofs=0;jdofs<ndofs_patch;jdofs++)
					{
					  printf("%8.4lf",f_patch[ider][jdofs]);
					}
				  printf("\n");
				}
			//kew*/

			/* according to the convention central node is always the first !!! */
			for (ider = 0; ider < nr_deriv; ider++) {
				patches[nno].deriv[ider] = f_patch[ider][0];
			}

			free(stiff_patch);
			for (k = 0; k < nr_deriv; k++) free(f_patch[k]);
			free(f_patch);

		} // end if patch not empty - real node
	} // end loop over patches

	/* having computed derivatives for real nodes, */
	/* now is the time for hanging nodes */
	nno = 0;
	while ((nno = mmr_get_next_node_all(mesh_id, nno)) != 0) {

		if (patches[nno].nr_elems == 0) {

			apr_get_constr_data(Field_id, nno, APC_VERTEX, constr, NULL);

			/*kbw
				  if(nno==848 || nno==848){
				printf("node %d is constrained by %d nodes:\n",nno,constr[0]);
				for(i=1;i<=constr[0];i++)
				  printf("%d\t",constr[i]);
				printf("\n");
				  }
			//kew*/

			for (jdofs = 1; jdofs <= constr[0]; jdofs++) {

				for (ider = 0; ider < 3; ider++) {

					patches[nno].deriv[ider] +=
							1. / constr[0] * patches[constr[jdofs]].deriv[ider];

				}
			}

			/*kbw
				  if(nno==848 || nno==848){
				double f,f_x,f_y,f_z,xcoor[3],daux,eaux; int j;
				printf("\nAt constrained node %d :\n",nno);
				printf("Real derivatives:\t");
				mmr_node_coor(mesh_id, nno, xcoor);
				pdr_exact_sol(iaux,xcoor[0],xcoor[1],xcoor[2],daux,
						  &f,&f_x,&f_y,&f_z,&eaux);
				printf("%8.4lf%8.4lf%8.4lf", f_x, f_y, f_z);
				printf("\nRecovered derivatives:\t");
				for(ider=0;ider<nr_deriv;ider++)
				  {
					printf("%8.4lf", patches[nno].deriv[ider]);
				  }
				printf("\n");
				for(i=1;i<=constr[0];i++){
				  printf("\nAt constraining node %d :\n",constr[i]);
				  printf("Real derivatives:\t");
				  mmr_node_coor(mesh_id, constr[i], xcoor);
				  pdr_exact_sol(iaux,xcoor[0],xcoor[1],xcoor[2],daux,
						&f,&f_x,&f_y,&f_z,&eaux);
				  printf("%8.4lf%8.4lf%8.4lf", f_x, f_y, f_z);
				  printf("\nRecovered derivatives:\t");
				  for(ider=0;ider<nr_deriv;ider++)
					{
					  printf("%8.4lf", patches[constr[i]].deriv[ider]);
					}

				}

				printf("\n");
				getchar();
				  }
			//kew*/

		} // end if patch empty - constrained node
	} // end loop over patches

	/*kbw
	  nno=0;
	  while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){
		if(nno==848 || nno==848){
		printf("\nRecovered derivatives for node %d :\n",nno);
		for(ider=0;ider<nr_deriv;ider++)
		  {
		printf("%8.4lf", patches[nno].deriv[ider]);
		  }
		printf("\n");
		getchar();
		}
	  }
	//kew*/

	free(list_dof_ent_type);
	free(list_dof_ent_id);
	free(list_dof_ent_nrdofs);

	return (1);
}

/*---------------------------------------------------------
 utr_create_patches_small - to create patches of elements containing a given node
  (uses lists of elements nodes instead of lists of elements dof entities)
---------------------------------------------------------*/
int utr_create_patches_small( /* returns: >0 - Nr_patches, <=0 - failure */
		int Field_id,       /* in: field data structure to be used */
		utt_patches ** patches_p /* in/out - array of patches for real nodes */
)
{
	int mesh_id, nreq, nr_patches;
	int i, ino, nno, jno, node_id, nel, nr_dof_ent_loc, ipos; //j,
	//int l_dof_ent_types[300], l_dof_ent_ids[300], l_dof_ent_nrdofs[300];
	int el_nodes[MMC_MAXELVNO + 1];      /* list of nodes of El */
	utt_patches * patches;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	mesh_id = apr_get_mesh_id(Field_id);

	nno = mmr_get_max_node_id(mesh_id);
	patches = (utt_patches *)malloc((nno + 2) * sizeof(utt_patches));
	nr_patches = nno + 1;
	for (ino = 0; ino <= nr_patches; ino++) patches[ino].deriv = NULL;
	for (ino = 1; ino <= nr_patches; ino++) { //!!!ino=0
		// hanging nodes have nr_elems > 0 (added below - different than in standard version)
		// real nodes have nr_elems > 0 (added below)
		patches[ino].nr_elems = 0;
		// free spaces are identified by nr_elems = -1
		if (mmr_node_status(mesh_id, ino) == MMC_FREE) patches[ino].nr_elems = -1;
	}

	/* initialize patches - central node is always the first !!! */
	ino = 0;
	while ((ino = mmr_get_next_node_all(mesh_id, ino)) != 0) {
		// patches[ino].nr_elems=0; - already done above
		for (i = 0; i < UTC_MAXEL_PATCH; i++)
			patches[ino].elems[i] = UTC_LIST_END_MARK;
		for (i = 0; i < UTC_MAXNO_PATCH; i++)
			patches[ino].nodes[i] = UTC_LIST_END_MARK;
		patches[ino].nr_nodes = 1;
		patches[ino].nodes[0] = ino;
	}

	nel = 0;
	while ((nel = mmr_get_next_act_elem(mesh_id, nel)) != 0) {

		// we get vertices - both real and hanging nodes
		mmr_el_node_coor(mesh_id, nel, el_nodes, NULL);

		//pdr_comp_stiff_mat(Field_id, PDC_ELEMENT,nel, PDC_NO_COMP, NULL,
		//		       &nr_dof_ent_loc, l_dof_ent_types, l_dof_ent_ids,
		//		       l_dof_ent_nrdofs, NULL, NULL, NULL, NULL);

		//for(idofent=0; idofent<nr_dof_ent_loc; idofent++){

		for (jno = 1; jno <= el_nodes[0]; jno++) {
			node_id = el_nodes[jno];
			ipos = utr_put_list(nel, patches[node_id].elems, UTC_MAXEL_PATCH);

#ifdef DEBUG
			if ( ipos == 0 || ipos == -UTC_MAXEL_PATCH) {
				printf("Number of elements %d in patch %d > maximum allowed %d\n",
						ipos, node_id, UTC_MAXEL_PATCH);
				printf("Recompile code with greater UTC_MAXEL_PATCH. Exiting!\n");
				exit(0);
			}
#endif

			if (ipos < 0) patches[node_id].nr_elems++;

			for (ino = 1; ino <= el_nodes[0]; ino++) {

				nno =  el_nodes[ino];
				if (nno != node_id) {

					ipos = utr_put_list(nno, patches[node_id].nodes, UTC_MAXNO_PATCH);

#ifdef DEBUG
					if ( ipos == 0 || ipos == -UTC_MAXNO_PATCH) {
						printf("Number of nodes %d in patch %d > maximum allowed %d\n",
								ipos, node_id, UTC_MAXNO_PATCH);
						printf("Recompile code with greater UTC_MAXNO_PATCH. Exiting!\n");
						exit(0);
					}
#endif
					if (ipos < 0) patches[node_id].nr_nodes++;

				}
			}

		} // end loop over vertices (nodes) of a given element
	} //end loop over elements

#ifdef DEBUG
	nno = 0;
	while ((nno = mmr_get_next_node_all(mesh_id, nno)) != 0) {

		int ino, iel, jno;
		int nr_alien = 0;

		for (ino = 0; ino < patches[nno].nr_nodes; ino++) {

			int ifound = 0;

			for (iel = 0; iel < patches[nno].nr_elems; iel++) {

				mmr_el_node_coor(mesh_id, patches[nno].elems[iel], el_nodes, NULL);
				for (jno = 1; jno <= el_nodes[0]; jno++) {

					if (patches[nno].nodes[ino] == el_nodes[jno]) {

						ifound = 1;
						break;

					}

				}

			}

			if (patches[nno].nr_elems > 0 && ifound == 0) {

				printf("Found alien node !!!\n");
				printf("Central node %d, alien node %d\n", nno, patches[nno].nodes[ino]);
				nr_alien++;

			}

		}

		if (patches[nno].nr_elems > 0 && nr_alien > 1) {

			printf("Patch with alien nodes !!!\n");
			//exit(0);

			/*kbw
				  printf("Patch for node %d, nr_elems %d, nr_nodes %d\n",
					 nno, patches[nno].nr_elems,patches[nno].nr_nodes);
				  for(i=0;i<patches[nno].nr_elems;i++)
				printf("%d\t", patches[nno].elems[i]);
				  printf("\n");

				  for(i=0;i<patches[nno].nr_nodes;i++)
				printf("%d\t", patches[nno].nodes[i]);
				  printf("\n");

				  for(iel=0;iel<patches[nno].nr_elems;iel++){
				mmr_el_node_coor(mesh_id,patches[nno].elems[iel],el_nodes,NULL);
				printf("Elem %d nodes: %d %d %d %d\n", patches[nno].elems[iel],
					   el_nodes[1], el_nodes[2], el_nodes[3], el_nodes[4]);
				  }


				  getchar();
			//kew*/


		}

	}
#endif

	/* nno=0; */
	/* while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){ */

	/*   //utr_insert_sort(patches[nno].nodes, 0, patches[nno].nr_nodes-1); */

	/* 	/\*kbw */
	/* 	printf("Patch for node %d, nr_elems %d, nr_nodes %d\n",  */
	/* 	   nno, patches[nno].nr_elems,patches[nno].nr_nodes); */
	/* 	for(i=0;i<patches[nno].nr_elems;i++) */
	/* 	  printf("%d\t", patches[nno].elems[i]); */
	/* 	printf("\n"); */

	/* 	for(i=0;i<patches[nno].nr_nodes;i++) */
	/* 	  printf("%d\t", patches[nno].nodes[i]); */
	/* 	printf("\n"); */
	/* 	//kew*\/ */


	/* } */

	*patches_p = patches;
	return (nr_patches);

}

/*---------------------------------------------------------
  utr_recover_derivatives_small-to recover derivatives using patches of elements
	(does not use apr_get_stiff_mat_data for getting lists of element nodes)
(procedure creates systems for all nodes, including hanging nodes)
  ---------------------------------------------------------*/
int utr_recover_derivatives_small( /* returns: >0 - success, <=0 - failure */
		int Field_id,         /* in: field data structure to be used */
		int Sol_vec_id,       /*in: which solution vector to take into account */
		int Nr_patches,
		utt_patches * patches /* in - array of patches for real nodes */
		/* out: array of patches with computed derivatives for ALL nodes */
)
{

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	{	// scope of variables

		int mesh_id = apr_get_mesh_id(Field_id);

		int nreq = apr_get_nreq(Field_id);
		int mixed_nreq = apr_get_mixed_nreq(Field_id);

		// in our 3D code number of derivatives is equal 3 for each solution component
		int nr_deriv = 3 * nreq;

		int ino;
		int i;
		for (ino = 1; ino <= Nr_patches; ino++) {
			if (patches[ino].deriv != NULL) free(patches[ino].deriv);
			patches[ino].deriv = (double *)malloc(nr_deriv * sizeof(double));
			for (i = 0; i < nr_deriv; i++) patches[ino].deriv[i] = 0.0;

			//?????????????????????????????????????????????????????????????
			/* if(patches[ino].nr_elems == 0 && ino != Nr_patches){ */
			/*   printf("patch %d with zero elements in small patches\n", ino); exit(-1); */
			/* } */

#ifdef DEBUG
			if (patches[ino].nr_elems < 0) {
				printf("free space %d with -1 elements in small patches\n", ino);
			}
			/*kbw
			if(patches[ino].nr_elems > 0){
			  printf("normal %d patch with %d elements in small patches\n",
				 ino, patches[ino].nr_elems);
			}
			//kew*/
#endif

		}

	}
	// Those lines below are here because of some strange openmp compiler behavior in win env.
	//
	char N_string[] = "N";
	char str1[] = "\nPatch for node %d, ndofs_patch %d:\n";
	char str2[] = "utr_recover_derivatives_small: DGETRF exited with error %d";
	char str3[] = "utr_recover_derivatives_small: derivative %d, DGETRS exited with error %d";
#if WIN32
	#pragma omp parallel default(none) firstprivate(Field_id, Nr_patches, Sol_vec_id, patches) shared(N_string,str1,str2,str3)
#else
	#pragma omp parallel default(none) firstprivate(Field_id, Nr_patches, Sol_vec_id, patches, N_string) shared(utv_log_out,str1,str2,str3)
#endif // WIN32
	{

		int k, l, nno, iel, nel, ider, mesh_id; //i,j,ino,
		int base;  /* type of basis functions  */
		int old_pdeg, pdeg, geo_pdeg, mixed_pdeg, nreq, mixed_nreq, sol_vec_id, ki, iaux;
		int el_nodes[APC_MAX_GEO_DOFS + 1];      /* list of nodes of El */
		int el_vertices = 0;   /* number of vertices nodes of El */
		int el_nodes_type[APC_MAX_GEO_DOFS + 1]; /* list of nodes type of El */
		double node_coor[3 * APC_MAX_GEO_DOFS]; /* coord of nodes of El */
		double dofs_loc[APC_MAXELSD]; /* element solution dofs */
		double xcoor[3];      /* global coord of gauss point */
		double u_val[UTC_MAXEQ]; /* computed solution */
		double u_x[UTC_MAXEQ];   /* gradient of computed solution */
		double u_y[UTC_MAXEQ];   /* gradient of computed solution */
		double u_z[UTC_MAXEQ];   /* gradient of computed solution */
		double base_phi[APC_MAXELVD];    /* basis functions */
		double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
		double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
		double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
		double determ;
		int ngauss;           /* number of gauss points */
		double xg[3000];   	/* coordinates of gauss points in 3D */
		double wg[1000];      /* gauss weights */
		double vol;           /* volume for integration rule */
		double stiff_loc[APC_MAXELVD * APC_MAXELVD];
		/* element stiffness matrix for local problems - single value projection */
		double f_loc[3 * UTC_MAXEQ][APC_MAXELVD]; /* rhs vector for local problems */
		/* for each derivative and each component - one vector */

		double * stiff_patch; /* patch stiffness matrix for local problems */
		double ** f_patch;   /* patch rhs vectors for local problems */

		int kk, ieq, idofs, jdofs, num_shap, ndofs_patch, ndofs_el, nr_deriv, ino, i;
		int comp_sm;
		int nr_dof_ent;
		int nrdofs_loc = 24;
		int * list_dof_ent_type;
		/* list of types of dofs for element (here only nodes-vertices) */
		int * list_dof_ent_id;  /* list of nodes' ids for element */
		int * list_dof_ent_nrdofs; /* list of no of dofs for nodes - here ieq or 1 */
		int * tabtrans; /* table to rewrite data for hanging nodes */
		int constr[5]; /* data  to rewrite data for hanging nodes */
		int ione = 1;

		int ipiv[APC_MAXELVD]; /* pivot information for gauss elimination */


		mesh_id = apr_get_mesh_id(Field_id);
		nreq = apr_get_nreq(Field_id);
		mixed_nreq = apr_get_mixed_nreq(Field_id);

		old_pdeg = -1;

		// in our 3D code number of derivatives is equal 3 for each solution component
		nr_deriv = 3 * nreq;

		// PDC_MAX_DOF_PER_INT - maximal number of dof entities (here called nodes)
		//                       per single element
		// for hexahedra: 8 vertices + 12 edges + 6 faces + 1 interior = 27
		list_dof_ent_type = (int *)malloc(PDC_MAX_DOF_PER_INT * sizeof(int));
		list_dof_ent_id = (int *)malloc(PDC_MAX_DOF_PER_INT * sizeof(int));
		list_dof_ent_nrdofs = (int *)malloc(PDC_MAX_DOF_PER_INT * sizeof(int));


		// old version - iterator is not parallelized
		//nno=0;
		//while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){

		#pragma omp for
		for (nno = 1; nno <= Nr_patches; nno++) {

			if (patches[nno].nr_elems > 0) {

				// allocate space for patch stiffness matrix - one dof per one node
				ndofs_patch = patches[nno].nr_nodes;
				stiff_patch = (double *)malloc(ndofs_patch * ndofs_patch * sizeof(double));

				for (k = 0; k < ndofs_patch * ndofs_patch; k++) stiff_patch[k] = 0.0;

				// allocate space for right hand sides - one RHS for each component
				// and each derivative

				f_patch = (double **)malloc(nr_deriv * sizeof(double *));
				for (k = 0; k < nr_deriv; k++)
					f_patch[k] = (double *)malloc(ndofs_patch * sizeof(double));

				for (l = 0; l < nr_deriv; l++)
					for (k = 0; k < ndofs_patch; k++)
						f_patch[l][k] = 0.0;

				// assemble patch stiffness matrix and load vectors
				for (iel = 0; iel < patches[nno].nr_elems; iel++) {

					nel = patches[nno].elems[iel];

					/* find degree of polynomial and number of element scalar dofs */
					apr_get_el_pdeg(Field_id, nel, &pdeg);
					apr_get_el_geo_pdeg(Field_id, nel, &geo_pdeg);

					/* get the coordinates of the nodes of El in the right order */
					if (pdeg == APC_QUADRATIC_HIERACHICAL_APPROXIMATION_PDEG) {
						mmr_el_node_coor(mesh_id, nel, el_nodes, node_coor);
						old_pdeg = pdeg;
						pdeg = APC_LINEAR_APPROXIMATION_PDEG;
					} else if (pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
						apr_get_el_geo_dofs(Field_id, nel, el_nodes, el_nodes_type, node_coor);
						el_vertices = 0;
						for (i = 1; i < el_nodes[0]; i++) {
							if (el_nodes_type[i] == APC_VERTEX || el_nodes_type[i] == APC_MIXED_VERTEX) {
								el_vertices++;
							}
						}
						old_pdeg = pdeg;
						pdeg = APC_LINEAR_APPROXIMATION_PDEG;
					} else if (pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
							pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
							pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG) {

						apr_get_mixed_pdeg(Field_id, &mixed_pdeg);
						mixed_nreq = apr_get_mixed_nreq(Field_id);
						if (mixed_pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
							apr_get_el_geo_dofs(Field_id, nel, el_nodes, el_nodes_type, node_coor);
							el_vertices = 0;
							for (i = 1; i < el_nodes[0]; i++) {
								if (el_nodes_type[i] == APC_VERTEX || el_nodes_type[i] == APC_MIXED_VERTEX) {
									el_vertices++;
								}
							}
							old_pdeg = mixed_pdeg;
							mixed_pdeg = APC_LINEAR_APPROXIMATION_PDEG;
						} else {
							apr_get_el_geo_dofs(Field_id, nel, el_nodes, el_nodes_type, node_coor);
						}

					} else {
						mmr_el_node_coor(mesh_id, nel, el_nodes, node_coor);
					}

					// get number of shape functions and element DOFs
					//num_shap = apr_get_el_pdeg_numshap(Field_id, nel, &pdeg);
					if (old_pdeg == APC_QUADRATIC_HIERACHICAL_APPROXIMATION_PDEG ||
							old_pdeg == APC_QUADRATIC_APPROXIMATION_PDEG ||
							old_pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
							old_pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
							old_pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG)
					{
						num_shap = apr_get_el_pdeg_numshap_related_with_dof_type(Field_id, nel, APC_VERTEX, &old_pdeg);
					}
					else
					{
						num_shap = apr_get_el_pdeg_numshap(Field_id, nel, &pdeg);
					}
					ndofs_el = num_shap;

					/* initialize the matrices to zero */
					for (k = 0; k < ndofs_el * ndofs_el; k++) stiff_loc[k] = 0.0;
					for (l = 0; l < nr_deriv; l++)
						for (k = 0; k < ndofs_el; k++) f_loc[l][k] = 0.0;

					base = apr_get_base_type(Field_id, nel);
					if (base == APC_BASE_PRISM_QUAD_STD) {
						base = APC_BASE_PRISM_STD;
					} else if (base == APC_BASE_TETRA_QUAD_STD) {
						base = APC_BASE_TETRA_STD;
					}

					/* prepare data for gaussian integration */
					if (pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
							pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
							pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG) {
						apr_set_quadr_3D(base, &mixed_pdeg, &ngauss, xg, wg);
					} else {
						apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
					}

					/* get the indicated solution degrees of freedom */
					apr_get_el_dofs(Field_id, nel, Sol_vec_id, dofs_loc);

					for (ki = 0; ki < ngauss; ki++) {

						/* at gauss point, compute basis functions, determinant etc*/
						iaux = 2; /* calculations with jacobian but not on the boundary */

						if (pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
								pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
								pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG) {

							determ = apr_elem_calc_3D(iaux, nreq, &mixed_pdeg, base,
											&xg[3 * ki], node_coor, dofs_loc,
											base_phi, base_dphix, base_dphiy, base_dphiz,
											xcoor, u_val, u_x, u_y, u_z, NULL);
						} else {

							determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base,
											&xg[3 * ki], node_coor, dofs_loc,
											base_phi, base_dphix, base_dphiy, base_dphiz,
											xcoor, u_val, u_x, u_y, u_z, NULL);
						}

						/*kbw
						if(ki==0){
						  printf("Patch element %d (%d) : node_coor[0] = %lf \n", iel, nel, node_coor[0]);
						  for(ieq=0;ieq<nreq;ieq++){
							printf("%d unknown: value %lf, derivatives: x = %lf, ,y = %lf, ,z = %lf\n",
						ieq, u_val[ieq], u_x[ieq], u_y[ieq], u_z[ieq]);

						  }
						}
						/*kew*/

						vol = determ * wg[ki];

						kk = 0;

						for (jdofs = 0; jdofs < num_shap; jdofs++) {
							for (idofs = 0; idofs < num_shap; idofs++) {

								/* mass matrix */
								stiff_loc[kk + idofs] += base_phi[jdofs] * base_phi[idofs] * vol;

							}/* idofs */
							kk += ndofs_el;

						} /* jdofs */

						for (ieq = 0; ieq < nreq; ieq++) {

							kk = 0;
							for (idofs = 0; idofs < num_shap; idofs++) {
								f_loc[3 * ieq][kk] +=  u_x[ieq] * base_phi[idofs] * vol;
								kk++;
							}
							kk = 0;
							for (idofs = 0; idofs < num_shap; idofs++) {
								f_loc[3 * ieq + 1][kk] += u_y[ieq] * base_phi[idofs] * vol;
								kk++;
							}
							kk = 0;
							for (idofs = 0; idofs < num_shap; idofs++) {
								f_loc[3 * ieq + 2][kk] += u_z[ieq] * base_phi[idofs] * vol;
								kk++;
							}

						}

					} // end loop over Gauss points

					/*kbw
						if(nno==848 || nno==848){
						  printf("\nEl %d, unconstrained nodes:", nel);
						  for(i=1; i<=el_nodes[0]; i++) printf("  %d,",el_nodes[i]);
						  printf("\n");
						  kk=0;
						  printf("\nStiff_mat for el %d:\n",nel);
						  for (jdofs=0;jdofs<ndofs_el;jdofs++)
							{
							  for (idofs=0;idofs<ndofs_el;idofs++)
							{
							  printf(" %15.10lf",stiff_loc[kk+idofs]);
							}
							  kk+=ndofs_el;
							  printf("\n");
							}

						  printf("\nRHS vectors:\n");
						  for(ider=0; ider < nr_deriv; ider++){
							for (idofs=0;idofs<ndofs_el;idofs++)
							  {
							printf(" %15.10lf",f_loc[ider][idofs]);
							  }
							printf("\n");
						  }
						  printf("\n\n");

						  getchar();
						}
					//kew*/

					// initialize table for rewriting values for constrained nodes
					if (geo_pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
						nrdofs_loc = el_vertices;
					} else {
						nrdofs_loc = el_nodes[0];
					}

					tabtrans = (int *)malloc((nrdofs_loc + 1) * sizeof(int));

					if (geo_pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
						for (k = 1; k <= el_vertices; k++) {
							for (l = 0; l < ndofs_patch; l++) {
								if (el_nodes[k] == patches[nno].nodes[l])
									tabtrans[k - 1] = l;
							}
						}
					} else {
						for (k = 1; k <= el_nodes[0]; k++) {
							for (l = 0; l < ndofs_patch; l++) {
								if (el_nodes[k] == patches[nno].nodes[l])
									tabtrans[k - 1] = l;
							}
						}
					}

					//assembling
					for (idofs = 0; idofs < nrdofs_loc; idofs++) {
						for (jdofs = 0; jdofs < nrdofs_loc; jdofs++) {
							stiff_patch[tabtrans[idofs] + ndofs_patch * tabtrans[jdofs]]
							+= stiff_loc[idofs + jdofs * nrdofs_loc];
						}
					}

					for (ider = 0; ider < nr_deriv; ider++) {
						for (idofs = 0; idofs < nrdofs_loc; idofs++) {
							f_patch[ider][tabtrans[idofs]] += f_loc[ider][idofs];
						}
					}



					/*kbw
					if(nno==848 || nno==848){
					printf("\nPatch for node %d - %d nodes\n", nno, patches[nno].nr_nodes);

					for(k=0;k<ndofs_patch;k++)
					printf("%d\t", patches[nno].nodes[k]);
					printf("\n");

					printf("nrdofs_loc=%d\nlist_dof_ent_id:\n",nrdofs_loc);
					for (idofs=0;idofs<nrdofs_loc;idofs++) { //for each row!!!!
					printf("%d\t",list_dof_ent_id[idofs]);
					}
					printf("\n");

					printf("Transfer table\n:");
					for (idofs=0;idofs<nrdofs_loc;idofs++)
					printf("%d\t",tabtrans[idofs]);
					printf("\n");

					printf("\nElement %d: Modified stiffness matrix:\n",nel);
					for (idofs=0;idofs<nrdofs_loc;idofs++) { //for each row!!!!
					for (jdofs=0;jdofs<nrdofs_loc;jdofs++) { // for each element in row !!!
					 printf("%15.10lf",stiff_loc[idofs+jdofs*(nrdofs_loc)]);
					}
					printf("\n");
					}
					getchar();

					printf("Element %d: Rhs_vect:\n",nel);
					for(ider=0; ider < nr_deriv; ider++){
					for (idofs=0;idofs<nrdofs_loc;idofs++)
					 printf("%15.10lf",f_loc[ider][idofs]);
					printf("\n");
					}
					printf("\nPatch for node %d stiffness matrix:\n",nno);
					for (idofs=0;idofs<ndofs_patch;idofs++) { //for each row!!!!
					for (jdofs=0;jdofs<ndofs_patch;jdofs++) { // for each element in row !!!
					 printf("%12.8lf",stiff_patch[idofs+jdofs*ndofs_patch]);
					}
					printf("\n");
					}


					getchar();
					}
					//kew*/

					free(tabtrans);

				} //petla po el w patchu


				/*kbw
					if(nno==848 || nno==848){
				 printf("\nPatch for node %d stiffness matrix:\n",nno);
				 for (idofs=0;idofs<ndofs_patch;idofs++) { //for each row!!!!
				   for (jdofs=0;jdofs<ndofs_patch;jdofs++) { // for each element in row !!!
					 printf("%12.8lf",stiff_patch[idofs+jdofs*ndofs_patch]);
				   }
				   printf("\n");
				 }

				 printf("Rhs_vect:\n");
				 for(ider=0; ider < nr_deriv; ider++){
				   for (idofs=0;idofs<ndofs_patch;idofs++)
					 printf("%12.8lf",f_patch[ider][idofs]);
				   printf("\n");
				 }
				 getchar();
					}
				//kew*/


				dgetrf_(&ndofs_patch, &ndofs_patch, stiff_patch, &ndofs_patch, ipiv, &iaux);
				if (iaux != 0) {
					printf(str1, nno, ndofs_patch);
					printf(str2, iaux);
					exit(-1);
				}
				for (ider = 0; ider < nr_deriv; ider++) {
					dgetrs_(N_string, &ndofs_patch, &ione, stiff_patch, &ndofs_patch, ipiv,
							&f_patch[ider][0], &ndofs_patch, &iaux);
					if (iaux != 0) {

						printf(str3,
								ider, iaux);
						exit(-1);
					}
				}

				/*kbw
				//KBAS
					  printf("\nRecovered derivatives for node %d :\n",nno);
					  for(ider=0;ider<nr_deriv;ider++)
					{
					  for (jdofs=0;jdofs<ndofs_patch;jdofs++)
						{
						  printf("%8.4lf",f_patch[ider][jdofs]);
						}
					  printf("\n");
					}
					getchar();
				//kew*/

				/* according to the convention central node is always the first !!! */
				for (ider = 0; ider < nr_deriv; ider++) {
					patches[nno].deriv[ider] = f_patch[ider][0];
				}

				free(stiff_patch);
				for (k = 0; k < nr_deriv; k++) free(f_patch[k]);
				free(f_patch);

			} // end if patch not empty - real node
		} // end loop over patches

		free(list_dof_ent_type);
		free(list_dof_ent_id);
		free(list_dof_ent_nrdofs);

	} // end parallel region


	return (1);
}


enum UTE_ADAPT_POLICY { DEREF_ONLY = -1, REF_ONLY = 1};
/*---------------------------------------------------------
utr_adapt_core - to perform actual refinement/unrefinement
---------------------------------------------------------*/
int utr_adapt_core(const int Problem_id,    // in
		const char * Work_dir,   // in
		const int * Elems,       // in
		const int   nElems,      // in
		const int  Adapt_policy, // in
		const FILE * Interactive_input,  // in
		FILE * Interactive_output, // in
		const int  Print_info)   // in
{

	if ( nElems < 1 ) {
		mf_check_info(Elems != NULL, "Empty array passed to utr_adapt_core!");
		return 0;
	}


	//time_init();

	// This is needed to allow reallocation this array on heap.
	int * elems = (int *) calloc(nElems, sizeof(int));
	mf_check_mem(elems);
	memcpy(elems, Elems, nElems * sizeof(int));


	int nrno, nmno, nred, nmed, nrfa, nmfa, nrel, nmel;
	int field_id, mesh_id, i, result = 0, massive_indicator = elems[0];

	i = 3;
	field_id = pdr_ctrl_i_params(Problem_id, i);
	mesh_id = apr_get_mesh_id(field_id);


	//utr_io_write_mesh_info_to_PAM(mesh_id,Work_dir,"elems_gen_prev","test comment",Interactive_output);

	nrno = mmr_get_nr_node(mesh_id);
	nmno = mmr_get_max_node_id(mesh_id);
	nred = mmr_get_nr_edge(mesh_id);
	nmed = mmr_get_max_edge_id(mesh_id);
	nrfa = mmr_get_nr_face(mesh_id);
	nmfa = mmr_get_max_face_id(mesh_id);
	nrel = mmr_get_nr_elem(mesh_id);
	nmel = mmr_get_max_elem_id(mesh_id);

#ifdef PARALLEL
	if (pcr_my_proc_id() == pcr_print_master()) {
#endif


		if (Print_info >= 1) {
			utr_print_adapt_info(field_id, Interactive_output);
		}

		if (Elems[0] < 0) {
			mf_log_info("utr_adapt_core:Uniform refinement");
		}
		else if (Adapt_policy == DEREF_ONLY) {
			mf_log_info("utr_adapt_core: number of elems to derefine = %d.", nElems);
		}
		else if (Adapt_policy == REF_ONLY) {
			mf_log_info("utr_adapt_core: number of elems to refine = %d.", nElems);
		}
		else {
			mf_log_info("utr_adapt_core: Unknown refinement policy!");
		}

#ifdef PARALLEL
	}

	mmpr_init_ref(mesh_id);
#endif



	mmr_init_ref(mesh_id);

#ifdef PARALLEL
	{
		int prev_nElems = 0;
		do {
			mfp_debug("----- New iteration of updating ref list (size=%d)!", nElems);
			prev_nElems = nElems;
			mmpr_update_ref_list(mesh_id, &nElems, &elems);
		} while (prev_nElems != nElems);
		mfp_debug("----- End of updating ref list (size=%d)!", nElems);
	}
#endif

	switch (Adapt_policy) {
		case DEREF_ONLY:
			for (i = 0; i < nElems; ++i) {
				int can_deref = APC_DEREF_ALLOWED;
				if (elems[i] > 0) {
					if (MMC_ACTIVE != mmr_el_status(mesh_id, elems[i])
							|| APC_DEREF_ALLOWED != apr_limit_deref(field_id, elems[i]) ) {
						can_deref = APC_DEREF_DENIED;
					}
				}

				if (can_deref == APC_DEREF_ALLOWED) {
					mfp_debug("\nutr_adapt_core: derefining elem %d", elems[i]);

					result = apr_derefine(field_id, elems[i]);

				}
				else { // APC_DEREF_DENIED
					mfp_debug("\nutr_adapt_core: derefining elem %d denied!", elems[i]);
					result = -2;
				}
			}
			break;

		case REF_ONLY:
			for (i = 0; i < nElems; ++i) {
				int can_ref = APC_REF_ALLOWED;
				if (elems[i] > 0) {
					if (MMC_ACTIVE != mmr_el_status(mesh_id, elems[i])
							|| APC_REF_ALLOWED != apr_limit_ref(field_id, elems[i])) {
						can_ref = APC_REF_DENIED;
					}
				}

				if ( can_ref == APC_REF_ALLOWED ) {
					mfp_debug("\nutr_adapt_core: refining elem %d", elems[i]);

					result = apr_refine(field_id, elems[i]);

				}
				else { //APC_REF_DENIED
					mfp_debug("\nutr_adapt_core: refining elem %d denied!", elems[i]);
					result = -1;
				}
			}
			break;

		default:
			mf_log_info("\nUknown adaptation policy, ignoring.");
			result = -3;
			break;
	}


	if ((nElems > 1) || (elems[0] < 0)) {
		massive_indicator = -3;
	}


#ifdef PARALLEL
	if (mmpr_is_ready_for_proj_dof_ref(mesh_id) < 0) {
		mf_fatal_err("Parallel mesh not ready for dofs projection!");
	}
#endif

	apr_proj_dof_ref(field_id, massive_indicator, nmel, nmfa, nmed, nmno);

	mmr_final_ref(mesh_id);
	apr_check_field(field_id);

#ifdef PARALLEL
	mmpr_final_ref(mesh_id);

	time_print();

	if (pcr_my_proc_id() == pcr_print_master()) {
#endif
		if (Print_info >= 1) {
			utr_print_adapt_info(field_id, Interactive_output);
		}
#ifdef PARALLEL
	}
#endif

	//utr_io_write_mesh_info_to_PAM(mesh_id,Work_dir,"elems_gen","test comment",Interactive_output);

	free(elems);

	return result;
}


/*---------------------------------------------------------
utr_adapt - to enforce default adaptation strategy based on
		element error indicators provided by pdr_err_indi
---------------------------------------------------------*/
int utr_adapt( /* returns: >0 - success, <=0 - failure */
		int Problem_id,       /* in: problem data structure to be used */
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
)
{

	/* local variables */
	int mesh_id, field_id;
	int	 type;		/* strategy number for adaptation */
	double eps;		/* eps is assumed to be a global tolerance level */
	double el_eps;	/* coefficient for choosing elements to adapt */
	double fam_error;	/* family error */
	double ratio;		/* ratio of errors for derefinements */
	double * error_indi;	/* array with elements errors */
	double sum_error, average_error, max_error;
	int nr_elem;
	int nr_ref, *list_ref;/* number and list of elements to refine */
	int nr_deref, *list_deref;/* number and list of elements to derefine */
	int nmel, nrel; /* number of elements: total, active */
	int nmfa, nrfa; /* number of faces: total, active */
	int nmed, nred; /* number of edges: total, active */
	int nmno, nrno; /* number of nodes: total, active */
	int father, elsons[MMC_MAXELSONS + 1]; /* family information */

	/* auxiliary variables */

	int i, iel, iaux, name, ison;

	int iprint;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* get formulation parameters */
	//i=1; name=pdr_ctrl_i_params(Problem_id,i); changed for string
	i = 2;
	mesh_id = pdr_ctrl_i_params(Problem_id, i);
	i = 3;
	field_id = pdr_ctrl_i_params(Problem_id, i);

	/* get adaptation parameters */
	i = 1;
	type = pdr_adapt_i_params(Problem_id, i);
	i = 5;
	eps = pdr_adapt_d_params(Problem_id, i);
	i = 6;
	ratio = pdr_adapt_d_params(Problem_id, i);
	i = 7;
	iprint = pdr_adapt_i_params(Problem_id, i);

	/* get necessary mesh parameters */
	nrno = mmr_get_nr_node(mesh_id);
	nmno = mmr_get_max_node_id(mesh_id);
	nred = mmr_get_nr_edge(mesh_id);
	nmed = mmr_get_max_edge_id(mesh_id);
	nrfa = mmr_get_nr_face(mesh_id);
	nmfa = mmr_get_max_face_id(mesh_id);
	nrel = mmr_get_nr_elem(mesh_id);
	nmel = mmr_get_max_elem_id(mesh_id);


	if (iprint > 1) {
		//utr_print_adapt_info(field_id,Interactive_output);
	}

	/* uniform refinement */
	if (type == -1) {

		i = utr_adapt_core(Problem_id, Work_dir, &type, 1, REF_ONLY,
						Interactive_input, Interactive_output, iprint);

//    /* initialize refinement data structures */
//    mmr_init_ref(mesh_id);

//    /* perform uniform refinement */
//    i=apr_refine(mesh_id,MMC_DO_UNI_REF);

//    /* project DOFS between generations */
//    iaux=-1;
//    /* project DOFS between generations */
//    apr_proj_dof_ref(field_id, iaux, nmel, nmfa, nmed, nmno);

//    /* restore consistency of mesh data structure and free space */
//    mmr_final_ref(mesh_id);

//    apr_check_field(field_id);

	}/* end uniform refinement */
	else {/* adaptive refinement with error indicators */


		/* allocate space for array with error indicators */
		error_indi = utr_dvector(nmel + 1, "error_indi in adapt");
		nr_ref = 0;
		list_ref = utr_ivector(nmel + 1, "list_ref in adapt");
		nr_deref = 0;
		list_deref = utr_ivector(nmel + 1, "list_deref in adapt");

		/*kbw
		printf("In adapt: type %d, eps %lf, ratio %lf\n",
		type,  eps, ratio);
		//kew*/

		sum_error = 0.0;
		nr_elem = 0;
		max_error = 0.0;
		/* compute error indicators for all active elements */
		iel = 0;
		while ((iel = mmr_get_next_act_elem(mesh_id, iel)) != 0) {

			error_indi[iel] = pdr_err_indi(Problem_id, type, iel);
			sum_error += error_indi[iel];
			if (error_indi[iel] > max_error) max_error = error_indi[iel];
			nr_elem++;

			/*kbw
			printf("in active element %d - error_indi %20.15lf \n",
			iel,error_indi[iel]);
			//kew*/

		}

		average_error = sum_error / nr_elem;
		/*kbw
			printf("In utr_adapt: total error %lf (sqrt %lf), average error %lf\n",
			   sum_error, sqrt(sum_error), average_error);
		//kew*/


		/* here: eps<0.1 denotes global tolerance */
		if (eps < 0.1) el_eps = eps;
		/* set element limit for errors (for equidistribution principle) */
		else if (eps < 1.0) el_eps = eps * max_error;
		else el_eps = eps * average_error;

		/* create lists of elements for refinements and derefinements */
		/* loop over all active elements */
		iel = 0;
		while ((iel = mmr_get_next_act_elem(mesh_id, iel)) != 0) {

			/*kbw
			printf("in active element %d - error_indi %20.15lf <> %20.15lf\n",
			iel,error_indi[iel],el_eps);
			//kew*/

			/* if error small enough and the family not yet considered */
			if (mmr_el_fam(mesh_id, iel, NULL, NULL) > 0 &&
					error_indi[iel] < ratio * el_eps && error_indi[iel] > -1) {

				/* find father and sons */
				father = mmr_el_fam(mesh_id, iel, NULL, NULL);
				mmr_el_fam(mesh_id, father, elsons, NULL);

				/* check whether all sons are active and compute the error for the family*/
				iaux = 0;
				fam_error = 0.0;
				for (ison = 1; ison <= elsons[0]; ison++) {
					if (mmr_el_status(mesh_id, elsons[ison]) <= 0) {
						iaux = 1;
						break;
					}
					else {
						fam_error += error_indi[elsons[ison]];
						/* indicate son should not be considered again */
						if (error_indi[elsons[ison]] < el_eps)
							error_indi[elsons[ison]] = -2;

						/*kbw
						printf("in active element %d, father %d, son %d (%d)\n",
						iel,father,elsons[ison],ison);
						printf("in active element %d - error_indi %20.15lf <> %20.15lf\n",
						elsons[ison],error_indi[elsons[ison]],el_eps);
						//kew*/

					}
				}

				/* if all sons active and family error is smaller than the limit */
				if (iaux == 0 && fam_error < ratio * el_eps) {

					/* if derefinement is not excluded because of irregularity constrained */
// WARNING: we check only iel because in apr_limit_deref we use simplified
//          strategy based on considering the edges of father element
//          hence allowing one son allows all sons to be derefined !!!
					if (apr_limit_deref(field_id, iel) == APC_DEREF_ALLOWED) {

						/*kbw
								printf("element %d (%d) to derefine: fam_error %15.12lf < %15.12lf\n",
								   elsons[1],nr_deref,fam_error,ratio*el_eps);
						//kew*/


						/* put first son (elsons[1]) on list to derefine */
						list_deref[nr_deref] = elsons[1];
						nr_deref++;
					}

				} /* end if all sons active and family error smaller than the limit */

			} /* end if error small in element */
			else if (error_indi[iel] > el_eps) {

				/*kbw
				printf("element %d (%d) to refine: error_indi %15.12lf > %15.12lf\n",
				iel,nr_ref,error_indi[iel],el_eps);
				//kew*/

				if (apr_limit_ref(field_id, iel) == APC_REF_ALLOWED) {
					/* put iel on list to refine */
					list_ref[nr_ref] = iel;
					nr_ref++;

				}

			}

		} /* end loop over active elements */

		/*kbw
			printf("List to derefine (%d elements): \n",nr_deref);
			for(iel=0;iel<nr_deref;iel++) printf("%d ",list_deref[iel]);
			printf("\n");
			//getchar();
		//kew*/

		/* DEREFINEMENTS */


		/*kbw
			printf("List to derefine after exchange (%d elements): \n",nr_deref);
			for(iel=0;iel<nr_deref;iel++) printf("%d ",list_deref[iel]);
			printf("\n");
			//getchar();
		//kew*/

		if (nr_deref > 0) {

			iaux = utr_adapt_core(Problem_id, Work_dir,
							list_deref, nr_deref, DEREF_ONLY,
							Interactive_input, Interactive_output, iprint);

		}
//    /* initialize refinement data structures */
//    mmr_init_ref(mesh_id);

//    /* perform derefinements */
//    for(iel=0;iel<nr_deref;iel++) {
//      if(mmr_el_status(mesh_id, list_deref[iel])==MMC_ACTIVE){
//	if(apr_limit_deref(field_id, list_deref[iel])==APC_DEREF_ALLOWED){
//	  iaux=apr_derefine(mesh_id,list_deref[iel]);
//	  if(iaux<0) {
//	    printf("Unsuccessful derefinement of El %d. Exiting!\n",
//		   list_deref[iel]);
//	    exit(-1);
//	  }
//	}
//      }
//    }


//    /* project DOFS between generations */
//    iaux=-1;
//    /* project DOFS between generations */
//    apr_proj_dof_ref(field_id, iaux, nmel, nmfa, nmed, nmno);

//    /* restore consistency of mesh data structure and free space */
//    mmr_final_ref(mesh_id);


		/* REFINEMENTS */

		/*kbw
			printf("List to refine before update (%d elements): \n",nr_ref);
			for(iel=0;iel<nr_ref;iel++) printf("%d ",list_ref[iel]);
			printf("\n");
			getchar();
		//kew*/


		/*kbw
			  printf("List to refine after exchange (%d elements): \n",nr_ref);
			  for(iel=0;iel<nr_ref;iel++) printf("%d ",list_ref[iel]);
			  printf("\n");
			  //getchar();
		//kew*/

		/* update list of refined elements due to irregularity constraint
		mmr_update_ref_list(mesh_id, &nr_ref, list_ref);
		*/

		/*kbw
			printf("List to refine after update (%d elements): \n",nr_ref);
			for(iel=0;iel<nr_ref;iel++) printf("%d ",list_ref[iel]);
			printf("\n");
			//getchar();
		//kew*/

		if (nr_ref > 0) {

			iaux = utr_adapt_core(Problem_id, Work_dir,
							list_ref, nr_ref, REF_ONLY,
							Interactive_input, Interactive_output, iprint);

		}

//    /* initialize refinement data structures */
//    mmr_init_ref(mesh_id);

//    /* next perform refinements */
//    for(iel=0;iel<nr_ref;iel++) {
//      if(mmr_el_status(mesh_id, list_ref[iel])==MMC_ACTIVE){
//	if(apr_limit_ref(field_id, list_ref[iel])==APC_REF_ALLOWED){
//	  iaux=apr_refine(mesh_id,list_ref[iel]);
//	  if(iaux<0) {
//	    printf("Unsuccessful derefinement of El %d. Exiting!\n",
//           list_ref[iel]);
//	    exit(-1);
//	  }
//	}
//      }
//    }


//    /* project DOFS between generations */
//    iaux=-1;
//    apr_proj_dof_ref(field_id, iaux, nmel, nmfa, nmed, nmno);

//    /* restore consistency of mesh data structure and free space */
//    mmr_final_ref(mesh_id);

//    apr_check_field(field_id);


		/* free the space */
		free(error_indi);
		free(list_ref);
		free(list_deref);


	} /* end if not uniform refinements */

	/* get necessary mesh parameters */

	//if(iprint>1){
	//utr_print_adapt_info(field_id,Interactive_output);
	//}


	return (1);

}

/*---------------------------------------------------------
utr_manual_refinement - to perform manual refinement/unrefinement
---------------------------------------------------------*/
int utr_manual_refinement(  /* returns: >0 - success, <=0 - failure */
		int Problem_id,       /* in: problem data structure to be used */
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
)
{

	/* local variables */
	int i, iaux = 0; //, jaux, iprob, icount
//  int lel[100], nreq, nr_sol, iel, nr_mat;
	int field_id, mesh_id;
	int policy;

#ifdef PARALLEL
	if (pcr_my_proc_id() == pcr_print_master()) {
#endif
		if (Interactive_input == stdin) {

			/* refine or derefine the mesh */
			printf("set element number ( >0 - refine, <0 - derefine, 0 - stop\n");
			printf("\t\t     -1 - uniform refinement, -2 - uniform derefinement):");
			scanf("%d", &iaux);

		} else {

			// read control data from Interactive_input file
			fscanf(Interactive_input, "%d\n", &iaux);

		}
#ifdef PARALLEL
	}

	pcr_bcast_int(pcr_print_master(), 1, &iaux);
	//printf("After BCAST %d\n",iaux);

#endif


	if (iaux == 0) {
		return 1;
	}
	/* initialize refinement data structures */
	// mmr_init_ref(mesh_id);

//  if(iaux==-1) i=apr_refine(mesh_id,MMC_DO_UNI_REF);
//  else if(iaux > 0 && iaux<1e10) {
//    i=apr_refine(mesh_id,iaux);
//#ifdef DEBUG
//    if(i<0) printf("Unsuccessful refinement of element %d\n",iaux);
//#endif
//  }
//  else if(iaux == -2) i=apr_derefine(mesh_id,MMC_DO_UNI_DEREF);
//  else if(iaux < -2 && iaux>-1e10) i=apr_derefine(mesh_id,-iaux);
//  else iaux=0;

//  /* project DOFS between generations */
//  if(iaux<-2) iaux=-iaux;
//  apr_proj_dof_ref(field_id, iaux, nmel, nmfa, nmed, nmno);


	switch (iaux) {
		case -1:
			policy = REF_ONLY;
			iaux = MMC_DO_UNI_REF;
			break;
		case -2:
			policy = DEREF_ONLY;
			iaux = MMC_DO_UNI_DEREF;
			break;
		default:
			policy = (iaux > 0) ? REF_ONLY
					:    DEREF_ONLY;
			iaux = (iaux > 0) ? iaux : -iaux;
			break;
	}

	i = utr_adapt_core(Problem_id, Work_dir, &iaux, 1, policy,
					Interactive_input, Interactive_output, 1);

	return (i);
}

/*---------------------------------------------------------
utr_test_refinements - to perform sequence of random refinements/unrefinements
---------------------------------------------------------*/
int utr_test_refinements(  /* returns: >0 - success, <=0 - failure */
		int Problem_id,       /* in: problem data structure to be used */
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
)
{

	/* local variables */
	int iaux, jaux, icount; //i, , iprob
	int iel;
	int field_id, mesh_id;
	double daux;//, t_wall;
	int nrno, nmno, nred, nmed, nrfa, nmfa, nrel, nmel;
	int success;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	field_id = Problem_id;
	mesh_id = apr_get_mesh_id(field_id);

	if (Interactive_input == stdin) {

		printf("set number of random refinements: ");
		scanf("%d", &icount);

	} else {

		// read control data from interactive_input file
		fscanf(Interactive_input, "%d\n", &icount);

	}

	nrno = mmr_get_nr_node(mesh_id);
	nmno = mmr_get_max_node_id(mesh_id);
	nred = mmr_get_nr_edge(mesh_id);
	nmed = mmr_get_max_edge_id(mesh_id);
	nrfa = mmr_get_nr_face(mesh_id);
	nmfa = mmr_get_max_face_id(mesh_id);
	nrel = mmr_get_nr_elem(mesh_id);
	nmel = mmr_get_max_elem_id(mesh_id);


	fprintf(Interactive_output, "\nBefore refinement test.\n");
	utr_print_adapt_info(field_id, Interactive_output);


	for (jaux = 0; jaux < abs(icount); jaux++) {

		/* int nr_patches, ino; */
		/* int nr_sol = 0; // consider first vector of unknowns */
		/* nr_patches = utr_create_patches_small(field_id,  &pdv_patches); */
		/* utr_recover_derivatives_small(field_id, nr_sol, nr_patches, pdv_patches); */
		/* for(ino=1; ino<=nr_patches; ino++){  */
		/*   if(pdv_patches[ino].deriv!=NULL) free(pdv_patches[ino].deriv); */
		/* } */
		/* free(pdv_patches); */

		fprintf(Interactive_output, "\n");

		nrno = mmr_get_nr_node(mesh_id);
		nmno = mmr_get_max_node_id(mesh_id);
		nred = mmr_get_nr_edge(mesh_id);
		nmed = mmr_get_max_edge_id(mesh_id);
		nrfa = mmr_get_nr_face(mesh_id);
		nmfa = mmr_get_max_face_id(mesh_id);
		nrel = mmr_get_nr_elem(mesh_id);
		nmel = mmr_get_max_elem_id(mesh_id);


//    mmr_init_ref(mesh_id);

		success = -1;
		while (success < 0) {

			/*  random number from 1 to nmel: */
			daux = (double)rand() / (double)RAND_MAX; //drand48();
			//daux=drand48();
			iaux = (int) (nmel * daux);
			if (iaux < 1) iaux = 1;
			if (iaux > nmel) iaux = nmel;

			success = utr_adapt_core(Problem_id, Work_dir, &iaux, 1, REF_ONLY,
							Interactive_input, Interactive_output, 0);

//      if( mmr_el_status(mesh_id, iaux)==MMC_ACTIVE &&
//      mmr_el_gen(mesh_id, iaux) < mmr_get_max_gen(mesh_id) ){
//	success = apr_refine(mesh_id,iaux);

			if (success < 0)
				fprintf(Interactive_output,
						"Trial %d: nmel %d, element to refine %d - failure\n",
						jaux, nmel, iaux);
			else
				fprintf(Interactive_output,
						"Trial %d: nmel %d, element to refine %d - success\n",
						jaux, nmel, iaux);

		}

//    }

//    /* update dofs structure */
//    iaux = -1; // indicate rewriting dofs after massive derefinements
//    apr_proj_dof_ref(field_id, iaux, nmel, nmfa, nmed, nmno);

//    /* clean space */
//    mmr_final_ref(mesh_id);
//    apr_check_field(field_id);

		if (icount > 0) {
			int attempt = 0;

			//   mmr_init_ref(mesh_id);

			success = -1;
			while (success < 0 && attempt < 100) {

				/*  random number from 1 to nmel: */
				daux = (double)rand() / (double)RAND_MAX; //drand48();
				//daux=drand48();
				iaux = (int) (nmel * daux);
				if (iaux < 3) iaux = 3;
				if (iaux > nmel) iaux = nmel;
				iel = iaux;

				attempt++;
				if (mmr_el_status(mesh_id, iel) == MMC_ACTIVE) {

					int ison;//, ifa, face, son, son_faces[MMC_MAXELFAC+1];
					int father, elsons[MMC_MAXELSONS + 1]; /* family information */
					/* find father and sons */
					//father=mmr_el_fam(mesh_id,iel,elsons,NULL);
					father = mmr_el_fam(mesh_id, iel, NULL, NULL);

					if (father > 0) {
						// Moved here by KM 08.2011
						// If father is 0, it throws an error!
						// Instead of moving the call here, mmr_el_fam should be modified
						// to return elsons[0]=0 - KB 08.2011
						mmr_el_fam(mesh_id, father, elsons, NULL);

						/* check whether all sons are active and compute the error for the family*/
						iaux = 0;
						for (ison = 1; ison <= elsons[0]; ison++) {
							/*kbw
							printf("in active element %d, father %d, son %d (%d)\n",
							iel,father,elsons[ison],ison);
							//kew*/
							if (mmr_el_status(mesh_id, elsons[ison]) <= 0) {
								iaux = 1;
								break;
							}


						}
					} else {
						iaux = 1;
					}

					/* if all sons active */
					if (iaux == 0) {

						/* if derefinement is not excluded because of irregularity constrained */
						if (apr_limit_deref(field_id, iel) == APC_DEREF_ALLOWED) {

							/*kbw
							printf("element %d (%d) to derefine: fam_error %15.12lf < %15.12lf\n",
							elsons[1],nr_deref,fam_error,ratio*el_eps);
							//kew*/


							/* put first son (elsons[1]) on list to derefine */

//	      success = apr_derefine(mesh_id,iel);
							success = utr_adapt_core(Problem_id, Work_dir, &iel, 1, DEREF_ONLY,
											Interactive_input, Interactive_output, 0);

							if (success < 0) {
								fprintf(Interactive_output,
										"Trial %d, attempt %d: nmel %d, element to derefine %d - failure\n",
										jaux, attempt, nmel, iel);
								exit(-1);
							} else {
								fprintf(Interactive_output,
										"Trial %d, attempt %d: nmel %d, element to derefine %d - success\n",
										jaux, attempt, nmel, iel);
							}

//	      /* update dofs structure */
//	      iaux = iel;
//	      // indicate rewriting dofs after derefinement of one element
//	      apr_proj_dof_ref(field_id, iaux, nmel, nmfa, nmed, nmno);

//	      /* clean space */
//	      mmr_final_ref(mesh_id);
//	      apr_check_field(field_id);

						}

					} /* end if all sons active and family error smaller than the limit*/
				}	 /* end if element active */
			} /* until succeded */
		}  /* end if derefinements requested */
	} /* end loop over trials */
	fprintf(Interactive_output, "\nAfter refinement test.\n");
	utr_print_adapt_info(field_id, Interactive_output);

	return (1);

}


