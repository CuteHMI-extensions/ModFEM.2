/**********************************************************************
File sis_lapack_intf.c - implementation of the linear solver interface
							for the dense LAPACK solver

Contains:

Interface routines:
  sir_direct_module_introduce - to return the solver name
  sir_direct_solve_lin_sys - to perform the five steps below in one call
  sir_direct_init - to create new solver instance and read its control parameters
  sir_direct_create - to create and initialize solver data structure
  sir_direct_solve - to solve the system for a given data
  sir_direct_free - to free memory for a given solver instance
  sir_direct_destroy - to destroy a solver instance

Local procedures:
  sir_direct_assemble_stiff_mat - to assemble entries to the global stiffness
						   matrixand the global load vector using the provided
						   stiffness matrix and load vector local
------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>

/* problem dependent interface between approximation and solver modules */
#include <modfem/pdh_intf.h>

/* interface for general purpose linear solver modules */
#include <modfem/sih_intf.h>

/* internal information for the LAPACK linear solver module */
#include <modfem/si_lapack/sih_lapack.h>

/* LAPACK procedures */
#include <modfem/lin_alg_intf.h>

#include <modfem/uth_log.h>

/*#define DEBUG_SIM 1 */

/*** GLOBAL VARIABLES (for the solver module) ***/
int   siv_direct_nr_solvers = 0;    /* the number of solvers in the problem */
int   siv_direct_cur_solver_id;                /* ID of the current problem */
sit_direct_solver siv_direct_solvers[SIC_MAX_NUM_SOLV]; /* array of solvers */


/* utility procedure */
int sir_direct_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
		/*  <0 - position at which put on the list */
		/*   0 - list full, not found on the list */
		int Num, 	/* in: number to put on the list */
		int * List, 	/* in: list */
		int Ll		/* in: total list's lengths */
);

/*------------------------------------------------------------
  sir_direct_module_introduce - to return the solver name
------------------------------------------------------------*/
int sir_direct_module_introduce( /* returns: >=0 - success code, <0 - error code */
		char * Solver_name /* out: the name of the solver */
)
{
	char * string = "LAPACK_DENSE";

	strcpy(Solver_name, string);

	return (1);
}

/*------------------------------------------------------------
  sir_direct_solve_lin_sys - to solve the system of linear equations for the current
		  setting of parameters (with initialization and finalization phases)
------------------------------------------------------------*/
int sir_direct_solve_lin_sys( /* returns: >=0 - success code, <0 - error code */
		int Problem_id,    /* ID of the problem associated with the solver */
		int Parallel,      /* parameter specifying sequential (SIC_SEQUENTIAL) */
		/* or parallel (SIC_PARALLEL) execution */
		char * Filename /* in: name of the file with control parameters */
)
{


	int comp_type, solver_id, monitor, ini_guess;
	//char solver_name[100];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/*kbw
	printf("In sir_direct_solve_lin_sys before init: problem_id %d\n",
	 problem_id);
	/*kew*/

	/* initialize the solver - empty data, nothing is read by LAPACK module */
	solver_id = sir_direct_init(Parallel, Filename);

	/*kbw
	printf("In sir_direct_solve_lin_sys before create: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
	/*kew*/

	/* create the solver data structure and asociate it with a given problem */
	sir_direct_create(solver_id, Problem_id);

	/*kbw
	printf("In sir_direct_solve_lin_sys before solve: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
	/*kew*/

	/* launch the solver */
	comp_type = SIC_SOLVE;
	monitor = SIC_PRINT_INFO;
	ini_guess = 0;
	sir_direct_solve(solver_id, comp_type, monitor);

	/*kbw
	printf("In sir_direct_solve_lin_sys before free: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
	/*kew*/

	/* free the solver data structure - together with renumbering data */
	sir_direct_free(solver_id);

	/*kbw
	printf("In sir_direct_solve_lin_sys before destroy: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
	/*kew*/

	/* destroy the solver */
	sir_direct_destroy(solver_id);

	return (0);
}


/*------------------------------------------------------------
  sir_direct_init - to create a new solver and read its control parameters
------------------------------------------------------------*/
int sir_direct_init( /* returns: >0 - solver ID, <0 - error code */
		int Parallel,      /* parameter specifying sequential (SIC_SEQUENTIAL) */
		/* or parallel (SIC_PARALLEL) execution */
		char * Filename /* in: name of the file with control parameters */
)
{
	/* increase the counter for solvers */
	siv_direct_nr_solvers++;

	/* set the current solver ID */
	siv_direct_cur_solver_id = siv_direct_nr_solvers;

	return (siv_direct_cur_solver_id);
}


/*------------------------------------------------------------
  sir_direct_create - to create solver's data structure
------------------------------------------------------------*/
int sir_direct_create( /* returns: >0 - solver ID, <0 - error code */
		int Solver_id,    /* in: solver identification */
		int Problem_id    /* ID of the problem associated with the solver */
)
{

	/* pointer to solver structure */
	sit_direct_solver * solver_p;
	/* pointer to dofs structure */
	sit_direct_dof_struct * dof_struct_p = NULL;

	/* the number of (different) mesh entities for which entries to the global */
	/* stiffness matrix and load vector will be provided */
	int nr_int_ent, nr_dof_ent, max_dofs_per_dof_ent, max_dofs_int_ent;
	/* the global number of degrees of freedom */
	int nrdofs_glob, int_ent_id, nr_dof_ent_loc, idofent, dof_ent_id = 0, pos_glob;
	int nr_dof_struct, dof_struct_id, dof_ent_type;
	int * temp_list_dof_type;
	int * temp_list_dof_id;
	int * temp_list_dof_nrdofs;
	int l_dof_ent_types[SIC_MAX_DOF_PER_INT], l_dof_ent_ids[SIC_MAX_DOF_PER_INT];
	int l_dof_ent_nrdofs[SIC_MAX_DOF_PER_INT];

	/* auxiliary variables */
	int i, idof, jdof, iaux, ient, ineig;

	/*++++++++++++++++ executable statements ++++++++++++++++*/


	/* get pointer to solver structure */
	solver_p = sir_direct_select_solver(Solver_id);
	solver_p->problem_id = Problem_id;

	/*kbw
	printf("In sir_direct_create: solver_id %d, problem_id %d\n",
	 Solver_id, solver_p->problem_id);
	/*kew*/

	/* prepare renumbering arrays */
	/*------------------------------------------------------------
	  pdr_get_list_int_ent - to return to the solver module :
							  1. the list of integration entities - entities
								 for which stiffness matrices and load vectors are
								 provided by the FEM code
							  2. the list of DOF entities - entities with which
								 there are dofs associated by the given approximation
	------------------------------------------------------------*/
	pdr_get_list_ent(solver_p->problem_id, &nr_int_ent,
			&solver_p->l_int_ent_type, &solver_p->l_int_ent_id,
			&nr_dof_ent, &temp_list_dof_type,
			&temp_list_dof_id, &temp_list_dof_nrdofs,
			&nrdofs_glob, &max_dofs_per_dof_ent);


	/*kbw
		printf("In sir_direct_create after pdr_get_list_ent\n");
		printf("nr_int_ent %d\n", nr_int_ent);
		for(i=0;i<nr_int_ent;i++)  printf("type %d, id %d\n",
			solver_p->l_int_ent_type[i],solver_p->l_int_ent_id[i]);
		printf("\nNr_dof_ent %d, Nrdofs_glob %d, Max_dofs_per_ent %d\n",
		   nr_dof_ent, nrdofs_glob, max_dofs_per_dof_ent);
		for(i=0;i<nr_dof_ent;i++)  printf("type %d, id %d, nrdofs %d\n",
			temp_list_dof_type[i], temp_list_dof_id[i],
						  temp_list_dof_nrdofs[i]);
	/*kew*/

	solver_p->nr_int_ent = nr_int_ent;
	solver_p->nr_dof_ent = nr_dof_ent;
	solver_p->nrdofs_glob = nrdofs_glob;

	/* array of structures storing DOF data */
	solver_p->l_dof_struct =
			(sit_direct_dof_struct *)malloc( nr_dof_ent * sizeof(sit_direct_dof_struct) );

	/* dof_ent_index to dof_struct_index (based on which dof_ent_id and */
	/* dof_ent_type can be find */
	solver_p->max_dof_vert_id = -1;
	solver_p->max_dof_edge_id = -1;
	solver_p->max_dof_face_id = -1;
	solver_p->max_dof_elem_id = -1;
	for (i = 0; i < nr_dof_ent; i++) {

		if (temp_list_dof_type[i] == PDC_ELEMENT ) {
			if ( temp_list_dof_id[i] > solver_p->max_dof_elem_id )
				solver_p->max_dof_elem_id = temp_list_dof_id[i];
		} else if (temp_list_dof_type[i] == PDC_FACE ) {
			if ( temp_list_dof_id[i] > solver_p->max_dof_face_id )
				solver_p->max_dof_face_id = temp_list_dof_id[i];
		} else if (temp_list_dof_type[i] == PDC_EDGE ) {
			if ( temp_list_dof_id[i] > solver_p->max_dof_edge_id )
				solver_p->max_dof_edge_id = temp_list_dof_id[i];
		} else if (temp_list_dof_type[i] == PDC_VERTEX ) {
			if ( temp_list_dof_id[i] > solver_p->max_dof_vert_id )
				solver_p->max_dof_vert_id = temp_list_dof_id[i];
		} else {
			printf("Error 87373732 in sis_lapack_intf/sir_direct_create!!! Exiting\n");
			exit(-1);
		}
	}

	if (solver_p->max_dof_elem_id >= 0) {
		solver_p->l_dof_elem_to_struct =
				(int *)malloc((solver_p->max_dof_elem_id + 1) * sizeof(int));
		for (i = 0; i <= solver_p->max_dof_elem_id; i++)
			solver_p->l_dof_elem_to_struct[i] = -1;
	}
	if (solver_p->max_dof_face_id >= 0) {
		solver_p->l_dof_face_to_struct =
				(int *)malloc((solver_p->max_dof_face_id + 1) * sizeof(int));
		for (i = 0; i <= solver_p->max_dof_face_id; i++)
			solver_p->l_dof_face_to_struct[i] = -1;
	}
	if (solver_p->max_dof_edge_id >= 0) {
		solver_p->l_dof_edge_to_struct =
				(int *)malloc((solver_p->max_dof_edge_id + 1) * sizeof(int));
		for (i = 0; i <= solver_p->max_dof_edge_id; i++)
			solver_p->l_dof_edge_to_struct[i] = -1;
	}
	if (solver_p->max_dof_vert_id >= 0) {
		solver_p->l_dof_vert_to_struct =
				(int *)malloc((solver_p->max_dof_vert_id + 1) * sizeof(int));
		for (i = 0; i <= solver_p->max_dof_vert_id; i++)
			solver_p->l_dof_vert_to_struct[i] = -1;
	}


	/* in case of no renumbering global position of dof block is specified */
	pos_glob = 0;
	for (idof = 0; idof < nr_dof_ent; idof++) {

		dof_struct_p = &solver_p->l_dof_struct[idof];

		nr_dof_ent_loc = temp_list_dof_nrdofs[idof];
		dof_struct_p->dof_ent_type = temp_list_dof_type[idof];
		dof_struct_p->dof_ent_id = temp_list_dof_id[idof];
		dof_struct_p->nrdofs = nr_dof_ent_loc;
		dof_struct_p->posglob = pos_glob;
		pos_glob += nr_dof_ent_loc;

		/* initialize lists of integration entities and neighbouring dof_ent */
		dof_struct_p->nr_int_ent = 0;
		for (i = 0; i < SIC_MAX_INT_PER_DOF; i++)
			dof_struct_p->l_int_ent_index[i] = SIC_LIST_END_MARK;
		dof_struct_p->nrneig = 0;
		for (i = 0; i < SIC_MAX_DOF_STR_NGB; i++)
			dof_struct_p->l_neig[i] = SIC_LIST_END_MARK;

		/* arrays dof_ent to dof_struct */
		if (dof_struct_p->dof_ent_type == PDC_ELEMENT ) {
			solver_p->l_dof_elem_to_struct[dof_struct_p->dof_ent_id] = idof;

#ifdef DEBUG_SIM
			if (dof_struct_p->dof_ent_id > solver_p->max_dof_elem_id) {
				printf("Error 84543732 in sis_lapack_intf/sir_direct_create!!! Exiting\n");
				printf("%d > %d\n", dof_struct_p->dof_ent_id, solver_p->max_dof_elem_id);
				exit(-1);
			}
#endif

		} else if (dof_struct_p->dof_ent_type == PDC_FACE ) {
			solver_p->l_dof_face_to_struct[dof_struct_p->dof_ent_id] = idof;

#ifdef DEBUG_SIM
			if (dof_struct_p->dof_ent_id > solver_p->max_dof_face_id) {
				printf("Error 84543732 in sis_lapack_intf/sir_direct_create!!! Exiting\n");
				printf("%d > %d\n", dof_struct_p->dof_ent_id, solver_p->max_dof_face_id);
				exit(-1);
			}
#endif

		} else if (dof_struct_p->dof_ent_type == PDC_EDGE ) {
			solver_p->l_dof_edge_to_struct[dof_struct_p->dof_ent_id] = idof;

#ifdef DEBUG_SIM
			if (dof_struct_p->dof_ent_id > solver_p->max_dof_edge_id) {
				printf("Error 84543732 in sis_lapack_intf/sir_direct_create!!! Exiting\n");
				printf("%d > %d\n", dof_struct_p->dof_ent_id, solver_p->max_dof_edge_id);
				exit(-1);
			}
#endif

		} else if (dof_struct_p->dof_ent_type == PDC_VERTEX ) {
			solver_p->l_dof_vert_to_struct[dof_struct_p->dof_ent_id] = idof;
#ifdef DEBUG_SIM
			if (dof_struct_p->dof_ent_id > solver_p->max_dof_vert_id) {
				printf("Error 84543732 in sis_lapack_intf/sir_direct_create!!! Exiting\n");
				printf("%d > %d\n", dof_struct_p->dof_ent_id, solver_p->max_dof_vert_id);
				exit(-1);
			}
#endif
		} else {
			printf("Error 8963732 in sis_lapack_intf/sir_direct_create!!! Exiting\n");
			exit(-1);
		}




		/*kbw
			printf("In sir_direct_create after filling dof_struct no %d\n", idof);
			printf("dof_ent_type %d, dof_ent_id %d, nrdofs %d, posglob %d\n",
			   dof_struct_p->dof_ent_type , dof_struct_p->dof_ent_id,
			   dof_struct_p->nrdofs, dof_struct_p->posglob);
			printf("Initialized lists of int_ent %d, neig %d\n",
			   dof_struct_p->l_int_ent_index[0], dof_struct_p->l_neig[0]);
		/*kew*/

	}

#ifdef DEBUG_SIM
	if ( solver_p->nrdofs_glob != pos_glob ) {
		printf("Error 843732 in sis_lapack_intf/sir_direct_create!!! Exiting\n");
		exit(-1);
	}
#endif

	/* getting information on the structure of the global stiffness matrix */
	nr_dof_struct = 0;
	max_dofs_int_ent = 0;
	for (ient = 0; ient < nr_int_ent; ient++) {

		int nrdofs_int_ent = 0;
		int_ent_id = solver_p->l_int_ent_id[ient];
		nr_dof_ent_loc = SIC_MAX_DOF_PER_INT;

		/* first call to pdr_comp_stiff_mat to get information on DOF entities */
		/*   associated with the given integration entity (types, IDs, nrdofs)
		/*   - stiffness matrix and load vector are NOT computed only the */
		/*   structures of the local and the global stiffness matrices are deduced */
		pdr_comp_stiff_mat(solver_p->problem_id, solver_p->l_int_ent_type[ient],
				solver_p->l_int_ent_id[ient], PDC_NO_COMP, NULL,
				&nr_dof_ent_loc, l_dof_ent_types,
				l_dof_ent_ids, l_dof_ent_nrdofs,
				NULL, NULL, NULL, NULL);


		/*kbw
		if(ient==3979){
		   printf(
		"in sir_direct_create after pdr_comp_stiff_mat for int_ent no %d (id %d), type %d:\n",
		ient, int_ent_id, solver_p->l_int_ent_type[ient]);
			printf("nr_dof_ent_loc %d, types, ids, nrdofs:\n", nr_dof_ent_loc);
			for(idofent=0; idofent<nr_dof_ent_loc; idofent++){
			  printf("%d %d %d\n", l_dof_ent_types[idofent],
				 l_dof_ent_ids[idofent], l_dof_ent_nrdofs[idofent]);
			}
			}
		/*kew*/

		for (idofent = 0; idofent < nr_dof_ent_loc; idofent++) {

			dof_ent_id = l_dof_ent_ids[idofent];
			dof_ent_type = l_dof_ent_types[idofent];

			if (dof_ent_type == PDC_ELEMENT ) {

				dof_struct_id = solver_p->l_dof_elem_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				assert(dof_struct_id > -1);
				if ( solver_p->l_dof_elem_to_struct[dof_ent_id] == -1) {
					printf("Error 3472941 in sir_direct_create!!! Exiting\n");
					exit(-1);
				}
#endif

			} else if (dof_ent_type == PDC_FACE ) {

				dof_struct_id = solver_p->l_dof_face_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				assert(dof_struct_id > -1);
				if ( solver_p->l_dof_face_to_struct[dof_ent_id] == -1) {
					printf("Error 3472942 in sir_direct_create!!! Exiting\n");
					exit(-1);
				}
#endif

			} else if (dof_ent_type == PDC_EDGE ) {

				dof_struct_id = solver_p->l_dof_edge_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				assert(dof_struct_id > -1);
				if ( solver_p->l_dof_edge_to_struct[dof_ent_id] == -1) {
					printf("Error 3472943 in sir_direct_create!!! Exiting\n");
					exit(-1);
				}
#endif

			} else if (dof_ent_type == PDC_VERTEX ) {

				dof_struct_id = solver_p->l_dof_vert_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				assert(dof_struct_id > -1);
				if ( solver_p->l_dof_vert_to_struct[dof_ent_id] == -1) {
					printf("in sir_direct_create after pdr_comp_stiff_mat for int_ent no %d (id %d):\n",
							ient, int_ent_id);
					printf("Error 3472944 in sir_direct_create for dof %d!!! Exiting\n",
							dof_ent_id);
					exit(-1);
				}
#endif

			} else {
				printf("Error 34fsf7294 in sir_direct_create!!! Exiting\n");
				exit(-1);
			}

			assert(dof_struct_id > -1);
			dof_struct_p = &solver_p->l_dof_struct[dof_struct_id];

			/*kbw
				printf("for int_type %d, int_id %d, dof_type %d, dof_id %d, struct %d\n",
					   solver_p->l_int_ent_type[ient], solver_p->l_int_ent_id[ient],
					   l_dof_ent_types[idofent], l_dof_ent_ids[idofent],dof_struct_id );
			/*kew*/

#ifdef DEBUG_SIM

			if ((dof_struct_p->dof_ent_id != dof_ent_id) ||
					(dof_struct_p->dof_ent_type != dof_ent_type) ||
					(dof_struct_p->nrdofs != l_dof_ent_nrdofs[idofent]) ) {
				printf("Error 3827 in sir_direct_create!!! Exiting");
				exit(-1);
			}
#endif

			nrdofs_int_ent += l_dof_ent_nrdofs[idofent];

			/*kbw
				  printf("putting int_ent no %d on the list of int_ent, nr_int_ent %d\n",
					 ient, nr_int_ent);
				  printf("before:");
				  for(i=0;i<SIC_MAX_INT_PER_DOF;i++) {
				printf("%d",dof_struct_p->l_int_ent_index[i]) ;
				  }
				  printf("\n");
			/*kew*/

			/*       iaux=sir_direct_put_list(ient,  */
			/* 			dof_struct_p->l_int_ent_index, SIC_MAX_INT_PER_DOF); */
			/*       if(iaux<0) dof_struct_p->nr_int_ent++; */

			/* #ifdef DEBUG_SIM */
			/* 	  if(iaux == 0){ // list full - increase SIC_MAX_INT_PER_DOF */
			/* 	    printf("Error 383627 in sir_direct_create!!! Exiting"); */
			/* 	    exit(-1); */
			/* 	  } */
			/* #endif       */

			/*kbw
				  printf("putting int_ent no %d on the list of int_ent, nr_int_ent %d\n",
					 ient, nr_int_ent);
				  printf("after:");
				  for(i=0;i<SIC_MAX_INT_PER_DOF;i++) {
				printf("%d",dof_struct_p->l_int_ent_index[i]) ;
				  }
				  printf("\n");
			/*kew*/

			for (ineig = 0; ineig < nr_dof_ent_loc; ineig++) {

				if (ineig != idofent) {

					int neig_id = l_dof_ent_ids[ineig];
					int neig_type = l_dof_ent_types[ineig];
					int neig_index = 0;
					if (neig_type == PDC_ELEMENT) {
						neig_index = solver_p->l_dof_elem_to_struct[neig_id];
					} else if (neig_type == PDC_FACE) {
						neig_index = solver_p->l_dof_face_to_struct[neig_id];
					} else if (neig_type == PDC_EDGE) {
						neig_index = solver_p->l_dof_edge_to_struct[neig_id];
					} else if (neig_type == PDC_VERTEX) {
						neig_index = solver_p->l_dof_vert_to_struct[neig_id];
					}

					/*kbw
						  printf("dof_ent %d: putting ineig no %d (id %d) on the list of neig, nrneig %d\n",
							 idofent, ineig, neig_id, dof_struct_p->nrneig);
						  printf("before:");
						  for(i=0;i<SIC_MAX_DOF_STR_NGB;i++) {
						printf("%d",dof_struct_p->l_neig[i]) ;
						  }
						  printf("\n");
					/*kew*/

					iaux = sir_direct_put_list(neig_index,
									dof_struct_p->l_neig, SIC_MAX_DOF_STR_NGB);
					if (iaux < 0) {

						dof_struct_p->nrneig++;

					}

#ifdef DEBUG_SIM
					if (iaux == 0) { // list full - increase SIC_MAX_DOF_STR_NGB
						printf("Error 385627 in sir_direct_create!!! Exiting");
						exit(-1);
					}
#endif

					/*kbw
						  printf("putting ineig no %d (id %d) on the list of neig, nrneig %d\n",
							 ineig, neig_id, dof_struct_p->nrneig);
						  printf("after:");
						  for(i=0;i<SIC_MAX_DOF_STR_NGB;i++) {
						printf("%d",dof_struct_p->l_neig[i]) ;
						  }
						  printf("\n");
					/*kew*/

				}

			}

		} /* end loop over dof_ents of a given int_ent */

		if (nrdofs_int_ent > max_dofs_int_ent) max_dofs_int_ent = nrdofs_int_ent;

	} /* end loop over int_ent */

	solver_p->max_dofs_int_ent = max_dofs_int_ent;

	// PLACE FOR CALL TO RENUMBERING PROCEDURE

#ifdef DEBUG_SIM
	printf("\nLAPACK solver module run in DEBUG mode!\n\n");
	printf("\nEntering LAPACK dense solver for linear equations\n");
	printf("Problem size %d degrees of freedom, storage %lf MBytes\n"
			, nrdofs_glob, (nrdofs_glob + 1)*nrdofs_glob * sizeof(double) / (1024.0 * 1024.0));
	if (nrdofs_glob > 10000) {
		//printf("Problem size too big. Exiting!\n"); exit(0);
		printf("Type any key to continue or CTRL C to quit\n");
		getchar();
		getchar();
	}
#endif
	if (nrdofs_glob > 5000) {
		printf("\nEntering LAPACK dense solver for linear equations\n");
		printf("Problem size %d degrees of freedom, storage %lf MBytes\n"
				, nrdofs_glob, (nrdofs_glob + 1)*nrdofs_glob * sizeof(double) / (1024.0 * 1024.0));
//    printf("Type any key to continue or CTRL C to quit\n"); getchar();getchar();
		//printf("Problem size [no dofs=%d] too big. Exiting!\n",nrdofs_glob); exit(0);
	}

#ifdef DEBUG_SIM
	/*kbw
	  printf("Before allocating global SM and RHS: nrdofs_glob %d\n",
	   nrdofs_glob);
	  printf("solver_p->stiff_mat %lu, solver_p->rhs_vect %lu, solver_p->sol_vect %lu, aux_vect %lu\n",
	   solver_p->stiff_mat,solver_p->rhs_vect,solver_p->sol_vect, solver_p->aux_vect);
	/*kew*/
#endif

	solver_p->stiff_mat = (double *)malloc(nrdofs_glob * nrdofs_glob * sizeof(double));
	mf_check(solver_p->stiff_mat != NULL, "Out of memory! Requested %lf MB (dofs %d)", (double)nrdofs_glob * (double)nrdofs_glob * (double) sizeof(double) * 0.000001, nrdofs_glob);
	solver_p->rhs_vect = (double *)malloc(nrdofs_glob * sizeof(double));
	mf_check_mem(solver_p->rhs_vect);
	solver_p->sol_vect = (double *)malloc(nrdofs_glob * sizeof(double));
	mf_check_mem(solver_p->rhs_vect);
	solver_p->aux_vect = (int *)malloc(nrdofs_glob * sizeof(int));
	mf_check_mem(solver_p->rhs_vect);

#ifdef DEBUG_SIM
	/*kbw
	  printf("After allocating global SM and RHS: nrdofs_glob %d\n",
	   nrdofs_glob);
	  printf("solver_p->stiff_mat %lu, solver_p->rhs_vect %lu, solver_p->sol_vect %lu, aux_vect %lu\n",
	   solver_p->stiff_mat,solver_p->rhs_vect,solver_p->sol_vect, solver_p->aux_vect);
	/*kew*/
#endif



	for (i = 0; i < nrdofs_glob * nrdofs_glob; i++) solver_p->stiff_mat[i] = 0.0;
	for (i = 0; i < nrdofs_glob; i++) solver_p->rhs_vect[i] = 0.0;
	for (i = 0; i < nrdofs_glob; i++) solver_p->sol_vect[i] = 0.0;
	for (i = 0; i < nrdofs_glob; i++) solver_p->aux_vect[i] = 0;

	if (solver_p->aux_vect == NULL) {
#ifdef DEBUG_SIM
		printf("Not enough memory for LAPACK dense solver\n");
		printf("Required storage > %d MBytes\n",
				(int)((nrdofs_glob + 2)*nrdofs_glob * sizeof(double) / 1024 / 1024));
#endif
		exit(-1);

	}

	/*kbw
	printf("In sir_direct_init befor leaving: nr_int_ent %d, max_nrdofs %d\n",
	 solver_p->nr_int_ent, solver_p->max_dofs_int_ent);
	for(i=0;i<solver_p->nr_dof_ent;i++){
	  printf("block %d, nrdofs %d, posglob %d\n",
	   i, solver_p->l_dof_struct[i].nrdofs,
	   solver_p->l_dof_struct[i].posglob);
	}
	/*kew*/


	/* no additional parameters read from input file Filename */

	free(temp_list_dof_type);
	free(temp_list_dof_id);
	free(temp_list_dof_nrdofs);

	return (0);


}


/*------------------------------------------------------------
sir_direct_solve - to solve the system for a given data
------------------------------------------------------------*/
int sir_direct_solve( /* returns: >=0 - success code, <0 - error code */
		int Solver_id,     /* in: solver identification */
		int Comp_type,     /* in: indicator for the scope of computations: */
		/*   SIC_SOLVE - solve the system */
		/*   SIC_RESOLVE - resolve for the new right hand side */
		int Monitor       /* in: monitoring flag with options: */
		/*   SIC_PRINT_NOT - do not print anything */
		/*   SIC_PRINT_ERRORS - print error messages only */
		/*   SIC_PRINT_INFO - print most important information */
		/*   SIC_PRINT_ALLINFO - print all available information */
)
{

	/* pointer to solver structure */
	sit_direct_solver * solver_p;

	/* auxiliary variables */
	int nrdofs_glob, max_nrdofs, nr_int_ent, nr_dof_ent, posglob;
	int l_dof_ent_id[SIC_MAX_DOF_PER_INT], l_dof_ent_nrdofs[SIC_MAX_DOF_PER_INT];
	int l_dof_ent_posglob[SIC_MAX_DOF_PER_INT];
	int l_dof_ent_type[SIC_MAX_DOF_PER_INT];
	int pdr_comp_type;
	double * stiff_mat;
	double * rhs_vect;
	int i, j, k, iaux, kaux, intent, idofent, ibl;
	const int level_id = 0;
	int nrdofs_int_ent = 0;

	/*++++++++++++++++ executable statements ++++++++++++++++*/


	/* get pointer to solver structure */
	solver_p = sir_direct_select_solver(Solver_id);

	/* global number of degrees of freedom */
	nrdofs_glob = solver_p->nrdofs_glob;
	max_nrdofs = solver_p->max_dofs_int_ent;
	nr_int_ent = solver_p->nr_int_ent;

	/*kbw
	  printf("in sir_direct_solve: Solver_id %d, Comp_type %d\n",
	   Solver_id, Comp_type);
	  printf("nrdofs_glob %d, max_nrdofs %d, nr_int_ent %d\n",
	   nrdofs_glob, max_nrdofs, nr_int_ent);
	/*kew*/

	/* storage for local stiffness matrices  */
	stiff_mat = (double *)malloc(max_nrdofs * max_nrdofs * sizeof(double));
	rhs_vect = (double *)malloc(max_nrdofs * sizeof(double));


	if (Comp_type == SIC_SOLVE) pdr_comp_type = PDC_COMP_BOTH;
	else pdr_comp_type = PDC_COMP_RHS;

	/* compute local stiffness matrices */
	for (intent = 0; intent < nr_int_ent; intent++) {
		nr_dof_ent = SIC_MAX_DOF_PER_INT; /* in: size of arrays, */
		/* out: no of filled entries */
		nrdofs_int_ent = max_nrdofs;   /* in: size of arrays, */
		/* out: no of filled entries */
		pdr_comp_stiff_mat(solver_p->problem_id, solver_p->l_int_ent_type[intent],
				solver_p->l_int_ent_id[intent], pdr_comp_type, NULL,
				&nr_dof_ent, l_dof_ent_type, l_dof_ent_id, l_dof_ent_nrdofs,
				&nrdofs_int_ent, stiff_mat, rhs_vect, NULL);


#ifdef DEBUG_SIM
		if (nrdofs_int_ent > max_nrdofs) {
			printf("Too small arrays stiff_mat and rhs_vect passed to comp_el_stiff_mat\n");
			printf("from sir_direct_create in lapack. %d < %d. Exiting !!!", max_nrdofs, nrdofs_int_ent);
			exit(-1);
		}
#endif


		/*kbw
		{
		  printf("In sir_direct_solve before assemble: Solver_id %d, level_id %d, sol_typ %d\n", Solver_id, level_id, SIC_SOLVE);
		  int ibl,jbl,pli,plj,nri,nrj,nrdofs,jaux;
		printf("intent %d (type %d, id %d), nr_dof_ent %d\n",
		intent, solver_p->l_int_ent_type[intent],
		solver_p->l_int_ent_id[intent], nr_dof_ent);
		  pli = 0; nrdofs=0;
		  for(ibl=0;ibl<nr_dof_ent; ibl++) nrdofs+=l_dof_ent_nrdofs[ibl];
		  for(ibl=0;ibl<nr_dof_ent; ibl++){
		int dof_struct_id = solver_p->l_dof_vert_to_struct[l_dof_ent_id[ibl]];
		l_dof_ent_posglob[ibl]=solver_p->l_dof_struct[dof_struct_id].posglob;
		printf("bl_id %d, dof_struct %d, dof_ent %d, bl_nrdofs %d, bl_posglob %d\n",
		   solver_p->l_dof_vert_to_struct[l_dof_ent_id[ibl]],dof_struct_id,l_dof_ent_id[ibl],
			   l_dof_ent_nrdofs[ibl],l_dof_ent_posglob[ibl]);
		nri = l_dof_ent_nrdofs[ibl];
		plj=0;
		for(jbl=0;jbl<nr_dof_ent;jbl++){
		  printf("Stiff_mat (blocks %d:%d)\n",ibl,jbl);
		  nrj = l_dof_ent_nrdofs[jbl];
		  for(i=0;i<nri;i++){
			jaux = plj+(pli+i)*nrdofs;
			for(j=0;j<nrj;j++){
			  printf("%20.15lf",stiff_mat[jaux+j]);
			}
			printf("\n");
		  }
		  plj += nrj;
		}
		printf("Rhs_vect:\n");
		for(i=0;i<nri;i++){
		  printf("%20.15lf",rhs_vect[pli+i]);
		}
		printf("\n");
		pli += nri;
		  }
		  getchar();
		}
		/*kew*/

		for (idofent = 0; idofent < nr_dof_ent; idofent++) {
			int dof_ent_id = l_dof_ent_id[idofent];
			int dof_ent_type = l_dof_ent_type[idofent];
			int dof_struct_id;

			if (dof_ent_type == PDC_ELEMENT ) {

				dof_struct_id = solver_p->l_dof_elem_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				if ( solver_p->l_dof_elem_to_struct[dof_ent_id] == -1) {
					printf("Error 347294 in sir_direct_create!!! Exiting\n");
					exit(-1);
				}
#endif

			} else if (dof_ent_type == PDC_FACE ) {

				dof_struct_id = solver_p->l_dof_face_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				if ( solver_p->l_dof_face_to_struct[dof_ent_id] == -1) {
					printf("Error 347294 in sir_direct_create!!! Exiting\n");
					exit(-1);
				}
#endif

			} else if (dof_ent_type == PDC_EDGE ) {

				dof_struct_id = solver_p->l_dof_edge_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				if ( solver_p->l_dof_edge_to_struct[dof_ent_id] == -1) {
					printf("Error 347294 in sir_direct_create!!! Exiting\n");
					exit(-1);
				}
#endif

			} else if (dof_ent_type == PDC_VERTEX ) {

				dof_struct_id = solver_p->l_dof_vert_to_struct[dof_ent_id];

#ifdef DEBUG_SIM
				if ( solver_p->l_dof_vert_to_struct[dof_ent_id] == -1) {
					printf("Error 347294 in sir_direct_create!!! Exiting\n");
					exit(-1);
				}
#endif

			}


			l_dof_ent_posglob[idofent] = solver_p->l_dof_struct[dof_struct_id].posglob;

		}

		sir_direct_assemble_stiff_mat(Solver_id, Comp_type, nr_dof_ent,
				l_dof_ent_nrdofs, l_dof_ent_posglob,
				stiff_mat, rhs_vect, NULL);

	}

	/*
	printf("\n\n\n\nMain diagonal & RHS:\n");
	for(i=0;i<nrdofs_glob; i++){
	  //printf("%5d%3d%10.4lg", j,i,solver_p->stiff_mat[i*nrdofs_glob+j]);
	  printf("%7.3lf\t%7.3lf\n",solver_p->stiff_mat[i*nrdofs_glob+i], solver_p->rhs_vect[i]);
	}
	printf("\n\n\n\n");
	*/


	/*
		printf("LV\n");
	  for(i=0;i<nrdofs_glob; i++){
		//printf("%8d%10.4lg", i,solver_p->rhs_vect[i]);
		printf("%10.4lg",solver_p->rhs_vect[i]);
	  }
		printf("\n");
		printf("\n");



		printf("SM (stored by columns, displayed by rows)\n");
	  for(j=0;j<27; j++){
	  for(i=0;i<27; i++){
		//printf("%5d%3d%10.4lg", j,i,solver_p->stiff_mat[i*nrdofs_glob+j]);
		printf("%7.3lf\t",solver_p->stiff_mat[i*nrdofs_glob+j]);
	  }
		printf("\n");

	}
		printf("\n");


			printf("\n");
		printf("\n");
		printf("SM (stored by columns, displayed by rows)\n");
	  for(j=50;j<69; j++){
	  for(i=50;i<69; i++){
		//printf("%5d%3d%10.4lg", j,i,solver_p->stiff_mat[i*nrdofs_glob+j]);
		printf("%7.3lf\t",solver_p->stiff_mat[i*nrdofs_glob+j]);
	  }
		printf("\n");

	}
		printf("\n");


			printf("\n");
		printf("\n");
		printf("SM (stored by columns, displayed by rows)\n");
	  for(j=200;j<219; j++){
	  for(i=200;i<219; i++){
		//printf("%5d%3d%10.4lg", j,i,solver_p->stiff_mat[i*nrdofs_glob+j]);
		printf("%7.3lf\t",solver_p->stiff_mat[i*nrdofs_glob+j]);
	  }
		printf("\n");

	}
		printf("\n");



		printf("\n");
		*/



//            printf("\n");
//    printf("\n");
//    printf("SM (stored by columns, displayed by rows)\n");
//  for(j=1200;j<1219/*nrdofs_glob*/; j++){
//  for(i=1000;i<1019/*nrdofs_glob*/; i++){
//    //printf("%5d%3d%10.4lg", j,i,solver_p->stiff_mat[i*nrdofs_glob+j]);
//    printf("%7.3lf\t",solver_p->stiff_mat[i*nrdofs_glob+j]);
//  }
//    printf("\n");
//
//}
//    printf("\n");
	/*kew*/


	if (Comp_type == SIC_SOLVE) {
		dgetrf_(&nrdofs_glob, &nrdofs_glob, solver_p->stiff_mat, &nrdofs_glob,
				solver_p->aux_vect, &iaux);
		mf_check(iaux == 0, "DGETRF exited with error %d", iaux);
	}

	kaux = 1;
	dgetrs_("N", &nrdofs_glob, &kaux, solver_p->stiff_mat, &nrdofs_glob,
			solver_p->aux_vect, solver_p->rhs_vect, &nrdofs_glob, &iaux);

	/* rewrite the solution */
	for (ibl = 0; ibl < solver_p->nr_dof_ent; ibl++) {
		sit_direct_dof_struct dof_struct = solver_p->l_dof_struct[ibl];
		pdr_write_sol_dofs(solver_p->problem_id,
				dof_struct.dof_ent_type,
				dof_struct.dof_ent_id,
				dof_struct.nrdofs,
				&solver_p->rhs_vect[dof_struct.posglob]);



		/* print out solution
			printf("Solution in block %d, posglob %d, elem %d \n",ibl, dof_struct.posglob, dof_struct.dof_ent_id);
			for (i=0;i<dof_struct.nrdofs;i++)
			  printf("%20.10lf",solver_p->rhs_vect[dof_struct.posglob+i]);
			printf("\n");
		/**/
	}

	/* free the storage space */
	free(stiff_mat);
	free(rhs_vect);

	return (1);
}


/*------------------------------------------------------------
  sir_direct_assemble_stiff_mat - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
------------------------------------------------------------*/
int sir_direct_assemble_stiff_mat(
		/* returns: >=0 - success code, <0 - error code */
		int Solver_id,         /* in: solver ID (used to identify the subproblem) */
		int Comp_type,         /* in: indicator for the scope of computations: */
		/*   SIC_SOLVE - solve the system */
		/*   SIC_RESOLVE - resolve for the new rhs vector */
		int Nr_dof_bl,         /* in: number of global dof entities (blocks) */
		/*     associated with the local stiffness matrix */
		int * L_bl_nrdofs,      /* in: list of dof blocks' nr dofs */
		int * L_bl_posglob,    /* in: list of blocks' global positions */
		double * Stiff_mat,    /* in: stiffness matrix stored columnwise */
		double * Rhs_vect,     /* in: rhs vector */
		char * Rewr_dofs        /* in: flag to rewrite or sum up entries */
		/*   'T' - true, rewrite entries when assembling */
		/*   'F' - false, sum up entries when assembling */
)
{

	/* pointers to solver structure and stiffness matrix blocks' info */
	sit_direct_solver * solver_p;

	/* auxiliary variables */
	int iblock, jblock, nrdofs_i, nrdofs_j, posglob_i, posglob_j;
	int posloc_i, posloc_j, nrdofs, nrdofs_glob;
	int i, j, k, iaux, jaux;

	/*kbw
	static FILE *fp=NULL;
	/*kew*/

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/*kbw
	if(fp==NULL) fp=fopen("stiff_mat","a");
	/*kew*/

	/* get pointer to solver structure */
	solver_p = sir_direct_select_solver(Solver_id);

	mf_check_mem(solver_p->stiff_mat);

	/* global number of degrees of freedom */
	nrdofs_glob = solver_p->nrdofs_glob;

	/* compute local stiffness matrix nrdofs */
	nrdofs = 0;
	for (iblock = 0; iblock < Nr_dof_bl; iblock++) {
		nrdofs += L_bl_nrdofs[iblock];
	}

	/* loop over stiffness matrix blocks */
	posloc_i = 0;
	for (iblock = 0; iblock < Nr_dof_bl; iblock++) {

		/* global position of the first entry and number of dofs for a block */
		posglob_i = L_bl_posglob[iblock];
		nrdofs_i = L_bl_nrdofs[iblock];


		/*kbw
		//if(L_bl_posglob[iblock]==0){
		printf("in horizontal block_i %d, posglob_i %d, nrdofs_i %d\n",
		   iblock, L_bl_posglob[iblock], L_bl_nrdofs[iblock]);
		//}
		/*kew*/


		posloc_j = 0;
		for (jblock = 0; jblock < Nr_dof_bl; jblock++) {

			/* global position of the first entry and number of dofs for a block */
			posglob_j = L_bl_posglob[jblock];
			nrdofs_j = L_bl_nrdofs[jblock];

			/*kbw
			if(L_bl_posglob[jblock]==0&&L_bl_posglob[iblock]==11){
			printf("block_i %d, posglob_i %d, nrdofs_i %d\n",
			   iblock, L_bl_posglob[iblock], L_bl_nrdofs[iblock]);
			printf("block_j %d, posglob_i %d, nrdofs_i %d\n",
			   jblock, L_bl_posglob[jblock], L_bl_nrdofs[jblock]);
			}
			/*kew*/

			for (i = 0; i < nrdofs_i; i++) {

				/* global and local positions of the first entry in the column */
				iaux = posglob_j + (posglob_i + i) * nrdofs_glob;
				jaux = posloc_j + (posloc_i + i) * nrdofs;

				/*kbw
				if(L_bl_posglob[jblock]==0&&L_bl_posglob[iblock]==11){
				printf("assembling sm entries \nfrom local : %d - %d ",
					   jaux,jaux+nrdofs_j-1 );
				printf("to global %d - %d\nvalues (local  global): ",
					   iaux,iaux+nrdofs_j-1);
				}
				/*kew*/

				for (j = 0; j < nrdofs_j; j++) {

					/* assemble stiffness matrix block's entries */
					solver_p->stiff_mat[iaux + j] += Stiff_mat[jaux + j];

					/*kbw
					if(L_bl_posglob[jblock]==0&&L_bl_posglob[iblock]==11){
					printf("%20.15lf   %20.15lf\n",
					 Stiff_mat[jaux+j], solver_p->stiff_mat[iaux+j]);
					}
					/*kew*/

				}
			}

			/*kbw
			if(L_bl_posglob[jblock]==0&&L_bl_posglob[iblock]==11){
			printf("\n");
			}
			/*kew*/

			posloc_j += nrdofs_j;
		}

		/*kbw
		if(posglob_i==-1){
		printf("assembling rhsv entries \nfrom local :  %d - %d ",
		   posloc_i,posloc_i+nrdofs_i-1);
		printf("to global %d - %d\nvalues (local  global): ",
			   posglob_i,posglob_i+nrdofs_i-1);
		}
		/*kew*/

		for (i = 0; i < nrdofs_i; i++) {

			/* assemble right hand side block's entries */
			solver_p->rhs_vect[posglob_i + i] += Rhs_vect[posloc_i + i];

			/*kbw
			if(posglob_i==-1){
			printf("%20.15lf   %20.15lf\n",
			   Rhs_vect[posloc_i+i], solver_p->rhs_vect[posglob_i+i]);
			}
			/*kew*/

		}

		/*kbw
		if(posglob_i==-1){
		printf("\n");
		}
		/*kew*/

		posloc_i += nrdofs_i;
	}

	/*kbw
	getchar();
	/*kew*/


	return (1);
}


/*------------------------------------------------------------
sir_direct_free - to  free memory
------------------------------------------------------------*/
int sir_direct_free( /* returns: >=0 - success code, <0 - error code */
		int Solver_id   /* in: solver identification */
)
{

	/* pointers to the solver and its levels' structures - to simplify */
	sit_direct_solver * solver_p;


	/*++++++++++++++++ executable statements ++++++++++++++++*/

	solver_p = sir_direct_select_solver(Solver_id);

	free(solver_p->l_int_ent_type);
	free(solver_p->l_int_ent_id);
	if (solver_p->max_dof_elem_id >= 0) free(solver_p->l_dof_elem_to_struct);
	if (solver_p->max_dof_face_id >= 0) free(solver_p->l_dof_face_to_struct);
	if (solver_p->max_dof_edge_id >= 0) free(solver_p->l_dof_edge_to_struct);
	if (solver_p->max_dof_vert_id >= 0) free(solver_p->l_dof_vert_to_struct);
	free(solver_p->l_dof_struct);

	free(solver_p->stiff_mat);
	free(solver_p->rhs_vect);
	free(solver_p->sol_vect);
	free(solver_p->aux_vect);

	return (0);
}

/*------------------------------------------------------------
sir_direct_destroy - to destroy the solver
------------------------------------------------------------*/
int sir_direct_destroy( /* returns: >=0 - success code, <0 - error code */
		int Solver_id   /* in: solver identification */
)
{

	/* set the current solver ID */
	if (siv_direct_cur_solver_id == siv_direct_nr_solvers) siv_direct_cur_solver_id--;

	/* decrease the counter for solvers */
	siv_direct_nr_solvers--;

	return (1);
}


/*** AUXILIARY LOCAL PROCEDURES ***/

/*------------------------------------------------------------
  sir_direct_select_solver - to return the pointer to a given solver
------------------------------------------------------------*/
sit_direct_solver * sir_direct_select_solver( /* returns: pointer to solver */
		int Solver_id    /* in: solver identification */
)
{
	if (Solver_id > 0 && Solver_id <= siv_direct_nr_solvers) siv_direct_cur_solver_id = Solver_id;
	return (&siv_direct_solvers[siv_direct_cur_solver_id - 1]);
}


/*---------------------------------------------------------
sir_direct_put_list - to put Num on the list List with length Ll
	(filled with numbers and SIC_LIST_END_MARK at the end)
---------------------------------------------------------*/
int sir_direct_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
		/*  <0 - position at which put on the list */
		/*   0 - list full, not found on the list */
		int Num, 	/* in: number to put on the list */
		int * List, 	/* in: list */
		int Ll		/* in: total list's lengths */
)
{

	int i, il, pos;

	for (i = 0; i < Ll; i++) {
		pos = i;
		if ((il = List[pos]) == SIC_LIST_END_MARK) break;
		/* found on the list on (i+1) position; 1 offset - position 1 = List[0] */
		if (Num == il) return (++pos);
	}
	/* if list is full return error message */
	if (pos == Ll - 1) return (0);
	/* update the list and return*/
	List[pos] = Num;
	return (-(++pos));
}
