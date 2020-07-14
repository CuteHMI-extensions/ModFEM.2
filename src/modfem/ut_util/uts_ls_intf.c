/***********************************************************************
File uts_ls_intf - utility routines for interactions with linear solver
				   (common to all problem modules,
		   possibly used also by other modules)

Contains definitions of routines:

from main include/uth_intf.h
  utr_get_list_ent - to return the list of:
			  1. integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module,
		  2. DOF entities - entities with which there are dofs
			 associated by the given approximation

  utr_get_list_ent_coarse - to return the list of integration entities - entities
						 for which stiffness matrices and load vectors are
						 provided by the FEM code to the solver module,
						 and DOF entities - entities with which there are dofs
						 associated by the given approximation for COARSE level
						 given the corresponding lists from the fine level
			  (it calls separate implementations for dg and std?)

  utr_get_max_num_grid_levels - to limit the number of levels in multigrid
						 based on mesh and field data

  utr_dof_ent_sons - to return info whether the entity is owned
					 and a list of dof entity sons for owned entity


  utr_comp_stiff_mat - to create a stiffness matrix
					  and a load vector corresponding to the specified
					  mesh entity, together with information on how to
					  assemble entries into the global stiffness matrix
					  and the global load vector

  utr_create_assemble_stiff_mat - to create element stiffness matrices
								 and assemble them to the global SM
			  (it calls implementations for particular platforms)
  utr_rewr_sol - to rewrite solution from one vector to another

  utr_renumber - to renumber (permute) graph vertices (used for bandwidth reduction)
  utr_color_int_ent_for_assembly_with_graph_creation - coloring elements for lock free
					  assembly (with creation of neighbourhood graphs)
  utr_color_int_ent_for_assembly - coloring elements for lock free assembly


------------------------------
History:
	08.2008 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<signal.h>

#ifdef _OPENMP
#include<omp.h>
#endif

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>

/* interface for all approximation modules */
#include <modfem/aph_intf.h>

#include <modfem/uth_system.h>

#ifdef PARALLEL
/* interface for parallel mesh manipulation modules */
#include <modfem/mmph_intf.h>

/* interface with parallel communication library */
#include <modfem/pch_intf.h>
#endif

/* interface for general purpose utilities - for all problem dependent modules*/
#include <modfem/uth_intf.h>

#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)
  // interface for accelerator utilities - for all problem dependent modules*/
  #include "./uth_accel_intf.h"
#endif

/* interface for linear algebra packages */
#include <modfem/lin_alg_intf.h>

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>

#include <modfem/pdh_intf.h>
#include <modfem/sih_intf.h>
#include <modfem/uth_log.h>

// maximal number of solution components
#define UTC_MAXEQ PDC_MAXEQ


/**-----------------------------------------------------------
  utr_rewr_sol - to rewrite solution from one vector to another
------------------------------------------------------------*/
int utr_rewr_sol(
		 int Field_id,      /** in: data structure to be used  */
		 int Sol_from,      /** in: ID of vector to read solution from */
		 int Sol_to         /** in: ID of vector to write solution to */
		 )
{
  // Rewrite solution on host side
  apr_rewr_sol(Field_id,Sol_from,Sol_to);

  // Rewrite solution on device side
  #if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)



  #endif

}


/*------------------------------------------------------------
  utr_get_list_int_ent - to return the list of integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module
------------------------------------------------------------*/
int utr_get_list_int_ent( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in:  problem (and solver) identification */
  int* Nr_int_ent,    /* out: number of integration entitites */
  int** List_int_ent_type,/* out: list of types of integration entitites */
  int** List_int_ent_id  /* out: list of IDs of integration entitites */
			  )
{

  char field_module_name[100];
  int solver_id, field_id, mesh_id, nreq;
  int nel, nfa, nno, int_ent, dof_ent, nr_elem, nr_face, nr_node, i, iaux, iedg;

  double time_begin = time_clock();

  /* associate the suitable field with the problem and the solver */
  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; field_id = pdr_ctrl_i_params(Problem_id,i);
  i=6; solver_id = pdr_ctrl_i_params(Problem_id,i);
  //mesh_id = apr_get_mesh_id(field_id);
  nreq=apr_get_nreq(field_id);

  // check the name of the field module
  apr_module_introduce(field_module_name);

  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) != 0 &&
	  strncmp(field_module_name,"STANDARD_QUADRATIC",18) != 0 ){
	printf("utr_get_list_int_ent works only for STANDARD (linear and quadratic) approximations!\n");
	printf("Exiting.\n");
	exit(-1);
  }

  /* prepare arrays */
  nr_elem = mmr_get_nr_elem(mesh_id);
  nr_face = mmr_get_nr_face(mesh_id);

  *List_int_ent_type = (int *) malloc( (nr_elem+nr_face+1)*sizeof(int) );
  *List_int_ent_id = (int *) malloc( (nr_elem+nr_face+1)*sizeof(int) );

  int_ent=0;

  /* loop over elements - integration entities*/
  nel=0;
  while((nel=mmr_get_next_elem_all(mesh_id, nel))!=0){

	if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE ) {

	  int_ent++;
	  (*List_int_ent_type)[int_ent-1]=PDC_ELEMENT;
	  (*List_int_ent_id)[int_ent-1]=nel;

	} // end if element included for integration

  } // end for all elements

  /*loop over faces*/
  nfa=0;
  while((nfa=mmr_get_next_face_all(mesh_id, nfa))!=0) {

	if(mmr_fa_status(mesh_id,nfa)==MMC_ACTIVE) {
	  if(mmr_fa_bc(mesh_id, nfa)>0) {

	int_ent++;
	(*List_int_ent_type)[int_ent-1]=PDC_FACE;
	(*List_int_ent_id)[int_ent-1]=nfa;

	  }
	}
  }

  *Nr_int_ent = int_ent;

  utv_time.create_int_ent_list += time_clock()-time_begin;

  return(0);
}

/*------------------------------------------------------------
  utr_get_list_dof_ent - to return the list of DOF entities
						 - entities with which there are dofs
			 associated by the given approximation
						 given the list of integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module
------------------------------------------------------------*/
int utr_get_list_dof_ent( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in:  problem (and solver) identification */
  int* Nr_int_ent,    /* in: number of integration entitites */
  int** List_int_ent_type,/* in: list of types of integration entitites */
  int** List_int_ent_id,  /* in: list of IDs of integration entitites */
	/* GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!! */
  int* Nr_dof_ent,    /* out: number of dof entities (entities with which there
			  are dofs associated by the given approximation) */
  int** List_dof_ent_type,/* out: list of types of integration entitites */
  int** List_dof_ent_id,  /* out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /* out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent/* out: maximal number of dofs per dof entity */
			  )
{

  char field_module_name[100];
  int solver_id, field_id, mesh_id, nreq;
  int nel, nfa, nno, ino, int_ent, dof_ent, nr_elem, nr_face, nr_edge, nr_node;
  int i, iaux, iedg, ifac, ielm;
  int nrdofsgl, nrdofent, face_neig[2];
  int nrdofsloc, max_dofs_ent_dof;
  // i.e. maximal number of DOF blocks (nodes) per element
  int l_dof_ent_types[SIC_MAX_DOF_PER_INT];
  int l_dof_ent_ids[SIC_MAX_DOF_PER_INT];
  int l_dof_ent_nrdofs[SIC_MAX_DOF_PER_INT];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  double time_begin = time_clock();

  /* associate the suitable field with the problem and the solver */
  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; field_id = pdr_ctrl_i_params(Problem_id,i);
  i=6; solver_id = pdr_ctrl_i_params(Problem_id,i);
  //mesh_id = apr_get_mesh_id(field_id);
  nreq=apr_get_nreq(field_id); // it is assumed that the first field has higher pdeg

  // check the name of the field module
  apr_module_introduce(field_module_name);

  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) != 0 &&
	  strncmp(field_module_name,"STANDARD_QUADRATIC",18) != 0 ){

	printf("utr_get_list_dof_ent works only for STANDARD (linear and quadratic) approximations!\n");
	printf("Exiting.\n");
	exit(-1);
  }

  // for standard continuous linear FEM approximation
  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0){

	// list to indicate whether the entity has been already taken into account
	int *temp_list_dofs;
	iaux = mmr_get_max_node_id(mesh_id);
	temp_list_dofs = (int *) malloc( (iaux+1)*sizeof(int) );
	for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;

	/* prepare arrays */
	nr_node= mmr_get_nr_node(mesh_id);


	*List_dof_ent_type = (int *) malloc( (nr_node+1)*sizeof(int) );
	*List_dof_ent_id = (int *) malloc( (nr_node+1)*sizeof(int) );
	*List_dof_ent_nrdofs = (int *) malloc( (nr_node+1)*sizeof(int) );

	if((*List_dof_ent_nrdofs)==NULL){
	  printf("Not enough space for allocating lists in utr_get_list_ent! Exiting\n");
	  exit(-1);
	}

	nrdofsgl = 0;

	dof_ent=0;

	/*loop over vertexes - dof entities*/
	//nno=0;
	//while((nno=mmr_get_next_node_all(mesh_id,nno))!=0) {
	// we consider entities that appear in integrals computed for integration entities'
	for(int_ent=0; int_ent<*Nr_int_ent;int_ent++){

	  // it is sufficient to take into account elements only
	  if((*List_int_ent_type)[int_ent]==PDC_ELEMENT){

	int nel = (*List_int_ent_id)[int_ent];
	int nr_dof_ent_loc = SIC_MAX_DOF_PER_INT; // maximal number of DOF ents per INT ent

	pdr_comp_stiff_mat(Problem_id, PDC_ELEMENT, nel, PDC_NO_COMP, NULL,
			   &nr_dof_ent_loc, l_dof_ent_types,
			   l_dof_ent_ids, l_dof_ent_nrdofs,
			   NULL, NULL, NULL, NULL);

	// check whether there are owned DOFs (vertices=nodes)
	int ino;
	for(ino=0;ino<nr_dof_ent_loc;ino++){

		int nno = l_dof_ent_ids[ino];

		// for active DOF entities
		if(apr_get_ent_pdeg(field_id,APC_VERTEX,nno)>0 ){

#ifdef PARALLEL
		  /* in first round only internal vertices=nodes are considered */
		  if(mmpr_ve_owner(mesh_id,nno)==pcr_my_proc_id()){
#endif

		// if not yet considered
		if(temp_list_dofs[nno]==0){

		  // mark as considered
		  temp_list_dofs[nno]=1;

		  (*List_dof_ent_type)[dof_ent] = PDC_VERTEX;
		  (*List_dof_ent_id)[dof_ent]=nno;
		  (*List_dof_ent_nrdofs)[dof_ent]=nreq;
		  dof_ent++;
		  nrdofsgl+=nreq;

/*kbw
  printf("In pdr_get_list_ent: nno %d, pdeg %d, dof_ent %d, nrdofsgl %d\n",
  nno,apr_get_ent_pdeg(field_id,APC_VERTEX,nno),dof_ent, nrdofsgl);
/*kew*/


		} // if entity not yet considered
#ifdef PARALLEL
		  } // end if owned entity
#endif
		} // end if active entity
	  } // end if DOF entity on a list from integration entity
	} // end if element
	  } // end loop over integration entities



/*kbw
  printf("After first round pdr_get_list_ent - internal nodes\n");
  printf("\nNr_dof_ent %d, Nrdofs_glob %d, Max_dofs_per_dof_ent %d\n",
  dof_ent, nrdofsgl, nreq);
  for(i=0;i<dof_ent;i++)  printf("type %d, id %d, nrdofs %d\n",
  (*List_dof_ent_type)[i],(*List_dof_ent_id)[i],(*List_dof_ent_nrdofs)[i]);
//kew*/

#ifdef PARALLEL
	// second round
	{

	  //nel=0;
	  //while((nel=mmr_get_next_act_elem(mesh_id, nel))!=0){

	  // for all elements (i.e. integration entities being elements)
	  for(int_ent=0; int_ent<*Nr_int_ent;int_ent++){

	if((*List_int_ent_type)[int_ent]==PDC_ELEMENT){

	  nel = (*List_int_ent_id)[int_ent];

	  int nr_dof_ent_loc = SIC_MAX_DOF_PER_INT; // maximal number of DOF ents per INT ent
	  int is_in_overlap=0;
	  pdr_comp_stiff_mat(Problem_id, PDC_ELEMENT, nel, PDC_NO_COMP, NULL,
				 &nr_dof_ent_loc, l_dof_ent_types,
				 l_dof_ent_ids, l_dof_ent_nrdofs,
				 NULL, NULL, NULL, NULL);

	  // check whether there are not owned DOFs (vertices=nodes)
	  for(ino=0;ino<nr_dof_ent_loc;ino++){

		nno = l_dof_ent_ids[ino];

		if(mmpr_ve_owner(mesh_id,nno)!=pcr_my_proc_id()){
		  is_in_overlap = 1;
		  break;
		}

	  }
	  if(is_in_overlap){

		for(ino=0;ino<nr_dof_ent_loc;ino++){

		  nno = l_dof_ent_ids[ino];

		  if(temp_list_dofs[nno]==0){

		temp_list_dofs[nno]=1;

		// GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!!
		//(*List_dof_ent_type)[dof_ent]=PDC_VERTEX;
		(*List_dof_ent_type)[dof_ent] = -PDC_VERTEX;
		(*List_dof_ent_id)[dof_ent]=nno;
		(*List_dof_ent_nrdofs)[dof_ent]=nreq;
		dof_ent++;
		nrdofsgl+=nreq;



		  }
		}
	  } // end if element in overlap
	} // end if element
	  } // end loop over integration entities
	} // end of second round

#endif

	*Nrdofs_glob = nrdofsgl;
	*Max_dofs_per_dof_ent = nreq;
	*Nr_dof_ent = dof_ent;

/*kbw
  printf("In pdr_get_list_ent for standard linear approximation\n");
  printf("nrelem %d, nrface %d, nr_int_ent %d\n",
  nr_elem, nr_face, *Nr_int_ent);
  for(i=0;i<int_ent;i++)  printf("type %d, id %d\n",
  (*List_int_ent_type)[i],(*List_int_ent_id)[i]);
  printf("\nNr_dof_ent %d, Nrdofs_glob %d, Max_dofs_per_dof_ent %d\n",
  *Nr_dof_ent, *Nrdofs_glob, *Max_dofs_per_dof_ent);
  for(i=0;i<dof_ent;i++)  printf("type %d, id %d, nrdofs %d\n",
  (*List_dof_ent_type)[i],(*List_dof_ent_id)[i],(*List_dof_ent_nrdofs)[i]);
  /*kew*/

	free(temp_list_dofs);

  } // end if standard continuous quadratic FEM approximation
  else if(!strncmp(field_module_name,"STANDARD_QUADRATIC",18)) {

	/* auxiliary variable */
	int nne, nnf, nnel;
	int pdeg, geo_pdeg;

	// get pdeg and geo_pdeg value
	pdeg = apr_get_el_pdeg(field_id,0,NULL);
	geo_pdeg = apr_get_el_geo_pdeg(field_id,0,NULL);

	// lists to indicate whether the entity has been already taken into account
	int* temp_list_dofs = NULL, * temp_list_edges = NULL, * temp_list_faces = NULL, * temp_list_elems = NULL;

	//
	// LINEAR APPROXIMATION (LINEAR OR QUADRATIC GEOMETRIC SHAPE FUNCTIONS)
	//
	if(pdeg == APC_LINEAR_APPROXIMATION_PDEG)
	  {
	/* prepare auxilliary structures */
	iaux = mmr_get_max_node_id(mesh_id);
	temp_list_dofs = (int *) malloc( (iaux+1)*sizeof(int) );
	for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;

	/* prepare arrays */
	nr_node = mmr_get_nr_node(mesh_id);

	*List_dof_ent_type = (int *) malloc( (nr_node+1)*sizeof(int) );
	*List_dof_ent_id = (int *) malloc( (nr_node+1)*sizeof(int) );
	*List_dof_ent_nrdofs = (int *) malloc( (nr_node+1)*sizeof(int) );
	  }

	//
	// HIERARCHICAL QUADRATIC APPROXIMATION (LINEAR GEOMETRIC SHAPE FUNCTIONS)
	//
	if(pdeg == APC_QUADRATIC_HIERACHICAL_APPROXIMATION_PDEG)
	  {
	/* prepare auxilliary structures */
	iaux = mmr_get_max_node_id(mesh_id);
	temp_list_dofs = (int *) malloc( (iaux+1)*sizeof(int) );
	for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;

	iedg = mmr_get_max_edge_id(mesh_id);
	temp_list_edges = (int *) malloc( (iedg+1)*sizeof(int) );
	for(i=0;i<=iedg;i++) temp_list_edges[i]=0;

	/* prepare arrays */
	nr_node = mmr_get_nr_node(mesh_id);
	nr_edge = mmr_get_nr_edge(mesh_id);

	*List_dof_ent_type = (int *) malloc( (nr_node+nr_edge+1)*sizeof(int) );
	*List_dof_ent_id = (int *) malloc( (nr_node+nr_edge+1)*sizeof(int) );
	*List_dof_ent_nrdofs = (int *) malloc( (nr_node+nr_edge+1)*sizeof(int) );
	  }



	//
	// QUADRATIC APPROXIMATION (QUADRATIC GEOMETRIC SHAPE FUNCTIONS)
	//
	if(pdeg == APC_QUADRATIC_APPROXIMATION_PDEG)
	  {
	/* prepare auxilliary structures */
	iaux = mmr_get_max_node_id(mesh_id);
	temp_list_dofs = (int *) malloc( (iaux+1)*sizeof(int) );
	for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;

	iedg = mmr_get_max_edge_id(mesh_id);
	temp_list_edges = (int *) malloc( (iedg+1)*sizeof(int) );
	for(i=0;i<=iedg;i++) temp_list_edges[i]=0;

	ifac = mmr_get_max_face_id(mesh_id);
	temp_list_faces = (int *) malloc( (ifac+1)*sizeof(int) );
	for(i=0;i<=ifac;i++) temp_list_faces[i]=0;

	/* prepare arrays */
	nr_node = mmr_get_nr_node(mesh_id);
	nr_edge = mmr_get_nr_edge(mesh_id);

	// WARNING: not all faces (only quadratic) have DOFs associated with them
	nr_face = mmr_get_nr_face(mesh_id);

	*List_dof_ent_type = (int *) malloc( (nr_node+nr_edge+nr_face+1)*sizeof(int) );
	*List_dof_ent_id = (int *) malloc( (nr_node+nr_edge+nr_face+1)*sizeof(int) );
	*List_dof_ent_nrdofs = (int *) malloc( (nr_node+nr_edge+nr_face+1)*sizeof(int) );

	  }

	//
	// MIXED APPROXIMATION
	// - LINEAR FIELD + CONSTANT FIELD (LINEAR GEOMETRIC SHAPE FUNCTIONS)
	// - QUADRATIC FIELD + LINEAR FIELD (QUADRATIC GEOMETRIC SHAPE FUNCTIONS)
	//

	// Second field ID
	int second_field_id = -1;

	// Mixed pdeg
	int mixed_pdeg_field_1 = -1; int mixed_pdeg_field_2 = -1;
	int mixed_nreq_field_1 = -1; int mixed_nreq_field_2 = -2;

	if(pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
	   pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
	   pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG)
	  {
	// Get second field ID
	second_field_id = apr_get_mixed_second_field_id(field_id);

	// Get mixed pdeg values
	apr_get_mixed_pdeg(field_id,&mixed_pdeg_field_1);
	apr_get_mixed_pdeg(second_field_id,&mixed_pdeg_field_2);

	// Get nreq values for mixed fields
	mixed_nreq_field_1 = apr_get_nreq(field_id);
	mixed_nreq_field_2 = apr_get_nreq(second_field_id);

	if(mixed_pdeg_field_2 > mixed_pdeg_field_1){
	  printf("higher pdeg in second field for mixed - check the code!\n");
	  exit(-1);
	}
	if(mixed_nreq_field_2 > mixed_nreq_field_1){
	  printf("higher nreq in second field for mixed - check the code!\n");
	  exit(-1);
	}


	/*jbw
	printf("mixed_nreq_field_1 = %d ; mixed_nreq_field_2 = %d\n",mixed_nreq_field_1,mixed_nreq_field_2);
	getchar();
	/*jbw*/

	//
	// MIXED APPROXIMATION - Field 1: LINEAR / Field 2: LINEAR
	// Very artificial combination for code testing
	//
	if(mixed_pdeg_field_1 == APC_LINEAR_APPROXIMATION_PDEG &&
	   mixed_pdeg_field_2 == APC_LINEAR_APPROXIMATION_PDEG)
	  {
		/* prepare auxilliary structures */
		iaux = mmr_get_max_node_id(mesh_id);
		iaux = iaux+iaux;
		temp_list_dofs = (int *) malloc((iaux+1)*sizeof(int) );
		for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;

		/* prepare arrays */
		nr_node = mmr_get_nr_node(mesh_id);

		*List_dof_ent_type = (int *) malloc( (nr_node+nr_node+1)*sizeof(int) );
		*List_dof_ent_id = (int *) malloc( (nr_node+nr_node+1)*sizeof(int) );
		*List_dof_ent_nrdofs = (int *) malloc( (nr_node+nr_node+1)*sizeof(int) );
	  }

	//
	// MIXED APPROXIMATION - Field 1: QUADRATIC / Field 2: LINEAR
	//
	if(mixed_pdeg_field_1 == APC_QUADRATIC_APPROXIMATION_PDEG &&
	   mixed_pdeg_field_2 == APC_LINEAR_APPROXIMATION_PDEG)
	  {
		/* prepare auxilliary structures */
		iaux = mmr_get_max_node_id(mesh_id);
		iaux = iaux+iaux;
		temp_list_dofs = (int *) malloc( (iaux+1)*sizeof(int) );
		for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;

		iedg = mmr_get_max_edge_id(mesh_id);
		temp_list_edges = (int *) malloc( (iedg+1)*sizeof(int) );
		for(i=0;i<=iedg;i++) temp_list_edges[i]=0;

		ifac = mmr_get_max_face_id(mesh_id);
		temp_list_faces = (int *) malloc( (ifac+1)*sizeof(int) );
		for(i=0;i<=ifac;i++) temp_list_faces[i]=0;

		/* prepare arrays */
		nr_node = mmr_get_nr_node(mesh_id);
		nr_edge = mmr_get_nr_edge(mesh_id);
		nr_face = mmr_get_nr_face(mesh_id);

		*List_dof_ent_type = (int *) malloc( (nr_node+nr_node+nr_edge+nr_face+1)*sizeof(int) );
		*List_dof_ent_id = (int *) malloc( (nr_node+nr_node+nr_edge+nr_face+1)*sizeof(int) );
		*List_dof_ent_nrdofs = (int *) malloc( (nr_node+nr_node+nr_edge+nr_face+1)*sizeof(int) );
	  }

	//
	// MIXED APPROXIMATION - Field 1: QUADRATIC / Field 2: QUADRATIC
	//
	if(mixed_pdeg_field_1 == APC_QUADRATIC_APPROXIMATION_PDEG &&
	   mixed_pdeg_field_2 == APC_QUADRATIC_APPROXIMATION_PDEG)
	  {
		/* prepare auxilliary structures */
		iaux = mmr_get_max_node_id(mesh_id);
		iaux = iaux+iaux;
		temp_list_dofs = (int *) malloc( (iaux+1)*sizeof(int) );
		for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;

		iedg = mmr_get_max_edge_id(mesh_id);
		iedg = iedg+iedg;
		temp_list_edges = (int *) malloc( (iedg+1)*sizeof(int) );
		for(i=0;i<=iedg;i++) temp_list_edges[i]=0;

		ifac = mmr_get_max_face_id(mesh_id);
		ifac = ifac+ifac;
		temp_list_faces = (int *) malloc( (ifac+1)*sizeof(int) );
		for(i=0;i<=ifac;i++) temp_list_faces[i]=0;

		/* prepare arrays */
		nr_node = mmr_get_nr_node(mesh_id);
		nr_edge = mmr_get_nr_edge(mesh_id);
		nr_face = mmr_get_nr_face(mesh_id);

		*List_dof_ent_type = (int *) malloc( (nr_node+nr_node+nr_edge+nr_edge+nr_face+nr_face+1)*sizeof(int) );
		*List_dof_ent_id = (int *) malloc( (nr_node+nr_node+nr_edge+nr_edge+nr_face+nr_face+1)*sizeof(int) );
		*List_dof_ent_nrdofs = (int *) malloc( (nr_node+nr_node+nr_edge+nr_edge+nr_face+nr_face+1)*sizeof(int) );
	  }
	  }

	// Check list of DOFS exit if list is NULL space
	if((*List_dof_ent_nrdofs)==NULL){
	  printf("Not enough space for allocating lists in utr_get_list_ent! Exiting\n");
	  exit(-1);
	}


	// ---------------------------------------------
	// -------- LOOP OVER DOFS ----------
	// ---------------------------------------------

	nrdofsgl = 0;
	dof_ent=0;

	// LINEAR AND QUADRATIC APPROXIMATION
	if(pdeg == APC_LINEAR_APPROXIMATION_PDEG ||
	   pdeg == APC_QUADRATIC_HIERACHICAL_APPROXIMATION_PDEG ||
	   pdeg == APC_QUADRATIC_APPROXIMATION_PDEG)
	  {

	/*loop over vertexes - dof vertex entities*/
	//nno=0;
	//while((nno=mmr_get_next_node_all(mesh_id,nno))!=0) {

	// we consider entities that appear in integrals computed for integration entities'
	for(int_ent=0; int_ent<*Nr_int_ent;int_ent++) {

	  // it is sufficient to take into account elements only
	  if((*List_int_ent_type)[int_ent]==PDC_ELEMENT){

		int nel = (*List_int_ent_id)[int_ent];
		int nr_dof_ent_loc = SIC_MAX_DOF_PER_INT; // maximal number of DOF ents per INT ent

		pdr_comp_stiff_mat(Problem_id, PDC_ELEMENT, nel, PDC_NO_COMP, NULL,
				   &nr_dof_ent_loc, l_dof_ent_types,
				   l_dof_ent_ids, l_dof_ent_nrdofs,
				   NULL, NULL, NULL, NULL);

		// check whether there are owned DOFs (vertices=nodes)
		int idofent;
		for(idofent=0;idofent<nr_dof_ent_loc;idofent++){

		  // if vertex
		  if(l_dof_ent_types[idofent]==PDC_VERTEX){

		int nno = l_dof_ent_ids[idofent];

		if(apr_get_ent_pdeg(field_id,APC_VERTEX,nno)>0 ){

#ifdef PARALLEL
		  /* in first round only internal vertices=nodes are considered */
		  if(mmpr_ve_owner(mesh_id,nno)==pcr_my_proc_id()){
#endif

			if(temp_list_dofs[nno]==0){

			  temp_list_dofs[nno]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_VERTEX;
			  (*List_dof_ent_id)[dof_ent]=nno;
			  (*List_dof_ent_nrdofs)[dof_ent]=nreq;
			  dof_ent++;
			  nrdofsgl+=nreq;

			  /* jbw
			 printf("utr_get_list_int_ent: dof_ent type -> PDC_VERTEX ; dof = %d ; dof_ent = %d\n",nno,dof_ent);
			 /* jbw */

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal vertex */
#endif
		} // end if active entity
		  } // end if vertex

		  // if edge
		  else if(l_dof_ent_types[idofent]==PDC_EDGE){

		if(pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
		  mf_fatal_err("In linear approximation should not be EDGE in DOF structure");
		}

		int nne = l_dof_ent_ids[idofent];

		if(apr_get_ent_pdeg(field_id,APC_EDGE,nne)>0){

#ifdef PARALLEL
		  /* in first round only internal edges are considered */
		  if(mmpr_ed_owner(mesh_id,nne)==pcr_my_proc_id()){
#endif

			if(temp_list_edges[nne]==0){

			  temp_list_edges[nne]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_EDGE;
			  (*List_dof_ent_id)[dof_ent]=nne;
			  (*List_dof_ent_nrdofs)[dof_ent]=nreq;
			  dof_ent++;
			  nrdofsgl+=nreq;

			  /* jbw
			 printf("utr_get_list_int_ent: dof_ent type -> APC_EDGE ; dof = %d ; dof_ent = %d\n",nne,dof_ent);
			 /* jbw */

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal edge */
#endif
		} // end if active entity
		  } // end if EDGE

		  // if face
		  else if(l_dof_ent_types[idofent]==PDC_FACE){

		if(pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
		  mf_fatal_err("In linear approximation should not be FACE in DOF structure");
		}

		int nnf = l_dof_ent_ids[idofent];

		if(apr_get_ent_pdeg(field_id,APC_FACE,nnf)>0){

#ifdef PARALLEL
		  /* in first round only internal faces are considered */
		  if(mmpr_ed_owner(mesh_id,nnf)==pcr_my_proc_id()){
#endif

			if(temp_list_faces[nnf]==0){

			  temp_list_faces[nnf]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_FACE;
			  (*List_dof_ent_id)[dof_ent]=nnf;
			  (*List_dof_ent_nrdofs)[dof_ent]=nreq;
			  dof_ent++;
			  nrdofsgl+=nreq;

			  /* jbw
			 printf("utr_get_list_int_ent: dof_ent type -> APC_FACE ; dof = %d ; dof_ent = %d\n",nnf,dof_ent);
			 /* jbw */

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal face */
#endif
		} // end if active entity
		  } // end if FACE

		} // end if DOF entity on a list from integration entity
	  } // end if element
	} // end loop over integration entities

	  } // END LINEAR AND QUADRATIC APPROXIMATION

	// MIXED APPROXIMATION
	else if(pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG ||
		pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
		pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG)
	  {

	// Configured in CMake
	//#define MIXED_FIELD_BY_FIELD
	//#define MIXED_DOF_BY_DOF
	#ifdef MIXED_FIELD_BY_FIELD
	printf("- MIXED APPROXIMATION storage FIELD by FIELD is ENABLED\n");
	#else
	printf("- MIXED APPROXIMATION storage DOF by DOF is ENABLED\n");
	#endif

	// we consider entities that appear in integrals computed for integration entities'
	for(int_ent=0; int_ent<*Nr_int_ent;int_ent++) {

	  // it is sufficient to take into account elements only
	  if((*List_int_ent_type)[int_ent]==PDC_ELEMENT){

		int nel = (*List_int_ent_id)[int_ent];
		int nr_dof_ent_loc = SIC_MAX_DOF_PER_INT; // maximal number of DOF ents per INT ent

		pdr_comp_stiff_mat(Problem_id, PDC_ELEMENT, nel, PDC_NO_COMP, NULL,
				   &nr_dof_ent_loc, l_dof_ent_types,
				   l_dof_ent_ids, l_dof_ent_nrdofs,
				   NULL, NULL, NULL, NULL);

		// check whether there are owned DOFs (vertices=nodes)
		int idofent;

		/*jbw
		for(idofent=0;idofent<nr_dof_ent_loc;idofent++) {
		  printf("(utr_get_list_int_ent) l_dof_ent_types[%d]=%d ; l_dof_ent_ids[%d]=%d ; l_dof_ent_nrdofs[%d]=%d\n",
			 idofent,l_dof_ent_types[idofent],
			 idofent,l_dof_ent_ids[idofent],
			 idofent,l_dof_ent_nrdofs[idofent]);
		}
		getchar();
		/*jbw*/

		for(idofent=0;idofent<nr_dof_ent_loc;idofent++) {

		  // if vertex
		  if(l_dof_ent_types[idofent] == PDC_VERTEX){

		int nno = l_dof_ent_ids[idofent];

		if(apr_get_ent_pdeg(field_id,APC_VERTEX,nno)>0 ){


#ifdef PARALLEL
		  /* in first round only internal vertices=nodes are considered */
		  if(mmpr_ve_owner(mesh_id,nno)==pcr_my_proc_id()){
#endif
				if(temp_list_dofs[nno]==0){

			  temp_list_dofs[nno]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_VERTEX;
			  (*List_dof_ent_id)[dof_ent]=nno;
			  (*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_1;
			  dof_ent++;
			  nrdofsgl+=mixed_nreq_field_1;

			  /* jbw
			  printf("utr_get_list_int_ent: dof_ent type -> PDC_VERTEX ; dof = %d ; [ dof_ent = %d ]\n",
				 nno,dof_ent);
			  /* jbw */

#ifdef MIXED_DOF_BY_DOF
			  // we add second vertex for second field
			  if(pdeg == APC_MIXED_P2_P1_APPROXIMATION_PDEG ||
			 pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG ||
			 pdeg == APC_MIXED_P1_P1_APPROXIMATION_PDEG){


			(*List_dof_ent_type)[dof_ent]=PDC_MIXED_VERTEX;
			(*List_dof_ent_id)[dof_ent]=nno;
			(*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_2;
			dof_ent++;
			nrdofsgl+=mixed_nreq_field_2;

			/* jbw
			printf("utr_get_list_int_ent: dof_ent type -> PDC_MIXED_VERTEX ; dof = %d ; [ dof_ent = %d ]\n",
				   nno,dof_ent);
			/* jbw */

					  }

#endif // end if DOF_BY_DOF

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal vertex */
#endif



		} // end if active entity
		  } // end if vertex

		  // if vertex mixed
		  else if(l_dof_ent_types[idofent] == PDC_MIXED_VERTEX ){

#ifdef MIXED_FIELD_BY_FIELD

		int nno_mixed_shift = mmr_get_max_node_id(mesh_id);
		int nno = l_dof_ent_ids[idofent];
		int loc_nno = l_dof_ent_ids[idofent]+nno_mixed_shift;

		if(apr_get_ent_pdeg(second_field_id,APC_MIXED_VERTEX,nno)>0 ){


#ifdef PARALLEL
		  /* in first round only internal vertices=nodes are considered */
		  if(mmpr_ve_owner(mesh_id,nno)==pcr_my_proc_id()){
#endif
			if(temp_list_dofs[loc_nno]==0){

			  temp_list_dofs[loc_nno]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_MIXED_VERTEX ;
			  (*List_dof_ent_id)[dof_ent]=nno;
			  (*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_2;
			  dof_ent++;
			  nrdofsgl+=mixed_nreq_field_2;

			  /* jbw
			  printf("utr_get_list_int_ent: dof_ent type -> APC_MIXED_VERTEX ; dof = %d ; dof_ent = %d\n",nno,dof_ent);
			  /* jbw */

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal vertex */
#endif



		} // end if active entity

#endif // endif FIELD_BY_FIELD

		  } // end if vertex mixed

		  // if edge
		  else if(l_dof_ent_types[idofent]==PDC_EDGE){

		if(pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
		  mf_fatal_err("In linear approximation should not be EDGE in DOF structure");
		}

		int nne = l_dof_ent_ids[idofent];

		if(apr_get_ent_pdeg(field_id,APC_EDGE,nne)>0){

#ifdef PARALLEL
		  /* in first round only internal edges are considered */
		  if(mmpr_ed_owner(mesh_id,nne)==pcr_my_proc_id()){
#endif

			if(temp_list_edges[nne]==0){

			  temp_list_edges[nne]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_EDGE;
			  (*List_dof_ent_id)[dof_ent]=nne;
			  (*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_1;
			  dof_ent++;
			  nrdofsgl+=mixed_nreq_field_1;

			  /* jbw
			  printf("utr_get_list_int_ent: dof_ent type -> APC_EDGE ; dof = %d ; dof_ent = %d\n",nne,dof_ent);
			  /* jbw */

#ifdef MIXED_DOF_BY_DOF

			  // we add second edge for mixed field p=2
			  if(pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG){

				(*List_dof_ent_type)[dof_ent]=PDC_MIXED_EDGE;
				(*List_dof_ent_id)[dof_ent]=nne;
				(*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_2;
				dof_ent++;
				nrdofsgl+=mixed_nreq_field_2;

			/* jbw
			printf("utr_get_list_int_ent: dof_ent type -> APC_MIXED_EDGE ; dof = %d ; dof_ent = %d\n",nne,dof_ent);
			/* jbw */

				  }

#endif // end if DOF_BY_DOF

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal edge */
#endif
		} // end if active entity
		  } // end if EDGE

		   // if edge mixed
		  else if(l_dof_ent_types[idofent]==PDC_MIXED_EDGE){

#ifdef MIXED_FIELD_BY_FIELD

		if(pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
		  mf_fatal_err("In linear approximation should not be EDGE in DOF structure");
		}

		int nne_mixed_shift = mmr_get_max_edge_id(mesh_id);
		int nne = l_dof_ent_ids[idofent];
		int loc_nne = l_dof_ent_ids[idofent]+nne_mixed_shift;

		if(apr_get_ent_pdeg(second_field_id,APC_MIXED_EDGE,nne)>0){

#ifdef PARALLEL
		  /* in first round only internal edges are considered */
		  if(mmpr_ed_owner(mesh_id,nne)==pcr_my_proc_id()){
#endif

			if(temp_list_edges[loc_nne]==0){

			  temp_list_edges[loc_nne]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_MIXED_EDGE;
			  (*List_dof_ent_id)[dof_ent]=nne;
			  (*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_2;
			  dof_ent++;
			  nrdofsgl+=mixed_nreq_field_2;

			  /*jbw
			  printf("utr_get_list_int_ent: dof_ent type -> APC_MIXED_EDGE ; dof = %d ; dof_ent = %d\n",nne,dof_ent);
			  /*jew*/

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal edge */
#endif
		} // end if active entity

#endif // endif FIELD_BY_FIELD

		  } // end if mixed EDGE

		  // if face
		  else if(l_dof_ent_types[idofent]==PDC_FACE){

		if(pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
		  mf_fatal_err("In linear approximation should not be FACE in DOF structure");
		}

		int nnf = l_dof_ent_ids[idofent];

		if(apr_get_ent_pdeg(field_id,APC_FACE,nnf)>0){

#ifdef PARALLEL
		  /* in first round only internal faces are considered */
		  if(mmpr_ed_owner(mesh_id,nnf)==pcr_my_proc_id()){
#endif

			if(temp_list_faces[nnf]==0){

			  temp_list_faces[nnf]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_FACE;
			  (*List_dof_ent_id)[dof_ent]=nnf;
			  (*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_1;
			  dof_ent++;
			  nrdofsgl+=mixed_nreq_field_1;

			  /* jbw
			  printf("utr_get_list_int_ent: dof_ent type -> APC_FACE ; dof = %d ; dof_ent = %d\n",nnf,dof_ent);
			  /* jbw */


#ifdef MIXED_DOF_BY_DOF

			  // we add second edge for mixed field p=2
			  if(pdeg == APC_MIXED_P2_P2_APPROXIMATION_PDEG){

				(*List_dof_ent_type)[dof_ent]=PDC_MIXED_FACE;
				(*List_dof_ent_id)[dof_ent]=nnf;
				(*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_2;
				dof_ent++;
				nrdofsgl+=mixed_nreq_field_2;

			/* jbw
			printf("utr_get_list_int_ent: dof_ent type -> APC_MIXED_FACE ; dof = %d ; dof_ent = %d\n",nnf,dof_ent);
			/* jbw */

				  }

#endif // end if DOF_BY_DOF

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal face */
#endif
		} // end if active entity
		  } // end if FACE

		  // if face mixed
		  else if(l_dof_ent_types[idofent]==PDC_MIXED_FACE){

#ifdef MIXED_FIELD_BY_FIELD

		if(pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
		  mf_fatal_err("In linear approximation should not be FACE in DOF structure");
		}

		int nnf_mixed_shift = mmr_get_max_face_id(mesh_id);
		int nnf = l_dof_ent_ids[idofent];
		int loc_nnf = l_dof_ent_ids[idofent]+nnf_mixed_shift;

		if(apr_get_ent_pdeg(second_field_id,APC_MIXED_FACE,nnf)>0){

#ifdef PARALLEL
		  /* in first round only internal faces are considered */
		  if(mmpr_ed_owner(mesh_id,nnf)==pcr_my_proc_id()){
#endif

			if(temp_list_faces[loc_nnf]==0){

			  temp_list_faces[loc_nnf]=1;

			  (*List_dof_ent_type)[dof_ent]=PDC_MIXED_FACE;
			  (*List_dof_ent_id)[dof_ent]=nnf;
			  (*List_dof_ent_nrdofs)[dof_ent]=mixed_nreq_field_2;
			  dof_ent++;
			  nrdofsgl+=mixed_nreq_field_2;

			  /* jbw
			  printf("utr_get_list_int_ent: dof_ent type -> APC_MIXED_FACE ; dof = %d ; dof_ent = %d\n",nnf,dof_ent);
			  /* jbw */

			} // if entity not yet considered
#ifdef PARALLEL
		  } /* end if internal face */
#endif
		} // end if active entity

#endif // endif FIELD_BY_FIELD

		  } // end if mixed FACE
		  else {
		mf_fatal_err("UNKNOWN TYPE!!!\n");
		exit(-1);
		  }

		} // end if DOF entity on a list from integration entity
	  } // end if element
	} // end loop over integration entities
	  }
	// END MIXED APPROXIMATION


/*kbw
  printf("After first round pdr_get_list_ent - internal nodes\n");
  printf("\nNr_dof_ent %d, Nrdofs_glob %d, Max_dofs_per_dof_ent %d\n",
  dof_ent, nrdofsgl, nreq);
  for(i=0;i<dof_ent;i++)  printf("type %d, id %d, nrdofs %d\n",
  (*List_dof_ent_type)[i],(*List_dof_ent_id)[i],(*List_dof_ent_nrdofs)[i]);
//kew*/

#ifdef PARALLEL

	printf("Quadratic approximation does not work with parallel solver yet!\n");
	exit(-1);

#endif

	*Nrdofs_glob = nrdofsgl;
	*Max_dofs_per_dof_ent = nreq;
	*Nr_dof_ent = dof_ent;

/*kbw
	  printf("In pdr_get_list_ent for standard quadratic approximation\n");
	  printf("nrelem %d, nrface %d, nr_int_ent %d\n",
	  nr_elem, nr_face, *Nr_int_ent);
	  for(i=0;i<int_ent;i++)  printf("type %d, id %d\n",(*List_int_ent_type)[i],(*List_int_ent_id)[i]);
	  printf("\nNr_dof_ent %d, Nrdofs_glob %d, Max_dofs_per_dof_ent %d\n",*Nr_dof_ent, *Nrdofs_glob, *Max_dofs_per_dof_ent);
	  for(i=0;i<dof_ent;i++)  printf("type %d, id %d, nrdofs %d\n",(*List_dof_ent_type)[i],(*List_dof_ent_id)[i],(*List_dof_ent_nrdofs)[i]);
/*kew*/

	free(temp_list_dofs);
	if(temp_list_edges != NULL) {
	  free(temp_list_edges);
	}
	if(temp_list_faces != NULL) {
	  free(temp_list_faces);
	}
	if(temp_list_elems != NULL) {
	  free(temp_list_elems);
	}

  } // end if standard quadratic approximation

  utv_time.create_dof_ent_list += time_clock()-time_begin;

  return(0);
}


/*------------------------------------------------------------
  utr_get_list_ent - to return the list of integration entities - entities
			 for which stiffness matrices and load vectors are
			 provided by the FEM code to the solver module,
			 and DOF entities - entities with which there are dofs
			 associated by the given approximation
------------------------------------------------------------*/
int utr_get_list_ent( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in:  problem (and solver) identification */
  int* Nr_int_ent,    /* out: number of integration entitites */
  int** List_int_ent_type,/* out: list of types of integration entitites */
  int** List_int_ent_id,  /* out: list of IDs of integration entitites */
  int* Nr_dof_ent,    /* out: number of dof entities (entities with which there
				are dofs associated by the given approximation) */
	/* GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!! */
  int** List_dof_ent_type,/* out: list of types of integration entitites */
  int** List_dof_ent_id,  /* out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /* out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent/* out: maximal number of dofs per dof entity */
  )
{

  char field_module_name[100];
  int solver_id, field_id, mesh_id, nreq;
  int nel, nfa, nno, int_ent, dof_ent, nr_elem, nr_face, nr_node, i, iaux, iedg;
  int nrdofsgl, nrdofent, face_neig[2];
  int nrdofsloc, max_dofs_ent_dof;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* associate the suitable field with the problem and the solver */
  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; field_id = pdr_ctrl_i_params(Problem_id,i);
  i=6; solver_id = pdr_ctrl_i_params(Problem_id,i);
  //mesh_id = apr_get_mesh_id(field_id);
  nreq=apr_get_nreq(field_id);

  // check the name of the field module
  apr_module_introduce(field_module_name);

  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0 ||
	  strncmp(field_module_name,"STANDARD_QUADRATIC",18)== 0 ){

	// default procedure for returning a list of integration entities
	// for standard and quadratic approximations
	utr_get_list_int_ent(Problem_id, Nr_int_ent,
			 List_int_ent_type, List_int_ent_id);

	// default procedure for returning lists of dof entities
	// for standard and quadratic approximations
	// given the list of integration entities as input
	utr_get_list_dof_ent(Problem_id, Nr_int_ent,
			 List_int_ent_type, List_int_ent_id,
			 Nr_dof_ent, List_dof_ent_type, List_dof_ent_id,
			 List_dof_ent_nrdofs, Nrdofs_glob, Max_dofs_per_dof_ent);


  } // end if standard linear and quadratic approximations
  else if( strncmp(field_module_name, "DG_SCALAR_PRISM", 15) == 0){

	// for discontinuous Galerkin approximation

	/* prepare arrays */
	nr_elem = mmr_get_nr_elem(mesh_id);
	nr_face = mmr_get_nr_face(mesh_id);

	*List_int_ent_type = (int *) malloc( (nr_elem+nr_face+1)*sizeof(int) );
	*List_int_ent_id = (int *) malloc( (nr_elem+nr_face+1)*sizeof(int) );
	*List_dof_ent_type = (int *) malloc( (nr_elem+1)*sizeof(int) );
	*List_dof_ent_id = (int *) malloc( (nr_elem+1)*sizeof(int) );
	*List_dof_ent_nrdofs = (int *) malloc( (nr_elem+1)*sizeof(int) );

#ifdef PARALLEL
	/* temporary list of elements */
	int * temp_list_dofs;
	iaux = mmr_get_max_elem_id(mesh_id);
	temp_list_dofs = (int *) malloc( (iaux+1)*sizeof(int) );
	for(i=0;i<=iaux;i++) temp_list_dofs[i]=0;
#endif

	if((*List_dof_ent_nrdofs)==NULL){
	  printf("Not enough space for allocating lists in pdr_get_list_ent! Exiting\n");
	  exit(-1);
	}

	nrdofsgl = 0; max_dofs_ent_dof = 0;

	int_ent=0; dof_ent=0;

	/* lop over elements - integration and dof entities in DG */
	nel=0;
	while((nel=mmr_get_next_elem_all(mesh_id, nel))!=0){

#ifdef PARALLEL
	  /* integration concerns only internal elements */
	  if(mmpr_el_owner(mesh_id,nel)==pcr_my_proc_id()){
#endif

	if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE ) {

	  int_ent++;
	  (*List_int_ent_type)[int_ent-1]=PDC_ELEMENT;
	  (*List_int_ent_id)[int_ent-1]=nel;

#ifdef PARALLEL
	  temp_list_dofs[nel]=1;
#endif

	  /* !!! only for DG - dofs are associated only with elements !!! */
	  dof_ent++;
	  nrdofsloc = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
	  (*List_dof_ent_type)[dof_ent-1]=PDC_ELEMENT;
	  (*List_dof_ent_id)[dof_ent-1]=nel;
	  (*List_dof_ent_nrdofs)[dof_ent-1]=nrdofsloc;

	  if(nrdofsloc>max_dofs_ent_dof) max_dofs_ent_dof = nrdofsloc;
	  nrdofsgl += nrdofsloc;

/*kbw
  printf(" element %d int_ent %d, dof_ent %d, nrdof %d\n",
  nel, int_ent, dof_ent, nrdofsgl);
/*kew*/


	} /* end if element included for integration */

#ifdef PARALLEL
	  } /* end if internal element */
#endif

	} /* end for all elements */


	/* loop over faces - integration entities in DG */
	nfa=0;
	while((nfa=mmr_get_next_face_all(mesh_id, nfa))!=0){

#ifdef PARALLEL
	  /* integration concerns only faces not on the inter-subdomain boundary */
	  if(!mmr_fa_sub_bnd(mesh_id, nfa)){
#endif

	if(mmr_fa_status(mesh_id,nfa)==MMC_ACTIVE ) {

#ifdef PARALLEL
	  /* integration concerns only faces adjacent to internal elements */
	  mmr_fa_neig(mesh_id, nfa, face_neig, NULL, NULL, NULL, NULL, NULL);
	  face_neig[0]=abs(face_neig[0]);
	  face_neig[1]=abs(face_neig[1]);
	  if(mmpr_el_owner(mesh_id,face_neig[0])==pcr_my_proc_id() ||
		 (face_neig[1]!=0 &&
		  mmpr_el_owner(mesh_id,face_neig[1])==pcr_my_proc_id())) {
#endif

		/* add face to the list of integration entities */
		int_ent++;
		(*List_int_ent_type)[int_ent-1]=PDC_FACE;
		(*List_int_ent_id)[int_ent-1]=nfa;

		/*kbw
		  printf("face %d int_ent %d\n", nfa, int_ent);
		  /*kew*/

#ifdef PARALLEL
		/* additionally to DOFs of internal elements also DOFs of neighbours
		   of considered faces must be taken into account */
		/* THESE GHOST ENTITIES HAVE NEGATIVE TYPE !!! */
		/* add first neighbor to the list of DOF entities */
		nel=face_neig[0];
		if(temp_list_dofs[nel]==0){
		  temp_list_dofs[nel]=1;
		  dof_ent++;
		  nrdofsloc = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
		  (*List_dof_ent_type)[dof_ent-1]=-PDC_ELEMENT;
		  (*List_dof_ent_id)[dof_ent-1]=nel;
		  (*List_dof_ent_nrdofs)[dof_ent-1]=nrdofsloc;

		  if(nrdofsloc>max_dofs_ent_dof) max_dofs_ent_dof = nrdofsloc;
		  nrdofsgl += nrdofsloc;

/*kbw
		printf(" element %d int_ent %d, dof_ent %d, nrdof %d\n",
		nel, int_ent, dof_ent, nrdofsgl);
/*kew*/
		}
		/* add second neighbor to the list of DOF entities */
		nel=face_neig[1];
		if(nel!=0 && temp_list_dofs[nel]==0){
		  temp_list_dofs[nel]=1;
		  dof_ent++;
		  nrdofsloc = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
		  /* GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!! */
		  (*List_dof_ent_type)[dof_ent-1] = -PDC_ELEMENT;
		  (*List_dof_ent_id)[dof_ent-1]=nel;
		  (*List_dof_ent_nrdofs)[dof_ent-1]=nrdofsloc;

		  if(nrdofsloc>max_dofs_ent_dof) max_dofs_ent_dof = nrdofsloc;
		  nrdofsgl += nrdofsloc;

/*kbw
		printf(" element %d int_ent %d, dof_ent %d, nrdof %d\n",
		nel, int_ent, dof_ent, nrdofsgl);
/*kew*/
		}

	  } /* end if face adjacent to internal element */
#endif

	}

#ifdef PARALLEL
	  } /* end if face not on the inter-subdomain boundary */
#endif

	} /* end loop over all faces: nfa */

#ifdef PARALLEL
	// for parallel DG we add all elements to dof entities
	/* loop over elements - dof entities in DG */
	nel=0;
	while((nel=mmr_get_next_elem_all(mesh_id, nel))!=0){

	  /* integration concerns only internal elements */
	  if(mmpr_el_owner(mesh_id,nel)!=pcr_my_proc_id()){

	if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE  && temp_list_dofs[nel]==0 ) {


	  temp_list_dofs[nel]=1;

	  /* !!! only for DG - dofs are associated only with elements !!! */
	  dof_ent++;
	  nrdofsloc = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
	  (*List_dof_ent_type)[dof_ent-1]= -PDC_ELEMENT;
	  (*List_dof_ent_id)[dof_ent-1]=nel;
	  (*List_dof_ent_nrdofs)[dof_ent-1]=nrdofsloc;

	  if(nrdofsloc>max_dofs_ent_dof) max_dofs_ent_dof = nrdofsloc;
	  nrdofsgl += nrdofsloc;

/*kbw
  printf(" element %d int_ent %d, dof_ent %d, nrdof %d\n",
  nel, int_ent, dof_ent, nrdofsgl);
/*kew*/


	} /* end if element included for integration */

	  } /* end if not internal element */

	} /* end for all elements */
#endif



	*Nr_int_ent = int_ent;
	*Nr_dof_ent = dof_ent;
	*Nrdofs_glob = nrdofsgl;
	*Max_dofs_per_dof_ent = max_dofs_ent_dof;

#ifdef PARALLEL
/*kbw
  sleep(pcr_my_proc_id());
  printf("Subdomain %d\n", pcr_my_proc_id());
/*kew*/
#endif

/*kbw
  printf("In pdr_get_list_ent for DG approximation\n");
  printf("nrelem %d, nrface %d, nr_int_ent %d\n",
  nr_elem, nr_face, *Nr_int_ent);
  for(i=0;i<int_ent;i++)  printf("type %d, id %d\n",
  (*List_int_ent_type)[i],(*List_int_ent_id)[i]);
  printf("\nNr_dof_ent %d, Nrdofs_glob %d, Max_dofs_per_dof_ent %d\n",
  *Nr_dof_ent, *Nrdofs_glob, *Max_dofs_per_dof_ent);
  for(i=0;i<dof_ent;i++)  printf("type %d, id %d, nrdof %d\n",
  (*List_dof_ent_type)[i],(*List_dof_ent_id)[i],(*List_dof_ent_nrdofs)[i]);
/*kew*/

#ifdef PARALLEL
	free(temp_list_dofs);
#endif

  } // end if DG approximation
  else{


	printf("utr_get_list_ent - unknown approximation type!\n");
	printf("Exiting.\n");
	exit(-1);
  }


/*kbw
	  printf("After pdr_get_list_ent: \n");
	  printf("nr_int_ent %d\n", *Nr_int_ent);
	  for(i=0;i<*Nr_int_ent;i++)
	printf("type %d, id %d\n",(*List_int_ent_type)[i],(*List_int_ent_id)[i]);
	  printf("\nNr_dof_ent %d, Nrdofs_glob %d, Max_dofs_per_dof_ent %d\n",
		 *Nr_dof_ent, *Nrdofs_glob, *Max_dofs_per_dof_ent);
	  for(i=0;i<*Nr_dof_ent;i++)
	printf("type %d, id %d, nrdofs %d\n",
		   (*List_dof_ent_type)[i],(*List_dof_ent_id)[i],(*List_dof_ent_nrdofs)[i]);
/*kew*/

  return(1);
}


/*------------------------------------------------------------
  utr_get_list_ent_coarse - the same as above but for COARSE level and
							given the corresponding lists from the fine level
------------------------------------------------------------*/
int utr_get_list_ent_coarse( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in:  problem (and solver) identification */
  int Nr_int_ent_fine,	/* in: number of integration entitites */
  int *List_int_ent_type_fine,	/* in: list of types of integration entitites */
  int *List_int_ent_id_fine,	/* in: list of IDs of integration entitites */
  int Nr_dof_ent_fine,	/* in: number of dof entities (entities with which there
			   are dofs associated by the given approximation) */
  int *List_dof_ent_type_fine,	/* in: list of types of integration entitites */
  int *List_dof_ent_id_fine,	/* in: list of IDs of integration entitites */
  int *List_dof_ent_nrdof_fine,	/* in: list of no of dofs for 'dof' entity */
  int Nrdof_glob_fine,	/* in: global number of degrees of freedom (unknowns) */
  int Max_dof_per_ent_fine,	/* in: maximal number of dofs per dof entity */
  int *Pdeg_coarse_p,	/* in: degree of approximation for coarse space */
  int *Nr_int_ent_p,	/* out: number of integration entitites */
  int **List_int_ent_type,	/* out: list of types of integration entitites */
  int **List_int_ent_id,	/* out: list of IDs of integration entitites */
  int *Nr_dof_ent_p,	/* out: number of dof entities (entities with which there
			   are dofs associated by the given approximation) */
  int **List_dof_ent_type,	/* out: list of types of integration entitites */
  int **List_dof_ent_id,	/* out: list of IDs of integration entitites */
  int **List_dof_ent_nrdofs,	/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob_p,    /* out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent_p /* out: maximal number of dofs per dof entity */
  )
{


  char field_module_name[100];
  int solver_id, field_id, mesh_id;
  int i, iaux, ient, nel, nfa, nreq;
  int nrdofgl, nr_int_ent, nr_dof_ent;
  int nrdofloc, max_dof_ent_dof;
  int ient_fine, int_ent_id, int_ent_type, dof_ent_id, dof_ent_type;
  int father, father_old, pdeg_coarse;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* associate the suitable field with the solver */
  /* associate the suitable field with the problem and the solver */
  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; field_id = pdr_ctrl_i_params(Problem_id,i);
  i=6; solver_id = pdr_ctrl_i_params(Problem_id,i);

  pdeg_coarse = *Pdeg_coarse_p;
  nreq = apr_get_nreq(field_id);

  // check the name of the field module
  apr_module_introduce(field_module_name);

  // for standard continuous linear FEM approximation
  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0 ||
	  strncmp(field_module_name,"STANDARD_QUADRATIC",18) == 0){

	printf("pdr_get_list_ent_coarse NOT IMPLEMENTED for STD and QUAD approximations!");
	exit (-1);

  }
  else if( strncmp(field_module_name,"DG_SCALAR_PRISM ",15) == 0) { // for DG approximation

/*ok_kbw*/
	printf("\nConstructing coarse mesh from the fine mesh: Problem_id %d\n",
	   Problem_id);
	printf("Nr_int_ent_fine %d, Nr_dof_ent_fine %d, Nrdof_glob_fine %d\n",
	   Nr_int_ent_fine, Nr_dof_ent_fine, Nrdof_glob_fine);
	printf("Max_dof_per_ent_fine %d, pdeg_coarse %d\n",
	   Max_dof_per_ent_fine, pdeg_coarse);
/*kew*/

	int nr_elem = mmr_get_nr_elem(mesh_id);
	int nr_face = mmr_get_nr_face(mesh_id);
	int max_elem_id = mmr_get_max_elem_id(mesh_id);

	// get generation level of first fine element
	int el_gen =0;
	for(ient_fine=0; ient_fine<Nr_int_ent_fine; ient_fine++){

	  int_ent_id = List_int_ent_id_fine[ient_fine];
	  int_ent_type = List_int_ent_type_fine[ient_fine];


	  if(int_ent_type == PDC_ELEMENT){

	el_gen = mmr_el_gen(mesh_id,int_ent_id);
/*kbw
	  printf("In int_ent fine: ient %d, id %d, type %d, el_gen %d\n",
		 ient_fine, int_ent_id, int_ent_type, el_gen);
/*kew*/
	break;

	  }
	}

	int gen_lev = el_gen-1;

	int* temp_list_int_ent_id = (int *) malloc( (Nr_int_ent_fine+1)*sizeof(int) );
	int* temp_list_int_ent_type = (int *) malloc( (Nr_int_ent_fine+1)*sizeof(int) );
	int* temp_list_dof_ent_id = (int *) malloc( (Nr_dof_ent_fine+1)*sizeof(int) );
	int* temp_list_dof_ent_type = (int *) malloc( (Nr_dof_ent_fine+1)*sizeof(int) );
	int* temp_list_dof_ent_nrdofs = (int *) malloc( (Nr_dof_ent_fine+1)*sizeof(int) );
	int* temp_list_bl_for_el = (int *) malloc( (max_elem_id+1)*sizeof(int) );
	for(i=0;i<=max_elem_id;i++) temp_list_bl_for_el[i] = -1;

	int nrdofs_glob = 0;
	int max_dofs_per_dof_ent=0;

	int int_ent=0; int dof_ent=0;
	int nel=0;
	while((nel=mmr_get_next_elem_all(mesh_id, nel))!=0){

#ifdef PARALLEL
	  /* first internal elements are considered */
	  if(mmpr_el_owner(mesh_id,nel)==pcr_my_proc_id()){
#endif

	/* get element generation level with respect to reference*/
	el_gen = mmr_el_gen(mesh_id,nel);

	/* enroll element on a list of mesh elements - no renumbering */
	if((el_gen==gen_lev) ||
	   (el_gen<gen_lev && mmr_el_status(mesh_id,nel)==MMC_ACTIVE )){

#ifdef DEBUG
	  // in the current setting there are no active elements in coarse meshes
	  if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE ){
		printf("Active element %d in coarse mesh! Exiting!\n", nel);
		exit(-1);
	  }
#endif


	  temp_list_dof_ent_id[dof_ent]=nel;
	  temp_list_dof_ent_type[dof_ent]=APC_ELEMENT;
	  temp_list_bl_for_el[nel]=dof_ent;

	  /* integration concerns only internal elements */
	  temp_list_int_ent_id[int_ent]=nel;
	  temp_list_int_ent_type[int_ent]=APC_ELEMENT;

	  int nrdofs_el;
	  if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE){
		nrdofs_el = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
	  }
	  else{
		nrdofs_el = nreq * apr_get_el_pdeg_numshap(field_id, nel, &pdeg_coarse);
	  }
	  nrdofs_glob += nrdofs_el;
	  temp_list_dof_ent_nrdofs[dof_ent]=nrdofs_el;
	  if(nrdofs_el > max_dofs_per_dof_ent) max_dofs_per_dof_ent = nrdofs_el;

/*kbw
	  printf(" element %d int_ent %d, dof_ent %d, nrdof %d\n",
		  nel, int_ent, dof_ent, nrdofs_glob);
/*kew*/

	  dof_ent++;
	  int_ent++;

	} /* end if element included for integration */

#ifdef PARALLEL
	  } /* end if internal element */
#endif

	} /* end for all elements */


	nfa=0;
	while((nfa=mmr_get_next_face_all(mesh_id, nfa))!=0){

#ifdef PARALLEL
/*||begin||*/
	  if(!mmr_fa_sub_bnd(mesh_id, nfa)){
/*||end||*/
#endif

	int face_neig[2];
	mmr_fa_neig(mesh_id, nfa, face_neig, NULL, NULL, NULL, NULL, NULL);

#ifdef DEBUG
	// in the current setting (uniform refinements) all neighbours must be
	// equal size neighbours
	if(face_neig[0]<0 ){
	  printf("Big neighbour element %d in coarse mesh! Exiting!\n",face_neig[0]);
	  exit(-1);
	}
	if(face_neig[1]<0 ){
	  printf("Big neighbour element %d in coarse mesh! Exiting!\n",face_neig[1]);
	  exit(-1);
	}
#endif

	face_neig[0]=abs(face_neig[0]);
	face_neig[1]=abs(face_neig[1]);


	iaux=0;

	/* get neighbors' generation levels with respect to reference*/
	int geneig1=-1;
	int geneig2=-1;

	if(face_neig[0]>0) geneig1 = mmr_el_gen(mesh_id,face_neig[0]);
	if(face_neig[1]>0) geneig2 = mmr_el_gen(mesh_id,face_neig[1]);

	/* check whether neighbors are active for a given mesh */
	if((geneig1==gen_lev) ){
	  // ||  (geneig1<gen_lev && mmr_el_status(mesh_id,face_neig[0])==MMC_ACTIVE )){

	  iaux = -1;

	}
	// if face on the boundary (face_neig[1]==0) or both neighbours have proper gen_lev
	if( face_neig[1]==0 || geneig2==gen_lev ){
	  // || (geneig2<gen_lev && mmr_el_status(mesh_id,face_neig[1])==MMC_ACTIVE )){

	  iaux *= -1;

	}

	if(iaux==1){

#ifdef DEBUG
	  // in the current setting there are no active elements in coarse meshes
	  if(face_neig[0]>0 && mmr_el_status(mesh_id,face_neig[0])==MMC_ACTIVE ){
		printf("Active element %d in coarse mesh! Exiting!\n",face_neig[0]);
		exit(-1);
	  }
	  if(face_neig[1]>0 && mmr_el_status(mesh_id,face_neig[1])==MMC_ACTIVE ){
		printf("Active element %d in coarse mesh! Exiting!\n",face_neig[1]);
		exit(-1);
	  }
#endif


#ifdef PARALLEL
	  /* integration concerns only faces adjacent to internal elements */
	  if(mmpr_el_owner(mesh_id,face_neig[0])==pcr_my_proc_id() ||
		 (face_neig[1]!=0 && mmpr_el_owner(mesh_id,face_neig[1])==pcr_my_proc_id() )) {
#endif

		/* add face to the list of integration entities */
		temp_list_int_ent_id[int_ent]=nfa;
		temp_list_int_ent_type[int_ent]=APC_FACE;

/*kbw
  printf("face %d int_ent %d", nfa, int_ent);
	  printf("   neig1 %d (%d), neig2 %d, (%d)\n",
		 abs(face_neig[0]), temp_list_bl_for_el[abs(face_neig[0])],
		 abs(face_neig[1]), temp_list_bl_for_el[abs(face_neig[1])]);
/*kew*/

	  /* additionally to DOFs of internal elements also DOFs of neighbours
		 of considered faces must be taken into account */
	  /* THESE GHOST ENTITIES HAVE NEGATIVE TYPE !!! */
	  /* add first neighbor to the list of DOF entities */

		/* add first neighbor to the list of DOF entities */
		nel=face_neig[0];
		if(temp_list_bl_for_el[nel]<0){

#ifdef PARALLEL
#ifdef DEBUG
		  // the added neighbour should be alien
		  if(mmpr_el_owner(mesh_id,nel)==pcr_my_proc_id()){
		printf("Owned element %d added as face neighbour! Exiting!\n", nel);
		exit(-1);
		  }
#endif
#endif

		  temp_list_bl_for_el[nel]=dof_ent;
		  temp_list_dof_ent_id[dof_ent]=nel;
		  temp_list_dof_ent_type[dof_ent]=-APC_ELEMENT;

		  int nrdofs_el;
		  if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE){
		nrdofs_el = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
		  }
		  else{
		nrdofs_el = nreq * apr_get_el_pdeg_numshap(field_id, nel, &pdeg_coarse);
		  }
		  nrdofs_glob += nrdofs_el;
		  temp_list_dof_ent_nrdofs[dof_ent]=nrdofs_el;
		  if(nrdofs_el > max_dofs_per_dof_ent) max_dofs_per_dof_ent = nrdofs_el;


/*kbw
  printf("element %d dof_ent %d, nrdof %d\n",
	 nel, dof_ent, nrdofs_el);
/*kew*/

		  dof_ent++;

#ifdef DEBUG
	for(dof_ent_id=0; dof_ent_id < dof_ent; dof_ent_id++){
	  if(abs(temp_list_dof_ent_type[dof_ent_id])!=APC_ELEMENT){
	printf("Dof_ent %d not element in DG! Exiting!\n", dof_ent_id);
	exit(-1);
	  }
	}
#endif



		}

		/* add second neighbor to the list of DOF entities */
		nel=face_neig[1];
		if(nel!=0 && temp_list_bl_for_el[nel]<0){

#ifdef PARALLEL
#ifdef DEBUG
		  // the added neighbour should be alien
		  if(mmpr_el_owner(mesh_id,nel)==pcr_my_proc_id()){
		printf("Owned element %d added as face neighbour! Exiting!\n", nel);
		exit(-1);
		  }
#endif
#endif

		  temp_list_bl_for_el[nel]=dof_ent;
		  temp_list_dof_ent_id[dof_ent]=nel;
		  temp_list_dof_ent_type[dof_ent]=-APC_ELEMENT;

		  int nrdofs_el;
		  if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE){
		nrdofs_el = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
		  }
		  else{
		nrdofs_el = nreq * apr_get_el_pdeg_numshap(field_id, nel, &pdeg_coarse);
		  }
		  nrdofs_glob += nrdofs_el;
		  temp_list_dof_ent_nrdofs[dof_ent]=nrdofs_el;
		  if(nrdofs_el > max_dofs_per_dof_ent) max_dofs_per_dof_ent = nrdofs_el;

/*kbw
		printf(" element %d  dof_ent %d, nrdof %d\n",
			nel, dof_ent, nrdofs_el);
		//getchar();
/*kew*/

		  dof_ent++;

#ifdef DEBUG
	for(dof_ent_id=0; dof_ent_id < dof_ent; dof_ent_id++){
	  if(abs(temp_list_dof_ent_type[dof_ent_id])!=APC_ELEMENT){
	printf("Dof_ent %d not element in DG! Exiting!\n", dof_ent_id);
	exit(-1);
	  }
	}
#endif


		}

		int_ent++;

#ifdef PARALLEL
	  } /* end if face adjacent to internal elements */
#endif

	}

#ifdef PARALLEL
/*||begin||*/
	  } /* end if face not on the inter-subdoamin boundary */
/*||end||*/
#endif
	} /* end loop over all faces: nfa */

#ifdef PARALLEL
	// for parallel DG we add all elements to dof entities
	/* loop over elements - dof entities in DG */
	nel=0;
	while((nel=mmr_get_next_elem_all(mesh_id, nel))!=0){

	  if(mmpr_el_owner(mesh_id,nel)!=pcr_my_proc_id()){

	/* get element generation level with respect to reference*/
	el_gen = mmr_el_gen(mesh_id,nel);

	/* enroll element on a list of mesh elements - no renumbering */
	if((el_gen==gen_lev) && temp_list_bl_for_el[nel] < 0 ){

#ifdef DEBUG
	  // in the current setting there are no active elements in coarse meshes
	  if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE ){
		printf("Active element %d in coarse mesh! Exiting!\n", nel);
		exit(-1);
	  }
#endif



	  temp_list_dof_ent_id[dof_ent]=nel;
	  temp_list_dof_ent_type[dof_ent]= -APC_ELEMENT;
	  temp_list_bl_for_el[nel]=dof_ent;


	  int nrdofs_el;
	  if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE){
		nrdofs_el = apr_get_ent_nrdofs(field_id, APC_ELEMENT, nel);
	  }
	  else{
		nrdofs_el = nreq * apr_get_el_pdeg_numshap(field_id, nel, &pdeg_coarse);
	  }
	  nrdofs_glob += nrdofs_el;
	  temp_list_dof_ent_nrdofs[dof_ent]=nrdofs_el;
	  if(nrdofs_el > max_dofs_per_dof_ent) max_dofs_per_dof_ent = nrdofs_el;

/*kbw
  printf(" element %d int_ent %d, dof_ent %d, nrdof %d\n",
  nel, int_ent, dof_ent, nrdofs_glob);
/*kew*/

	  dof_ent++;



	} /* end if element included for integration */

	  } /* end if not internal element */

	} /* end for all elements */
#endif


	*Nr_int_ent_p = int_ent;
	(*List_int_ent_type) = (int *) malloc((int_ent+1)*sizeof(int));
	(*List_int_ent_id) = (int *) malloc((int_ent+1)*sizeof(int));

	*Nr_dof_ent_p = dof_ent;
	(*List_dof_ent_type) = (int *) malloc((dof_ent+1)*sizeof(int));
	(*List_dof_ent_id) = (int *) malloc((dof_ent+1)*sizeof(int));
	(*List_dof_ent_nrdofs) = (int *) malloc((dof_ent+1)*sizeof(int));

	if((*List_dof_ent_nrdofs)==NULL){
	  printf("Not enough space for allocating lists in pdr_get_list_ent! Exiting\n");
	  exit(-1);
	}

	for(int_ent_id=0; int_ent_id < int_ent; int_ent_id++){
	  (*List_int_ent_type)[int_ent_id]=temp_list_int_ent_type[int_ent_id];
	  (*List_int_ent_id)[int_ent_id]=temp_list_int_ent_id[int_ent_id];
	}

	for(dof_ent_id=0; dof_ent_id < dof_ent; dof_ent_id++){
	  (*List_dof_ent_type)[dof_ent_id]=temp_list_dof_ent_type[dof_ent_id];

#ifdef DEBUG
	  //for(dof_ent_id=0; dof_ent_id < dof_ent; dof_ent_id++){
	  if(abs(temp_list_dof_ent_type[dof_ent_id])!=APC_ELEMENT){
	printf("Dof_ent %d not element in DG! Exiting!\n", dof_ent_id);
	exit(-1);
	  }
	  //}
#endif

	  (*List_dof_ent_id)[dof_ent_id]=temp_list_dof_ent_id[dof_ent_id];
	  (*List_dof_ent_nrdofs)[dof_ent_id]=temp_list_dof_ent_nrdofs[dof_ent_id];
	}

	*Nrdofs_glob_p = nrdofs_glob;
	*Max_dofs_per_dof_ent_p = max_dofs_per_dof_ent;
	*Pdeg_coarse_p = pdeg_coarse;

	free(temp_list_int_ent_id);
	free(temp_list_int_ent_type);
	free(temp_list_dof_ent_id);
	free(temp_list_dof_ent_type);
	free(temp_list_dof_ent_nrdofs);
	free(temp_list_bl_for_el);

  }
  else { // experimental for probable future use

	int* temp_list_elem, * temp_list_face;

	father_old = -1; nr_int_ent = 0;
	temp_list_face = (int *) malloc( (Nr_int_ent_fine+1)*sizeof(int) );
	/* in a loop over fine level integration and dof entities */
	/*    - create lists of integration and dof entities for coarse level */
	/*    - find data related to coarse entities */
	for(ient_fine=0; ient_fine<Nr_int_ent_fine; ient_fine++){

	  int_ent_id = List_int_ent_id_fine[ient_fine];
	  int_ent_type = List_int_ent_type_fine[ient_fine];

	  /*kbw
	printf("In int_ent fine: ient %d, id %d, type %d\n",
	ient_fine, int_ent_id, int_ent_type);
	/*kew*/

	  if(int_ent_type == PDC_ELEMENT){

	father=mmr_el_fam(mesh_id,int_ent_id,NULL,NULL);

	/*kbw
	  printf("father: %d (previous %d)\n", father, father_old);
	  /*kew*/

	if(father==0){

	  temp_list_face[nr_int_ent]=int_ent_id;
	  nr_int_ent++;

	}
	/* we assume families are listed together !!!!*/
	else if(father != father_old){

	  temp_list_face[nr_int_ent]=father;
	  nr_int_ent++;
	  father_old = father;

	}
	  }
	  else if(int_ent_type == PDC_FACE){

	father = mmr_fa_fam(mesh_id, int_ent_id, NULL, NULL);

/*kbw
	  printf("father: %d (previous %d)\n", father, father_old);
/*kew*/

	if(father==0){
	  int face_neig[2];
	  /* check whether neighbors are initial mesh elements */
	  mmr_fa_eq_neig(mesh_id, int_ent_id, face_neig, NULL, NULL);
	  if(mmr_el_fam(mesh_id, face_neig[0], NULL, NULL)==0){

/*kbw
	  printf("face %d, neig %d, father_neig %d\n", int_ent_id,
		 face_neig[0], mmr_el_fam(mesh_id, face_neig[0], NULL, NULL));
/*kew*/

		temp_list_face[nr_int_ent]=int_ent_id;
		nr_int_ent++;
	  }
	}
	else{

	  if(father != father_old){

		temp_list_face[nr_int_ent]=father;
		nr_int_ent++;
		father_old = father;

	  }


	}
	  }
	  else{

	printf("Unknown int_ent type %d in get_list_ent_coarse. Exiting!\n",
		   int_ent_type);
	exit(-1);

	  }
	}

	father_old = -1; nr_dof_ent = 0;
	temp_list_elem = (int *) malloc( (Nr_dof_ent_fine+1)*sizeof(int));
	/* in a loop over fine level integration and dof entities */
	/*    - create lists of integration and dof entities for coarse level */
	/*    - find data related to coarse entities */
	for(ient_fine=0; ient_fine<Nr_dof_ent_fine; ient_fine++){

	  dof_ent_id = List_dof_ent_id_fine[ient_fine];
	  dof_ent_type = List_dof_ent_type_fine[ient_fine];

/*kbw
	printf("In dofent fine: ient %d, id %d, type %d\n",
	   ient_fine, dof_ent_id, dof_ent_type);
/*kew*/

	  if(dof_ent_type == PDC_ELEMENT){

	father=mmr_el_fam(mesh_id,dof_ent_id,NULL,NULL);

/*kbw
	  printf("father: %d (previous %d)\n", father, father_old);
/*kew*/

	if(father==0){

	  temp_list_elem[nr_dof_ent]=dof_ent_id;
	  nr_dof_ent++;

	}
	/* we assume families are listed together !!!!*/
	else if(father != father_old){

	  temp_list_elem[nr_dof_ent]=father;
	  nr_dof_ent++;
	  father_old = father;

	}
	  }
	  else{

	printf("Unknown dof_ent type %d in get_list_ent_coarse. Exiting!\n",
		   dof_ent_type);
	exit(-1);

	  }

	}

/*kbw
  printf("%d coarse elements:\n", nr_dof_ent);
  for(i=0;i<nr_dof_ent; i++){
	printf("%d  ", temp_list_elem[i]);
	if(temp_list_elem[i] != temp_list_face[i]){
	  printf("Not DG - really?\n");
	}
  }
  printf("\n");
  printf("%d coarse faces:\n", nr_int_ent-nr_dof_ent);
  for(i=nr_dof_ent;i<nr_int_ent; i++){
	printf("%d  ", temp_list_face[i]);
  }
  printf("\n");
/*kew*/

	if(pdeg_coarse <=0){
	  iaux = temp_list_elem[0];
	  if(mmr_el_status(mesh_id,iaux)==MMC_ACTIVE){
	pdeg_coarse = apr_get_ent_pdeg(field_id, APC_ELEMENT, iaux);
	  }
	  else{
	int istop=0;
	int sons[10];
	while(istop==0){
	  mmr_el_fam(mesh_id,iaux,sons,NULL);
	  for(i=1;i<=sons[0];i++){
		if(mmr_el_status(mesh_id,sons[i])==MMC_ACTIVE){

		  pdeg_coarse = apr_get_ent_pdeg(field_id, APC_ELEMENT, sons[i]);;

/*kbw
		printf("original element %d, antecedent %d, pdeg %d, nrdofloc %d\n",
		   nel, sons[i], pdeg_coarse, nrdofloc);
/*kew*/

		  istop=1;
		  break;
		}
	  }
	  iaux = sons[1];
	}

	  }
	}


/*kbw
  printf("nr_int_ent %d, nr_dof_ent %d\n", nr_int_ent, nr_dof_ent);
/*kew*/

  /* prepare arrays */

	*List_int_ent_type = (int *) malloc( (nr_int_ent+1)*sizeof(int) );
	*List_int_ent_id = (int *) malloc( (nr_int_ent+1)*sizeof(int) );
	*List_dof_ent_type = (int *) malloc( (nr_dof_ent+1)*sizeof(int) );
	*List_dof_ent_id = (int *) malloc( (nr_dof_ent+1)*sizeof(int) );
	*List_dof_ent_nrdofs = (int *) malloc( (nr_dof_ent+1)*sizeof(int) );

	if((*List_dof_ent_nrdofs)==NULL){
	  printf("Not enough space for allocating lists in pdr_get_list_ent! Exiting\n");
	  exit(-1);
	}



	nrdofgl = 0; max_dof_ent_dof = 0;
	for(ient=0;ient<nr_dof_ent; ient++){

	  nel = temp_list_elem[ient];
	  (*List_int_ent_type)[ient]=PDC_ELEMENT;
	  (*List_int_ent_id)[ient]=nel;

	  nrdofloc = nreq*apr_get_el_pdeg_numshap(field_id, nel, &pdeg_coarse);
	  (*List_dof_ent_type)[ient]=PDC_ELEMENT;
	  (*List_dof_ent_id)[ient]=nel;
	  (*List_dof_ent_nrdofs)[ient]=nrdofloc;

	  if(nrdofloc>max_dof_ent_dof) max_dof_ent_dof = nrdofloc;
	  nrdofgl += nrdofloc;

/*kbw
	  printf(" element %d ient %d, nrdofloc %d, nrdofgl %d\n",
		 nel, ient, nrdofloc, nrdofgl);
/*kew*/


	} /* end for all elements */


	/* loop over faces - integration entities in DG */
	for(ient=nr_dof_ent;ient<nr_int_ent; ient++){
	  nfa = temp_list_face[ient];

	  /* add face to the list of integration entities */
	  (*List_int_ent_type)[ient]=PDC_FACE;
	  (*List_int_ent_id)[ient]=nfa;

/*kbw
	  printf("face %d ient %d\n", nfa, ient);
/*kew*/


	} /* end loop over all faces: nfa */


	*Nr_int_ent_p = nr_int_ent;
	*Nr_dof_ent_p = nr_dof_ent;
	*Nrdofs_glob_p = nrdofgl;
	*Max_dofs_per_dof_ent_p = max_dof_ent_dof;
	*Pdeg_coarse_p = pdeg_coarse;

/*kbw
	printf("In pdr_get_list_ent\n");
	printf("nr_int_ent %d\n",*Nr_int_ent_p);
	for(i=0;i<nr_int_ent;i++)  printf("type %d, id %d\n",
	 (*List_int_ent_type)[i],(*List_int_ent_id)[i]);
	printf("\nNr_dof_ent_p %d, Nrdof_glob %d, Max_dof_per_ent %d\n",
	   *Nr_dof_ent_p, *Nrdof_glob, *Max_dof_per_ent);
	for(i=0;i<nr_dof_ent;i++)  printf("type %d, id %d, nrdof %d\n",
	(*List_dof_ent_type)[i],(*List_dof_ent_id)[i],(*List_dof_ent_nrdof)[i]);

	getchar(); getchar();
/*kew*/

	free(temp_list_elem);
	free(temp_list_face);

  } // end if experimental based on fine level mesh

  return(1);


}

/*------------------------------------------------------------
utr_get_max_num_grid_levels - for limiting nr_levels in multigrid
							  based on mesh and field data
------------------------------------------------------------*/
int utr_get_max_num_grid_levels(
  int Problem_id
)
{

  int gen_max=0;
  char field_module_name[100];

  // currently supported fields:
  // STANDARD_LINEAR - for continuous, vector, linear approximations
  // DG_SCALAR_PRISM - for discontinuous, scalar, high order approximations
  apr_module_introduce(field_module_name);

  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0
	  || strncmp(field_module_name,"STANDARD_QUADRATIC",18) == 0){

	// single level solver for the time being
	gen_max = 0;

  }
  else if( strncmp(field_module_name,"DG_SCALAR_PRISM ",15) == 0) { // for DG approximation

	/* associate the suitable field with the problem and the solver */
	int i=2; int mesh_id=pdr_ctrl_i_params(Problem_id,i);
	i=3; int field_id = pdr_ctrl_i_params(Problem_id,i);
	/* compute the maximal generation level and the minimal p for the mesh */
	int nel=0; int p_min=100;
	while((nel=mmr_get_next_act_elem(mesh_id, nel))!=0){

	  int gen_el = mmr_el_gen(mesh_id,nel);
	  int p_el = apr_get_ent_pdeg(field_id, APC_ELEMENT, nel);

	  if(gen_max<gen_el) gen_max=gen_el;
	  if(p_min>p_el) p_min=p_el;
	}

  }

  return(gen_max+1);
}

/*---------------------------------------------------------
  utr_dof_ent_sons - to return a list of dof entity sons
---------------------------------------------------------*/
int utr_dof_ent_sons( /* returns: info whether the entity is owned */
  int Problem_id,     /* in: problem ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,        /* in: mesh entity ID */
  int *Ent_sons     /* out: list of dof entity sons (for owned entities) */
					 /* 	Ent_sons[0] - number of sons */
  )
{

  int field_id, mesh_id, i;

  int is_owned=1;

  field_id = pdr_ctrl_i_params(Problem_id,3);
  mesh_id = apr_get_mesh_id(field_id);

  if(Ent_type==PDC_ELEMENT){

/*kbw
	printf("utr_dof_ent_sons: element %d, mesh_id %d\n", Ent_id, mesh_id);
/*kew*/

#ifdef PARALLEL
	if(mmpr_el_owner(mesh_id,Ent_id)==pcr_my_proc_id()){
#endif

	  mmr_el_fam(mesh_id, Ent_id, Ent_sons, NULL);

/*kbw
	  printf("sons:");
	  for(i=0;i<Ent_sons[0];i++){
	printf("%d  ", Ent_sons[i+1]);
	  }
	  printf("\n");
/*kew*/

#ifdef PARALLEL
	} else if(mmpr_el_owner(mesh_id,Ent_id)!=pcr_my_proc_id()){
	  is_owned = 0;
	}
#endif

  }
  else if(Ent_type == PDC_FACE){

#ifdef PARALLEL
	if(mmpr_fa_owner(mesh_id,Ent_id)==pcr_my_proc_id()){
#endif

	  mmr_fa_fam(mesh_id, Ent_id, Ent_sons, NULL);

#ifdef PARALLEL
	} else if(mmpr_fa_owner(mesh_id,Ent_id)!=pcr_my_proc_id()){
	  is_owned = 0;
	}
#endif

  }
  else if(Ent_type == PDC_EDGE){

#ifdef PARALLEL
	if(mmpr_ed_owner(mesh_id,Ent_id)==pcr_my_proc_id()){
#endif

	  mmr_edge_sons(mesh_id, Ent_id, Ent_sons, NULL);

#ifdef PARALLEL
	} else if(mmpr_ed_owner(mesh_id,Ent_id)!=pcr_my_proc_id()){
	  is_owned = 0;
	}
#endif

  } else { // for vertices

#ifdef PARALLEL
	if(mmpr_ve_owner(mesh_id,Ent_id)==pcr_my_proc_id()){
#endif

	  Ent_sons[0] = 0;

#ifdef PARALLEL
	} else if(mmpr_ve_owner(mesh_id,Ent_id)!=pcr_my_proc_id()){
	  is_owned = 0;
	}
#endif

  }

  return(is_owned);
}


/*------------------------------------------------------------
  utr_comp_stiff_mat - to create a stiffness matrix
					  and a load vector corresponding to the specified
					  mesh entity, together with information on how to
					  assemble entries into the global stiffness matrix
					  and the global load vector
------------------------------------------------------------*/
int utr_comp_stiff_mat(	 /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in: approximation field ID  */
  int Int_ent_type,	/* in: unique identifier of the integration entity */
  int Int_ent_id,	/* in: unique identifier of the integration entity */
  int Comp_sm,	/* in: indicator for the scope of computations: */
  /*   PDC_NO_COMP  - do not compute anything */
  /*   PDC_COMP_SM - compute entries to stiff matrix only */
  /*   PDC_COMP_RHS - compute entries to rhs vector only */
  /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,	/* in: enforced degree of polynomial (if > 0 ) */
  int *Nr_dof_ent,	/* in: size of arrays, */
			/* out: number of mesh entities with which dofs and */
			/*      stiffness matrix blocks are associated */
  int *List_dof_ent_type,	/* out: list of types for 'dof' entities */
  int *List_dof_ent_id,	/* out: list of ids for 'dof' entities */
  int *List_dof_ent_nrdofs,	/* out: list of no of dofs for 'dof' entity */
  int *Nrdofs_loc,	/* in(optional): size of Stiff_mat and Rhs_vect */
		   /* out(optional): actual number of dofs per integration entity */
  double *Stiff_mat,	/* out(optional): stiffness matrix stored columnwise */
  double *Rhs_vect,	/* outpds_elast_ls_std_intf.c(optional): rhs vector */
  char *Rewr_dofs	/* out(optional): flag to rewrite or sum up entries */
			/*   'T' - true, rewrite entries when assembling */
			/*   'F' - false, sum up entries when assembling */
 )
{

  if(Comp_sm != PDC_NO_COMP){

	/*jbw
	printf("Comp_sm != PDC_NO_COMP\n"); getchar();
	/*jbw*/

	double time_begin = time_clock();

	// first construct unconstrained stiffness matrix
	pdr_comp_stiff_mat_uncon(Problem_id, Int_ent_type, Int_ent_id, Comp_sm,
				 Pdeg_vec, Nrdofs_loc,
				 Stiff_mat, Rhs_vect, Rewr_dofs);

	utv_time.create_sm_uncon += time_clock()-time_begin;

  }

  double time_begin = time_clock();

  // fill lists of DOF entities (taking constraints into account)

  /* change the option compute SM and RHSV to rewrite SM and RHSV */
  if (Comp_sm != PDC_NO_COMP)  Comp_sm += 3;

  int field_id = pdr_ctrl_i_params(Problem_id, 3);

  int elem, mesh_id, face_neig[2];
  if (Int_ent_type == PDC_ELEMENT) {

	elem = Int_ent_id;

  } else if (Int_ent_type == PDC_FACE) {

	mesh_id = apr_get_mesh_id(field_id);
	mmr_fa_neig(mesh_id, Int_ent_id, face_neig, NULL, NULL, NULL, NULL, NULL);

	// we get the first neighbour, since stiffness matrices are created only
	// for boundary faces - and these faces has only one neighbour, always the first
	elem = abs(face_neig[0]);

  }

  int pdeg;
  /* element degree of approximation for linear prisms is a single number */
  if (Pdeg_vec == NULL)
	pdeg = 0;
  else
	pdeg = Pdeg_vec[0];

/*kbw
	 printf("in utr_comp_stiff_mat: problem_id %d, int_ent: type %d, id %d, enforced pdeg %d\n",
	 Problem_id, Int_ent_type, Int_ent_id, pdeg);
/*kew */

	/*jbw
	printf("-> apr_get_stiff_mat_data\n");
	/*jbw*/


  /* obligatory procedure to fill Lists of dof_ents and rewite SM and RHSV */
  /* the reason is to take into account POSSIBLE CONSTRAINTS (HANGING NODES) */
  apr_get_stiff_mat_data(field_id, elem, Comp_sm, 'N',
			 pdeg, 0, Nr_dof_ent, List_dof_ent_type,
			 List_dof_ent_id, List_dof_ent_nrdofs,
			 Nrdofs_loc, Stiff_mat, Rhs_vect);

/*kbw
  if(Comp_sm!=PDC_NO_COMP)
	{
	  int i;
	  printf("In utr_comp_el_stiff_mat: field_id %d, El_id %d, Comp_sm %d, Nr_dof_ent %d\n",
		 field_id, elem, Comp_sm, *Nr_dof_ent);
	  printf("For each block: \ttype, \tid, \tnrdof\n");
	  for(i=0;i<*Nr_dof_ent;i++){
	printf("\t\t\t%d\t%d\t%d\n",
		   List_dof_ent_type[i],List_dof_ent_id[i],List_dof_ent_nrdofs[i]);
	  }
	  printf("\n\n");
	}
  //getchar();getchar();
/*kew */

  /* matrix displayed by rows, altghough stored by columns !!!!!!!!! */
/*kbw
  if(Comp_sm!=PDC_NO_COMP)
	{
	  int idofs, jdofs;
	  printf("\nElement %d: Modified stiffness matrix:\n",elem);
	  for (idofs=0;idofs<*Nrdofs_loc;idofs++) { //for each row!!!!
	for (jdofs=0;jdofs<*Nrdofs_loc;jdofs++) { // for each element in row !!!
	  printf("%7.3lf",Stiff_mat[idofs+jdofs*(*Nrdofs_loc)]);
	}
	printf("\n");
	  }
	  printf("Element %d: Rhs_vect:\n",elem);
	  for (idofs=0;idofs<*Nrdofs_loc;idofs++) {
	printf("%7.3lf",Rhs_vect[idofs]);
	  }
	  printf("\n\n");
	  //getchar();
	}
  /* */

  utv_time.adding_constraints_info += time_clock()-time_begin;

  return (1);
}




/*------------------------------------------------------------
 utr_create_assemble_stiff_mat - to create element stiffness matrices
								 and assemble them to the global SM
------------------------------------------------------------*/
int utr_create_assemble_stiff_mat(
  int Problem_id,
  int Level_id,
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  int Nr_int_ent,
  int* L_int_ent_type,
  int* L_int_ent_id,
  int Nr_colors_elems,
  int* L_color_index_elems,
  int Nr_colors_faces,
  int* L_color_index_faces,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int* Pos_first_dof_int_ent,
  int* Local_to_global,
  int Max_dofs_int_ent
)
{

int i,j,k;


/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kbw
  printf("In utr_create_assemble_stiff_mat: Problem_id %d, Level_id %d\n",
	 Problem_id, Level_id);
  printf("Comp_type %d, Nr_int_ent %d, Max_dofs_int_ent %d\n",
	 Comp_type, Nr_int_ent, Max_dofs_int_ent);
/*kew*/

#define TIME_TEST_INTEGRATION_ASSEMBLY
// we define certain time measurment variables but use them only if switched on above
  double time_integration=0.0;
  double time_assembly=0.0;
  double time_elements=0.0;
  double time_faces=0.0;
  double time_accelerator=0.0;
  double time_concurrent_CPU=0.0;
  double time_total=0.0;
  double time_begin=time_clock(); // we always perform some time measurments

  // divide list into elements and faces
  int nr_elems=0;
  int int_entity;
  int previous_type = PDC_ELEMENT;
  int ok = 1;
  for(int_entity=0;int_entity<Nr_int_ent;int_entity++){
	if(L_int_ent_type[int_entity]==PDC_ELEMENT){
	  nr_elems++;
	  if(previous_type != PDC_ELEMENT) ok = 0;
	}
	else previous_type = PDC_FACE;
  }
  int nr_faces =  Nr_int_ent - nr_elems;

  // divide elements into processed by accelerator and standard routines
  int nr_colors_accel=0;
#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)
#ifdef GPU_SMALL_EXAMPLE
  int nr_elems_accel = nr_elems;///2; // limit of elements for accelerator
#else
  int nr_elems_accel = 0.99*nr_elems; // limit of elements for accelerator
#endif
  for(i=0; i<Nr_colors_elems; i++){ // for all colors
	if(L_color_index_elems[i+1]>nr_elems_accel){ // if color ends after the limit
	  nr_colors_accel = i; // set the number of colors till the limit and break
	  break;
	}
  }
#ifdef GPU_SMALL_EXAMPLE
  nr_colors_accel = Nr_colors_elems;

#else
  if(Nr_colors_elems==1) nr_colors_accel=0; // make it explicit - for one color NO GPU computing
#endif
#endif

  /*ok_kbw*/
  if(nr_colors_accel>0){
	printf("\nStarting parallel OpenMP/OpenCL integration:\n");
  }
  else{
	printf("\nStarting parallel OpenMP only integration:\n");
  }
  printf("nr_elems %d (Nr_colors_elems %d, nr_colors for accelerator %d)\n",
	 nr_elems, Nr_colors_elems, nr_colors_accel);
  if(nr_colors_accel>0){
	printf("nr_elems per color (accel): ");
	for(i=0; i<nr_colors_accel; i++){
	  printf("%d  ",  L_color_index_elems[i+1] - L_color_index_elems[i]);
	}
	printf("\n");
  }
  printf("nr_elems per color (CPU): ");
  for(i=nr_colors_accel; i<Nr_colors_elems; i++){
	printf("%d  ",  L_color_index_elems[i+1] - L_color_index_elems[i]);
  }
  printf("\n");
  printf("nr_faces %d (Nr_colors_faces %d), ok %d\n",
	 nr_faces, Nr_colors_faces, ok);
  printf("nr_faces per color: ");
  for(i=0; i<Nr_colors_faces; i++){
	printf("%d  ",  L_color_index_faces[i+1] - L_color_index_faces[i]);
  }
  printf("\n");
  /*kew*/


  if(ok!=1){
	printf("Elements are not first on the list to integrate. \n");
#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)
	printf("Accelerators do not work in this case.\n");
	nr_colors_accel = 0;
#else
	printf("OpenMP is not optimized for that case.\n");
#endif
	exit(-1);
  }

  int number_of_threads; // number of threads performing CPU integration of elements

  // for OpenMP+OpenCL we create the first parallel region with two threads only
  // the first thread will create
#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)

  omp_set_nested(1);
#pragma omp parallel num_threads(2) if(Nr_colors_faces>1 && Nr_colors_elems>1)  default(none) \
  shared(nr_colors_accel, nr_elems, nr_faces, number_of_threads)	\
  shared(time_elements, time_faces, time_accelerator, time_concurrent_CPU, \
	 time_integration, time_assembly, utv_time)			\
  firstprivate(Nr_colors_elems, Nr_colors_faces, L_color_index_elems, L_color_index_faces) \
  firstprivate(Problem_id, Level_id, Comp_type, Nr_int_ent, Pdeg_coarse_p, \
		   L_int_ent_type, L_int_ent_id,  \
		   Asse_pos_first_dof_int_ent, Assembly_table, \
			   Pos_first_dof_int_ent, Local_to_global, \
			   Max_dofs_int_ent)
  {

#else // for OpenMP only

#pragma omp parallel if(Nr_colors_faces>1 && Nr_colors_elems>1)  default(none) \
  shared(nr_colors_accel, nr_elems, nr_faces, number_of_threads, PDC_COMP_BOTH, PDC_COMP_RHS, PDC_COMP_SM )	\
  firstprivate(Nr_colors_elems, Nr_colors_faces, L_color_index_elems, L_color_index_faces) \
  firstprivate(Problem_id, Level_id, Comp_type, Nr_int_ent, Pdeg_coarse_p, \
		   L_int_ent_type, L_int_ent_id,  \
		   Asse_pos_first_dof_int_ent, Assembly_table, \
			   Pos_first_dof_int_ent, Local_to_global, \
			   Max_dofs_int_ent)			\
  shared(time_elements, time_faces, time_accelerator, \
	time_integration, time_assembly, utv_time)
  {

#endif

	// after if (OpenMP+OpenCL) else (openMP only) endif
	// we are inside parallel region with:
	// OpenMP+OpenCL - two threads
	// openMP - a set of threads (according to OMP_NUM_THREADS values)


#ifdef _OPENMP
	int my_thread_id = omp_get_thread_num();
#else
	int my_thread_id = 0;
#endif

	// appeared after merging dev_forming_merge - to be deleted?
/* <<<<<<< HEAD */
/* ======= */

/*     if(my_thread_id==0) number_of_threads = num_threads; */

/*     int nrdofs_glob, max_nrdofs, nr_dof_ent, nr_levels, posglob, nrdofs_int_ent; */
/*     int l_dof_ent_id[PDC_MAX_DOF_PER_INT], l_dof_ent_nrdof[PDC_MAX_DOF_PER_INT]; */
/*     int l_dof_ent_posglob[PDC_MAX_DOF_PER_INT]; */
/*     int l_dof_ent_type[PDC_MAX_DOF_PER_INT]; */
/*     int i,j,k, iaux, kaux, int_ent, ibl, ient, ini_zero; */
/*     char rewrite; */
/*     int icolor; */

/*     /\* allocate memory for the stiffness matrices and RHS *\/ */
/*     max_nrdofs = Max_dofs_int_ent; */
/*     double *stiff_mat; */
/*     double *rhs_vect; */
/*     if(Comp_type == PDC_COMP_BOTH) { */
/*       stiff_mat = (double *)malloc(max_nrdofs*max_nrdofs*sizeof(double)); */
/*       rhs_vect = (double *)malloc(max_nrdofs*sizeof(double)); */
/*       memset(stiff_mat, 0, max_nrdofs*max_nrdofs*sizeof(double)); */
/*       memset(rhs_vect, 0, max_nrdofs*sizeof(double)); */
/*     } */
/*     else if(Comp_type == PDC_COMP_RHS) { */
/*       rhs_vect = (double *)malloc(max_nrdofs*sizeof(double)); */
/*       memset(rhs_vect, 0, max_nrdofs*sizeof(double)); */
/*     } */
/*     else if(Comp_type == PDC_COMP_SM) { */
/*       stiff_mat = (double *)malloc(max_nrdofs*max_nrdofs*sizeof(double)); */
/*       memset(stiff_mat, 0, max_nrdofs*max_nrdofs*sizeof(double)); */
/*     } */
/*     else{ */
/*       //printf("How the hell you want to solve something without SM and RHS?\n"); */
/*       exit(-1); */
/*     } */

/* /\*kbw */
/* //#pragma omp critical(printing) */
/* { */
/*   printf("In utr_create_assemble_stiff_mat before loop over entities\n"); */
/*   printf("thread_id %d, num_threads %d\n", */
/* 	 omp_get_thread_num(),	 omp_get_num_threads() ); */
/* } */
/* /\*kew*\/ */

/* // BECAUSE WE INTRODUCE COLORS, WE NEED TO SYNCHRONISE ASSEMBLY, HENCE: */
/* // ALL THREADS COMPUTE FACE STIFFNESS MATRICES AND NEXT ALL THREADS COMPUTE */
/* // ELEMENT STIFFNESS MATRICES (WITH SEVERAL COLORS POSSIBLY OFFLOADED TO ACCELERATOR) */
/*     //if( nr_colors_accel==0 || num_threads==1 || (num_threads>1 && my_thread_id!=0) ){ */
/* >>>>>>> dev_forming_merge */



	// WE DEVELOP A GENERIC INTEGRATION ROUTINE FOR DIFFERENT OPENCL DEVICES
	// DIFFERENT VERSIONS ARE CREATED USING SWITCHES (SEE BELOW)
	// PARTICULAR VERSIONS _GPU, _CPU, _PHI - are considered obsolete, but may
	// be used for testing !!!!

#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)

	// the master thread dispatches several colors to accelerator (if present)
	if(my_thread_id==0 && nr_colors_accel>0){
/*kbw
	  printf("In utr_create_assemble_stiff_mat before accel in thread %d\n", my_thread_id);
/*kew*/

	  double t_temp = time_clock();

	  utr_create_assemble_stiff_mat_elem_accel(Problem_id, Level_id, Comp_type,
						   Pdeg_coarse_p,
						   Nr_int_ent, L_int_ent_type, L_int_ent_id,
						   Nr_colors_elems, L_color_index_elems,
						   nr_colors_accel,
						   Asse_pos_first_dof_int_ent,
						   Assembly_table,
						   Pos_first_dof_int_ent,
						   Local_to_global,
						   Max_dofs_int_ent);

	  time_accelerator += time_clock() - t_temp;

	} // end if master thread (my_id==0) (controlling OpenCL create-assemble)

	  // if the second thread for CPU+GPU
	else{


	  double t_temp = time_clock();

	  // THE SECOND PARALLEL REGION!!!
#pragma omp parallel num_threads(omp_get_max_threads()-1) default(none) \
  shared(nr_colors_accel, nr_elems, nr_faces, number_of_threads )	\
  firstprivate( L_color_index_elems, L_color_index_faces)		\
  firstprivate(Problem_id, Level_id, Comp_type, Nr_int_ent, Pdeg_coarse_p, \
		   L_int_ent_type, L_int_ent_id,  \
		   Asse_pos_first_dof_int_ent, Assembly_table, \
			   Pos_first_dof_int_ent, Local_to_global, \
			   Max_dofs_int_ent) \
  shared( Nr_colors_elems, Nr_colors_faces)				\
  shared(time_integration,time_assembly,utv_time)
	  {

	  // end if OpenMP+OpenCL

#else // if OpenMP only

  // we need to start two regions:
	{ // the first corresponding to the "else" in (my_id==0) i.e. to the second
	  // thread in the first parallel region for OpenMP+OpenCL

	  { // the second corresponding to the opening bracket in the second
	// parallel region for OpenMP+OpenCL

#endif // end if OpenMP only

	// all the rest below for both cases OpenMP only and OpenMP+OpenCL
	// but for OpenMP only inside a single parallel region
	//          with two artificial {{ added to match the OpenMP+OpenCL case
	// and for OpenMP+OpenCL additionaly inside the second (nested) parallel region
	//          and the body of "else" for if(my_id==0) in the first region


	// for OpenMP only threads are numbered from 0 to  OMP_NUM_THREADS
	// for OpenMP+OpenCL  threads are numbered from 0 to  OMP_NUM_THREADS-1
	// ( the first thread in the first region has id==0 the second has id==1
	//   but the thread with id==1 starts the second region and gets id==0 inside
	//   hence we have two threads with id==0 .... )
#ifdef _OPENMP
	int my_thread_id = omp_get_thread_num();
	int num_threads = omp_get_num_threads();
#else
	int my_thread_id = 0;
	int num_threads = 1;
#endif

	if(my_thread_id==0) number_of_threads = num_threads;


/*kbw
#pragma omp critical(printing)
	{
	  printf("In utr_create_assemble_stiff_mat before OpenMP loop over entities\n");
	  printf("thread_id %d, num_threads %d\n",
		 my_thread_id,	 number_of_threads );
	}
/*kew*/


// all variables below are private to threads - we are inside parallel region
	int nrdofs_glob, max_nrdofs, nr_dof_ent, nr_levels, posglob, nrdofs_int_ent;
	int l_dof_ent_id[PDC_MAX_DOF_PER_INT], l_dof_ent_nrdof[PDC_MAX_DOF_PER_INT];
	int l_dof_ent_posglob[PDC_MAX_DOF_PER_INT];
	int l_dof_ent_type[PDC_MAX_DOF_PER_INT];
	int i,j,k, iaux, kaux, int_ent, ibl, ient, ini_zero;
	char rewrite;
	int icolor;

	/* allocate memory for the stiffness matrices and RHS */
	max_nrdofs = Max_dofs_int_ent;
	double *stiff_mat;
	double *rhs_vect;
	if(Comp_type == PDC_COMP_BOTH) {
	  stiff_mat = (double *)malloc(max_nrdofs*max_nrdofs*sizeof(double));
	  rhs_vect = (double *)malloc(max_nrdofs*sizeof(double));
	  memset(stiff_mat, 0, max_nrdofs*max_nrdofs*sizeof(double));
	  memset(rhs_vect, 0, max_nrdofs*sizeof(double));
	}
	else if(Comp_type == PDC_COMP_RHS) {
	  rhs_vect = (double *)malloc(max_nrdofs*sizeof(double));
	  memset(rhs_vect, 0, max_nrdofs*sizeof(double));
	}
	else if(Comp_type == PDC_COMP_SM) {
	  stiff_mat = (double *)malloc(max_nrdofs*max_nrdofs*sizeof(double));
	  memset(stiff_mat, 0, max_nrdofs*max_nrdofs*sizeof(double));
	}
	else{
	  printf("How the hell you want to solve something without SM and RHS?\n");
	  exit(-1);
	}


	for(icolor = 0; icolor<Nr_colors_faces; icolor++){


	  int nr_faces_this_color;

	  int flag = 0;
	  int faces_per_thread;
	  int my_first_face;
	  int my_last_face;



#ifdef TIME_TEST_INTEGRATION_ASSEMBLY
  #pragma omp for reduction(+:time_integration) reduction(+:time_assembly)
#else
  #pragma omp for
#endif
	  for(int_ent=L_color_index_faces[icolor];
		  int_ent<L_color_index_faces[icolor+1];
		  int_ent++){

		int nr_dof_ent = PDC_MAX_DOF_PER_INT;
		int nrdofs_int_ent = max_nrdofs;


		/* jbw
		printf("THREAD %d:\tint_ent = %d ; L_color_index_faces[%d] = %d\n",omp_get_thread_num(),int_ent,(icolor+1),L_color_index_faces[icolor+1]);
		/* jbw */

		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// ****************************************************************
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//#pragma omp critical(pdr_comp_stiff_mat)
		{
		  double time_tmp = time_clock();

		  pdr_comp_stiff_mat(Problem_id, L_int_ent_type[int_ent],
							 L_int_ent_id[int_ent], Comp_type,  Pdeg_coarse_p,
							 &nr_dof_ent,l_dof_ent_type,l_dof_ent_id,l_dof_ent_nrdof,
							 &nrdofs_int_ent, stiff_mat, rhs_vect, &rewrite);

/*kbw
#pragma omp critical(printing)
  {
	printf("After pdr_comp_stiff_mat in utr_create_assemble: type %d, id %d, nr_dof_ent %d\n",
	   L_int_ent_type[int_ent], L_int_ent_id[int_ent], nr_dof_ent);
	int idofent;
	for(idofent=0; idofent<nr_dof_ent; idofent++){
	  printf("dof_ent: id %d, type %d, nrdofs %d\n",
		 l_dof_ent_id[idofent],l_dof_ent_type[idofent],l_dof_ent_nrdof[idofent]);
	}
  }
/*kew*/

#pragma omp atomic
		  utv_time.CPU_integration_faces += time_clock() - time_tmp;

#ifdef TIME_TEST_INTEGRATION_ASSEMBLY
#pragma omp atomic
		  time_integration += time_clock() - time_tmp;
#endif
		}


#ifdef DEBUG_SIM
		if(nrdofs_int_ent>max_nrdofs){
		  printf("Too small arrays stiff_mat and rhs_vect passed to comp_el_stiff_mat\n");
		  printf("from sir_create in sis_mkb. %d < %d. Exiting !!!",
			   max_nrdofs, nrdofs_int_ent);
		  exit(-1);
		}
#endif


/*kbw
#pragma omp critical(printing)
	{
	  if(L_int_ent_type[int_ent]==PDC_ELEMENT) {
		//if(L_int_ent_id[int_ent]!=-1) {
		int i=6;
		int solver_id = pdr_ctrl_i_params(Problem_id,i);
		printf("utr_create_assemble: before assemble:\n");
		printf("Problem_id %d, Solver_id %d, Level_id %d, sol_typ %d\n",
		   Problem_id, solver_id, Level_id, Comp_type);
		int ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
		printf("ient %d, int_ent_id %d, int_ent_type %d, nr_dof_ent %d\n",
		   int_ent, L_int_ent_id[int_ent], L_int_ent_type[int_ent],
		   nr_dof_ent);
		if(Comp_type!=PDC_NO_COMP){
		  pli = 0; nrdof=0;
		  for(ibl=0;ibl<nr_dof_ent; ibl++) nrdof+=l_dof_ent_nrdof[ibl];
		  for(ibl=0;ibl<nr_dof_ent; ibl++){
		printf("bl_id %d, bl_nrdof %d\n",
			   l_dof_ent_id[ibl],l_dof_ent_nrdof[ibl]);
		nri = l_dof_ent_nrdof[ibl];
		plj=0;
		for(jbl=0;jbl<nr_dof_ent;jbl++){
		  //printf("Stiff_mat transposed! (blocks %d:%d)\n",jbl,ibl);
		  //	nrj = l_dof_ent_nrdof[jbl];
		  //	for(i=0;i<nri;i++){
		  //	  jaux = plj+(pli+i)*nrdof;
		  //	  for(j=0;j<nrj;j++){
		  //	    printf("%20.15lf",stiff_mat[jaux+j]);
		  //
		  //	    if(stiff_mat[jaux+j]< -1.e20){
		  //	      getchar(); getchar(); getchar();
		  //	    }
		  //
		  ///	  }
		  //	  printf("\n");
		  //	}
		  printf("Stiff_mat (blocks %d:%d)\n",jbl,ibl);
		  nrj = l_dof_ent_nrdof[jbl];
		  for(j=0;j<nrj;j++){
			for(i=0;i<nri;i++){
			  jaux = plj+(pli+i)*nrdof;
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
	  }
	}
	/*kew*/

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ****************************************************************
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		double time_tmp = time_clock();

		if(Nr_colors_faces <= 1){

#ifdef _OPENMP
		  if(omp_in_parallel()){
				printf("In OpenMP parallel region in utr_create_assemble for 1 color!\n. Exiting.\n");
			exit(0);
		  }
#endif

		  // critical section switched off since for one color we serialize the execution
		  //#pragma omp critical(assembling)
		  {

			if(Assembly_table != NULL){
			//if(level_p->l_int_ent_type[intent]==3){

			  int position = Asse_pos_first_dof_int_ent[int_ent];
		  int pos_global = Pos_first_dof_int_ent[int_ent];
			  pdr_assemble_local_stiff_mat_with_table(Problem_id, Level_id, Comp_type,
							nr_dof_ent,
							&Assembly_table[position],
							&Local_to_global[pos_global],
							stiff_mat, rhs_vect, &rewrite);


		}
			else{
			  pdr_assemble_local_stiff_mat(Problem_id, Level_id, Comp_type,
						 nr_dof_ent, l_dof_ent_type,
						 l_dof_ent_id,l_dof_ent_nrdof,
						 stiff_mat, rhs_vect, &rewrite);
			}
		  }

		}
		else{ // if more colors - no critical section !!!

//#pragma omp critical(assembling) //if something goes wrong
		  {
			if(Assembly_table != NULL){

			  //printf("OK - faces with table\n");

			  //if(level_p->l_int_ent_type[intent]==3){
			  int position = Asse_pos_first_dof_int_ent[int_ent];
		  int pos_global = Pos_first_dof_int_ent[int_ent];
			  pdr_assemble_local_stiff_mat_with_table(Problem_id, Level_id, Comp_type,
							nr_dof_ent,
							&Assembly_table[position],
							&Local_to_global[pos_global],
							stiff_mat, rhs_vect, &rewrite);
			}
			else{
			  pdr_assemble_local_stiff_mat(Problem_id, Level_id, Comp_type,
						 nr_dof_ent, l_dof_ent_type,
						 l_dof_ent_id,l_dof_ent_nrdof,
						 stiff_mat, rhs_vect, &rewrite);
			}
		  } // end critical

		}  // end if more colors

#pragma omp atomic
		utv_time.CPU_assembly_faces += time_clock() - time_tmp;

#ifdef TIME_TEST_INTEGRATION_ASSEMBLY
#pragma omp atomic
		time_assembly += time_clock() - time_tmp;
#endif

		  } /* end PARALLEL (OpenMP) loop over single color faces - implicit barrier */

#pragma omp barrier

		} // end loop over colors


		for(icolor = nr_colors_accel; icolor<Nr_colors_elems; icolor++){

	  // just parallel loop over elements of a single color - all threads compute
#ifdef TIME_TEST_INTEGRATION_ASSEMBLY
  #pragma omp for reduction(+:time_integration) reduction(+:time_assembly)
#else
  #pragma omp for
#endif
	  for(int_ent=L_color_index_elems[icolor];
		  int_ent<L_color_index_elems[icolor+1];
		  int_ent++)  {
		 /* jbw
		 printf("THREAD %d:\tint_ent = %d ; L_color_index_elems[%d] = %d\n",omp_get_thread_num(),int_ent,(icolor+1),L_color_index_elems[icolor+1]);
		 /* jbw */

		int nr_dof_ent = PDC_MAX_DOF_PER_INT;
		int nrdofs_int_ent = max_nrdofs;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ****************************************************************
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//#pragma omp critical (pdr_comp_stiff_mat)
		  {

			double time_tmp = time_clock();

			pdr_comp_stiff_mat(Problem_id, L_int_ent_type[int_ent],
						L_int_ent_id[int_ent], Comp_type, Pdeg_coarse_p,
					&nr_dof_ent,l_dof_ent_type,l_dof_ent_id,l_dof_ent_nrdof,
					&nrdofs_int_ent, stiff_mat, rhs_vect, &rewrite);


/*kbw
#pragma omp critical(printing)
  {
	printf("After pdr_comp_stiff_mat in utr_create_assemble: type %d, id %d, nr_dof_ent %d\n",
	   L_int_ent_type[int_ent], L_int_ent_id[int_ent], nr_dof_ent);
	int idofent;
	for(idofent=0; idofent<nr_dof_ent; idofent++){
	  printf("dof_ent: id %d, type %d, nrdofs %d\n",
		 l_dof_ent_id[idofent],l_dof_ent_type[idofent],l_dof_ent_nrdof[idofent]);
	}
  }
/*kew*/

#pragma omp atomic
			utv_time.CPU_integration_elems += time_clock() - time_tmp;

#ifdef TIME_TEST_INTEGRATION_ASSEMBLY
#pragma omp atomic
			time_integration += time_clock() - time_tmp;
#endif

			  }


#ifdef DEBUG_SIM
		  if(nrdofs_int_ent>max_nrdofs){
			printf("Too small arrays stiff_mat and rhs_vect passed to comp_el_stiff_mat\n");
			printf("from sir_create in sis_mkb. %d < %d. Exiting !!!",
		   max_nrdofs, nrdofs_int_ent);
			exit(-1);
		  }
#endif


/*kbw
#pragma omp critical(printing)
	{
	 if(L_int_ent_type[int_ent]==PDC_ELEMENT) {
	  int i=6;
	  int solver_id = pdr_ctrl_i_params(Problem_id,i);
	  printf("utr_create_assemble: before assemble:\n");
	  printf("Problem_id %d, Solver_id %d, Level_id %d, sol_typ %d\n",
		 Problem_id, solver_id, Level_id, Comp_type);
	  int ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
	  printf("ient %d, int_ent_id %d, int_ent_type %d, nr_dof_ent %d\n",
		 int_ent, L_int_ent_id[int_ent], L_int_ent_type[int_ent],
		 nr_dof_ent);
	  if(Comp_type!=PDC_NO_COMP){
		pli = 0; nrdof=0;
		for(ibl=0;ibl<nr_dof_ent; ibl++) nrdof+=l_dof_ent_nrdof[ibl];
		for(ibl=0;ibl<nr_dof_ent; ibl++){
		  printf("bl_id %d, bl_nrdof %d\n",
			 l_dof_ent_id[ibl],l_dof_ent_nrdof[ibl]);
		  nri = l_dof_ent_nrdof[ibl];
		  plj=0;
		  for(jbl=0;jbl<nr_dof_ent;jbl++){
	//printf("Stiff_mat transposed! (blocks %d:%d)\n",jbl,ibl);
	//	nrj = l_dof_ent_nrdof[jbl];
	//	for(i=0;i<nri;i++){
	//	  jaux = plj+(pli+i)*nrdof;
	//	  for(j=0;j<nrj;j++){
	//	    printf("%20.15lf",stiff_mat[jaux+j]);
	//
	//	    if(stiff_mat[jaux+j]< -1.e20){
	//	      getchar(); getchar(); getchar();
	//	    }
	//
	///	  }
	//	  printf("\n");
	//	}
		printf("Stiff_mat (blocks %d:%d)\n",jbl,ibl);
		nrj = l_dof_ent_nrdof[jbl];
		  for(j=0;j<nrj;j++){
		for(i=0;i<nri;i++){
		  jaux = plj+(pli+i)*nrdof;
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
	 }
	}
/*kew*/

		double time_tmp = time_clock();

		if(Nr_colors_elems<=1){

#ifdef _OPENMP
		  if(omp_in_parallel()){
			printf("In OpenMP parallel region in utr_create_assemble for 1 color!\n. Exiting.\n");
			exit(0);
		  }
#endif

	  // critical section switched off since for one color we serialize the execution
	  //#pragma omp critical(assembling)
		  {

			if(Assembly_table != NULL){
			  //if(level_p->l_int_ent_type[intent]==3){

			  int position = Asse_pos_first_dof_int_ent[int_ent];
		  int pos_global = Pos_first_dof_int_ent[int_ent];
			  pdr_assemble_local_stiff_mat_with_table(Problem_id, Level_id, Comp_type,
							nr_dof_ent,
							&Assembly_table[position],
							&Local_to_global[pos_global],
							stiff_mat, rhs_vect, &rewrite);


				}
			else{
			  pdr_assemble_local_stiff_mat(Problem_id, Level_id, Comp_type,
						 nr_dof_ent, l_dof_ent_type,
						 l_dof_ent_id,l_dof_ent_nrdof,
						 stiff_mat, rhs_vect, &rewrite);
			}
		  }

		}
		else{ // if more colors - no critical section !!!

//#pragma omp critical(assembling) //if something goes wrong
		  {
			if(Assembly_table != NULL){
			  //if(level_p->l_int_ent_type[intent]==3){

			  //printf("OK - elements with table\n");


			  int position = Asse_pos_first_dof_int_ent[int_ent];
		  int pos_global = Pos_first_dof_int_ent[int_ent];

/*kbw
			  printf("Assembling int_ent %d, using assembly_table at position %d\n",
		   int_ent, position);
/*kew*/

			  pdr_assemble_local_stiff_mat_with_table(Problem_id, Level_id, Comp_type,
							nr_dof_ent,
							&Assembly_table[position],
							&Local_to_global[pos_global],
							stiff_mat, rhs_vect, &rewrite);



			}
			else{

			pdr_assemble_local_stiff_mat(Problem_id, Level_id, Comp_type,
					   nr_dof_ent, l_dof_ent_type,
					   l_dof_ent_id,l_dof_ent_nrdof,
					   stiff_mat, rhs_vect, &rewrite);

			}
			  }

			} // end if more colors

#pragma omp atomic
		utv_time.CPU_assembly_elems += time_clock() - time_tmp;

#ifdef TIME_TEST_INTEGRATION_ASSEMBLY
#pragma omp atomic
		time_assembly += time_clock() - time_tmp;
#endif



	  } // end PARALLEL (OpenMP) loop over elements of a single color - implicit barrier

//#pragma omp barrier

		} // end loop over colors processed in standard way

		free(stiff_mat);
		free(rhs_vect);

	  }// the end of the second parallel region (OpenMP+OpenCL) or fake {

#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)
	  time_concurrent_CPU += time_clock() - t_temp;
#endif


	} // the end of else for if(my_id==0) (OpenMP+OpenCL) or fake {

  } // the end of the first parallel region
#pragma omp barrier

  // after creating CRS structures on CPU and GPU we get one CRS
#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)

  double t_temp = time_clock();

  utr_merge_solver_structures_accel();

  utv_time.CPU_GPU_solver_ds_merging += time_clock() - t_temp;

#endif


/*ok_kbw*/
#ifdef TIME_TEST_INTEGRATION_ASSEMBLY
  printf("\n<<<---***---!!!--- Performance summary data begin in utr_create_assemble ---!!!---***--->>>\n");
  time_total = time_clock() - time_begin;
  printf("Execution time (sum of all threads): integration %lf, assembly %lf\n",
	 time_integration, time_assembly);
  printf("Execution time: total (outside both parallel loops) %lf\n", time_total);
  printf("Execution time (perfect balance): integration %lf, assembly %lf, overhead %lf\n",
	 time_integration/number_of_threads, time_assembly/number_of_threads,
	 time_total - (time_integration + time_assembly)/number_of_threads );
  printf("(we take the number of CPU threads integrating elements not faces! - inaccuracy...)\n");
#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)
  printf("Concurrent execution times: GPU %lf, CPU %lf\n",
	 time_accelerator, time_concurrent_CPU);
#endif
  printf("<<<---***---!!!--- Performance summary data end in utr_create_assemble ---!!!---***-->>>\n\n");
#endif
/*kew*/


  return(1);
}



	/////////////////************ RENUMBER UTILITIES ***********//////////////

typedef struct list {
	int row;
	struct list *next;
}list_row;

/*------------------------------------------------------------
  utr_rcm_renumber - Reverse Cuthill McKee bandwidth reduction algorithm (renumbering)
------------------------------------------------------------*/
void utr_rcm_renumber(
  int ** l_neig,  /* in: adjacency list */
  int * nrneig,   /* in: list of numbers of neighbors */
  int nr_dof_ent,
  int *rpermutation  /* out: permutation vector */
			  );

/*------------------------------------------------------------
  utr_renumber - to renumber (permute) graph vertices
------------------------------------------------------------*/
void utr_renumber(
  int Nr_vertices, /* in: number of graph vertices */
  int * Nrneig,    /* in: list of numbers of neighbors */
  int ** L_neig,   /* in: adjacency list */
  int *Permutation_array  /* out: permutation vector */
			  )
{

  double time_begin = time_clock();

  // Reverse Cuthill McKee algorithm is used as default
  utr_rcm_renumber(L_neig, Nrneig, Nr_vertices, Permutation_array);

  utv_time.renumber_dofs += time_clock()-time_begin;


}

/*------------------------------------------------------------
  utr_rcm_renumber - Reverse Cuthill McKee bandwidth reduction algorithm (renumbering)
------------------------------------------------------------*/
void utr_rcm_renumber(
  int ** l_neig,  /* in: adjacency list */
  int * nrneig,   /* in: list of numbers of neighbors */
  int nr_dof_ent,
  int *rpermutation  /* out: permutation vector */
){

  list_row *head=NULL, *tail=NULL, *lhead=NULL, *ltail=NULL, *elem=NULL, *lelem=NULL;
  int sum_nreig, icrs=0, jcrs=0;

  /* array of visited DOF - 1 visited, 0 not visited */
  int *vt=NULL, *permutation=NULL;

  permutation =(int *)malloc( nr_dof_ent*sizeof(int) );
  if(permutation==NULL){printf("no memory in utr_reverse_cuthill_mckee permutation"); exit(0);}
  if(rpermutation==NULL){printf("no memory in utr_reverse_cuthill_mckee rpermutation"); exit(0);}

  vt =(int *)malloc( nr_dof_ent*sizeof(int));
  if(vt==NULL){printf("no memory in utr_reverse_cuthill_mckee vt"); exit(0);}

  int ibl,ineig,i,start_dof=0,iper;
  int iaux,idofent;

  sum_nreig=nrneig[start_dof];
  iper=nr_dof_ent-1;
  for(ibl = 0; ibl< nr_dof_ent; ibl++){vt[ibl]=0;(rpermutation)[ibl]=-1;permutation[iper]=-1;}


  for(ibl = 1; ibl< nr_dof_ent; ibl++){
	sum_nreig+=nrneig[ibl];
	if(nrneig[ibl] <= nrneig[start_dof]) {start_dof=ibl;}
  }


  /* start compute permutation vector */
  vt[start_dof]=1;
  // permutation[iper]=start_dof;
  // (rpermutation)[start_dof]=iper;
  (rpermutation)[iper]=start_dof;
  permutation[start_dof]=iper;
  iper=iper-1;
  while(iper>=0){
	for(ineig=0;ineig<nrneig[start_dof];ineig++){
	  iaux=l_neig[start_dof][ineig];


	  if(vt[iaux]==0){
  vt[iaux]+=1;
  elem =(list_row*)malloc(sizeof(list_row));
  elem->row=iaux;
  elem->next=NULL;
  if(lhead==NULL){
	lhead = elem;
	ltail = elem;
  }
  else{
	if(nrneig[iaux] < nrneig[lhead->row]){
	  elem->next = lhead;
	  lhead = elem;
	}else{
	  lelem = lhead;
	  while((lelem->next != NULL) &&
	  (nrneig[iaux] >= nrneig[lelem->row])){
		lelem = lelem->next;
	  }
	  elem->next = lelem->next;
	  lelem->next = elem;
	  if(lelem == ltail)ltail=elem;
	}
  }
	  }

	}

	if(head == NULL) { head=lhead; tail=ltail; lhead=NULL; ltail=NULL; }
	else if(lhead!=NULL){ tail->next = lhead; tail=ltail; lhead=NULL; ltail=NULL;}

	if(head == NULL){ //&& iper>=0){
	  for(start_dof=0;vt[start_dof]==1;++start_dof);
	  for(ibl = start_dof+1; ibl< nr_dof_ent; ibl++){
		if(vt[ibl]==0 && nrneig[ibl] <= nrneig[start_dof]) {start_dof=ibl;}
	  }
	  vt[start_dof]+=1;
	}
	else{
	  start_dof=head->row;
	  elem=head;
	  head=head->next;
	  if(head==NULL)tail=NULL;
	  free(elem);
	}

	// permutation[iper]=start_dof;
	// (rpermutation)[start_dof]=iper;
	(rpermutation)[iper]=start_dof;
	permutation[start_dof]=iper;
	iper=iper-1;

  }
  /*
  printf("permutacja:%d \n",nr_dof_ent);
  for(idofent = 0; idofent< nr_dof_ent; idofent++)printf("%d ",permutation[idofent]);
  printf("\nrpermutacja:\n");
  for(idofent = 0; idofent< nr_dof_ent; idofent++)printf("%d ",(rpermutation)[idofent]);
  /**/

  // setting block_id, pos_glob, l_neig_bl

  /*
  start_dof=0;
  for(idofent = 0; idofent< nr_dof_ent; idofent++){
	// block_id is offset 0 !!!
	l_dof_struct[rpermutation[idofent]].block_id = idofent;
	l_dof_struct[rpermutation[idofent]].posglob=start_dof;
   start_dof+=l_dof_struct[rpermutation[idofent]].nrdofs;

  }
  */


  free(vt);
  free(permutation);
  //free(rpermutation);

}

	/////////////////************ COLORING UTILITIES ***********//////////////

/*---------------------------------------------------------
  utr_color_int_ent_for_assembly_with_graph_creation - coloring elements for lock free
					  assembly (with creation of neighbourhood graphs)
----------------------------------------------------------*/
int utr_color_int_ent_for_assembly_with_graph_creation(
  int Problem_id,           /* in: Problem id */
  int Level_id,             /* in: Level id */
  int nr_int_ents,             /* in: number of integration entities (int_ents,faces,elems+faces) */
  int offset,               /* in: offset in L_int_ent... structures */
  int* L_int_ent_type,      /* in: integration entities type,
							   out: coloured  */
  int* L_int_ent_id,        /* in: integration entities id,
							   out: coloured  */
  int* Nr_colors_int_ents,     /* out: number of colours  */
  int** L_color_index_int_ents /* out: crs table with index of L_int_ent.. where colour starts and ends  */
				   )
{

  // i.e. maximal number of DOF blocks (nodes) per element
  int l_dof_ent_types[SIC_MAX_DOF_PER_INT];
  int l_dof_ent_ids[SIC_MAX_DOF_PER_INT];
  int l_dof_ent_nrdofs[SIC_MAX_DOF_PER_INT];
  int i,mesh_id,j,k,int_ent_current,node_current,iux,dux,jj,flag;


  const int maxnode_neig=4*PDC_MAX_DOF_PER_INT;  //one node belongs max to maxnode_neig integration entities

  //int nr_elemk = mmr_get_nr_elem(mesh_id);
  int nr_int_entk =nr_int_ents;
  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  int nr_node = mmr_get_nr_node(mesh_id) + mmr_get_nr_edge(mesh_id);

  int *int_ent_color;

  int *int_ent_list_id=(int*)malloc(nr_int_entk * sizeof(int));
  int *int_ent_list_type=(int*)malloc(nr_int_entk * sizeof(int));

  //printf("before alocation in coloring \n\n");

  int **int_ent2node=(int**)malloc(nr_node * sizeof(int*));
  if(int_ent2node==NULL){printf("no memory in utr_color_int_ent_for_assembly: int_ent2node\n"); exit(0);}
  for(i=0;i<nr_node;i++){
	int_ent2node[i]=(int*)malloc((maxnode_neig+1) * sizeof(int));
	if(int_ent2node[i]==NULL){printf("no memory in utr_color_int_ent_for_assembly: int_ent2node[%d]\n",i); exit(0);}
	int_ent2node[i][0]=0;
  }

  int **int_ent2int_ent=(int**)malloc(nr_int_entk * sizeof(int*));
  if(int_ent2int_ent==NULL){printf("no memory in utr_color_int_ent_for_assembly: int_ent2int_ent\n"); exit(0);}
  for(i=0;i<nr_int_entk;i++){
	int_ent2int_ent[i]=(int*)malloc(maxnode_neig * sizeof(int));
	if(int_ent2int_ent[i]==NULL){printf("no memory in utr_color_int_ent_for_assembly: int_ent2int_ent[%d]\n",i); exit(0);}
	for(j=0;j<maxnode_neig;j++)int_ent2int_ent[i][j]=0;
  }

  /*
	for(i=0;i<nr_int_ents;i++){
	printf("%d ",L_int_ent_id[i]);
	}
	printf("\n");
  */

  for(i=0;i<nr_int_entk;i++){

	int nr_dof_ent_loc = SIC_MAX_DOF_PER_INT; // maximal number of DOF ents per INT ent
	pdr_comp_stiff_mat(Problem_id, L_int_ent_type[i], L_int_ent_id[i], PDC_NO_COMP, NULL,
			   &nr_dof_ent_loc, l_dof_ent_types,
			   l_dof_ent_ids, l_dof_ent_nrdofs,
			   NULL, NULL, NULL, NULL);

	int_ent_list_id[i]=L_int_ent_id[i];
	int_ent_list_type[i]=L_int_ent_type[i];
	int_ent_current=i;
	for(j=0;j<nr_dof_ent_loc;j++){
	  //printf("%d -> %d\n",i,(l_dof_ent_ids[j]-1));
	  node_current=l_dof_ent_ids[j]-1;

	  //#ifdef DEBUG
	  if(node_current>=nr_node){
		printf("no memory allocated in utr_color_int_ent_for_assembly: int_ent2node\n wrong value of nr_node\n");
		exit(0);}
	  //#endif

	  ++int_ent2node[node_current][0];
	  iux=int_ent2node[node_current][0];

	  //#ifdef DEBUG
	  if(iux>=maxnode_neig){
		printf("no memory allocated in utr_color_int_ent_for_assembly: int_ent2node\n increase value of maxnode_neig\n");
		exit(0);}
	  if(int_ent_current>=nr_int_entk){
		printf("no memory allocated in utr_color_int_ent_for_assembly: int_ent2int_ent\n wrong value of nr_int_entk\n");
		exit(0);}
	  //#endif

	  int_ent2node[node_current][iux]=int_ent_current;
	  if(iux>1){
		for(k=1;k<iux;k++){
		  dux=int_ent2node[node_current][k];

		  for(jj=1,flag=0;jj<=int_ent2int_ent[int_ent_current][0];jj++){
			if(int_ent2int_ent[int_ent_current][jj]==dux){flag=1;break;}
		  }
		  if(flag==0){
			++int_ent2int_ent[int_ent_current][0];
			jj=int_ent2int_ent[int_ent_current][0];
			int_ent2int_ent[int_ent_current][jj]=dux;
			++int_ent2int_ent[dux][0];
			jj=int_ent2int_ent[dux][0];
			int_ent2int_ent[dux][jj]=int_ent_current;
		  }

		}
	  }

	} /* end loop over nodes in int_ent */
  }

  // for(i=0;i<nr_int_entk;i++){
  // for(j=0;j<nr_int_entk;j++)printf("%d ",int_ent2int_ent[i][j]);
  // printf("\n");

  // }

  for(i=0;i<nr_node;i++) free(int_ent2node[i]);
  free(int_ent2node);

  int_ent_color=(int*)malloc(nr_int_entk * sizeof(int));

  *Nr_colors_int_ents=utr_generate_int_ent_coloring(nr_int_entk,int_ent2int_ent,int_ent_color);

/*kbw
  printf("Nr_colors = %d\n", *Nr_colors_int_ents);
/*kew*/

  //for(i=0;i<nr_int_entk;i++)printf(" %d=%d,%d ",i,int_ent_color[i],int_ent_list_id[i]);printf("\n");

  for(i=0;i<nr_int_entk;i++) free(int_ent2int_ent[i]);
  free(int_ent2int_ent);

  //int *int_ent_color_crs_val=(int*)malloc(nr_int_entk * sizeof(int));
  *L_color_index_int_ents=(int*)malloc((*Nr_colors_int_ents+1) * sizeof(int));


  //convert color vector to crs
  (*L_color_index_int_ents)[0]=offset;
  for(i=0,k=0;i< (*Nr_colors_int_ents);i++){
	for(j=0;j<nr_int_entk;j++){
	  if(int_ent_color[j]==i){
		//int_ent_color_crs_val[k]=j;
		L_int_ent_id[k]=int_ent_list_id[j];
		L_int_ent_type[k]=int_ent_list_type[j];
		k++;
	  }
	}
	(*L_color_index_int_ents)[i+1]=k+offset;
  }
  /*kbw
	for(i=0;i< (*Nr_colors_int_ents);i++){
	for(j=(*L_color_index_int_ents)[i];j<(*L_color_index_int_ents)[i+1];j++)printf("%d ",L_int_ent_id[j]);
	printf("\n");
	}
	/*kew*/
  free(int_ent_list_id);
  free(int_ent_list_type);
  free(int_ent_color);

  return(0);
}


/*---------------------------------------------------------
  utr_color_int_ent_for_assembly - coloring elements for lock free assembly
----------------------------------------------------------*/
int utr_color_int_ent_for_assembly(
  int Problem_id,           /* in: Problem id */
  int Level_id,             /* in: Level id */
  int Nr_elems,             /* in: number of elements */
  int Nr_faces,             /* in: number of faces */
  int* L_int_ent_type,      /* in: integration entities type,
							   out: coloured  */
  int* L_int_ent_id,        /* in: integration entities id,
							   out: coloured  */
  int Nr_dof_ent,           /* number of dof entities dof_ents */
  int* Nr_int_ent_loc,      /* number of int_ents for each dof_ent */
  int** L_int_ent_index,    /* list of indices in int_ent table for each dof_ent */
  int* Nr_colors_elems,     /* out: number of colours for elements  */
  int** L_color_index_elems,/* out: crs table with index of L_int_ent.. where colour starts and ends  */
  int* Nr_colors_faces,     /* out: number of colours for faces  */
  int** L_color_index_faces /* out: crs table with index of L_int_ent.. where colour starts and ends  */
				   )
{


  char field_module_name[100];

  // currently supported fields:
  // STANDARD_LINEAR - for continuous, vector, linear approximations
  // STANDARD_QUADRATIC - for continuous, vector, linear and quadratic and mixed approximations
  // DG_SCALAR_PRISM - for discontinuous, scalar, high order approximations
  apr_module_introduce(field_module_name);

  /* new coloring */
  utr_color_int_ent_for_assembly_with_std_vector(
				   Problem_id, Level_id, Nr_elems, Nr_faces,
							   L_int_ent_type,
							   L_int_ent_id,
							   Nr_dof_ent,
							   Nr_int_ent_loc,
							   L_int_ent_index,
							   Nr_colors_elems,
							   L_color_index_elems,
							   Nr_colors_faces,
							   L_color_index_faces);
  return 0;
  /* new coloring */

  //
  // OLD/NEW COLORING IN ERROR CASE
  //

  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0 ){

  /**/
	utr_color_int_ent_for_assembly_with_std_vector(
				   Problem_id, Level_id, Nr_elems, Nr_faces,
							   L_int_ent_type,
							   L_int_ent_id,
							   Nr_dof_ent,
							   Nr_int_ent_loc,
							   L_int_ent_index,
							   Nr_colors_elems,
							   L_color_index_elems,
							   Nr_colors_faces,
							   L_color_index_faces);

  /**/
  }
  else if( strncmp(field_module_name,"STANDARD_QUADRATIC",18) == 0){

	int field_id = pdr_ctrl_i_params(Problem_id, 3);
	int pdeg = apr_get_el_pdeg(field_id,0,NULL);

	if(pdeg == APC_LINEAR_APPROXIMATION_PDEG){

  /**/
	utr_color_int_ent_for_assembly_with_std_vector(
				   Problem_id, Level_id, Nr_elems, Nr_faces,
							   L_int_ent_type,
							   L_int_ent_id,
							   Nr_dof_ent,
							   Nr_int_ent_loc,
							   L_int_ent_index,
							   Nr_colors_elems,
							   L_color_index_elems,
							   Nr_colors_faces,
							   L_color_index_faces);

  /**/


	   }
	else{

	double time_begin = time_clock();
	double time_tmp = time_clock();

  // obsolete version
  int i,j,k,l,flag;

  int nr_int_ent=Nr_elems+Nr_faces;

  int **int_ent_neig=(int**)malloc(nr_int_ent * sizeof(int*));
  if(int_ent_neig==NULL){
	printf("no memory in utr_color_int_ent_for_assembly: int_ent_neig\n");
	exit(0);
  }

  //#pragma omp parallel for // does not help - do not know why
  for(i=0;i<nr_int_ent;i++){

	int_ent_neig[i]=(int*)malloc((UTC_MAX_INT_ENT_NEIG+1) * sizeof(int));
	if(int_ent_neig[i]==NULL){
	  printf("no memory in utr_color_int_ent_for_assembly: faces_neig[%d]\n",i);
	  exit(0);
	}

	for(j=1;j<=UTC_MAX_INT_ENT_NEIG;j++) int_ent_neig[i][j]=UTC_LIST_END_MARK;
	int_ent_neig[i][0]=0;

  }

/*kbw*/
  printf("!!!!!!!!!!! Time for allocating graph for coloring: %lf\n",
	 time_clock()-time_tmp);
/*kew*/
  time_tmp = time_clock();

  //printf("\n Nr_elems=%d Nr_faces=%d\n",Nr_elems,Nr_faces);

  //#pragma omp parallel for // does not help - do not know why
  for(i=0;i<Nr_dof_ent;i++){ // for each dof entity

	for(j=0;j<Nr_int_ent_loc[i];j++){ // for each integration entity containing dof entity

	  int int_ent_current=L_int_ent_index[i][j]; // index of integration ent

	  //if(L_int_ent_type[int_ent_current]==PDC_ELEMENT){
	  if(int_ent_current<Nr_elems){ //element


	for(k=0;k<Nr_int_ent_loc[i];k++){
	//for(k=j+1;k<Nr_int_ent_loc[i];k++){


	  int neig_index = L_int_ent_index[i][k];

	  //if(k!=j && L_int_ent_type[neig_index]==PDC_ELEMENT){
	  //if(L_int_ent_type[neig_index]==PDC_ELEMENT){
	  if(neig_index<Nr_elems){ //element

		int iaux = utr_put_list(neig_index,
					&int_ent_neig[int_ent_current][1], UTC_MAX_INT_ENT_NEIG);

		if(iaux<0) int_ent_neig[int_ent_current][0]++;

		if(iaux==0){ // list full - increase UTC_MAX_INT_ENT_NEIG
		  printf(" list full - increase UTC_MAX_INT_ENT_NEIG!!! Exiting");
		  exit(-1);
		}

	  }

	}

	  }
	  else{ // if face


	for(k=0;k<Nr_int_ent_loc[i];k++){
	//for(k=j+1;k<Nr_int_ent_loc[i];k++){


	  int neig_index = L_int_ent_index[i][k];

	  //if(k!=j && L_int_ent_type[neig_index]==PDC_FACE){
	  //if(L_int_ent_type[neig_index]==PDC_FACE){
	  if(neig_index>=Nr_elems){ // face

		int iaux = utr_put_list(neig_index-Nr_elems,
					&int_ent_neig[int_ent_current][1], UTC_MAX_INT_ENT_NEIG);


		if(iaux<0) int_ent_neig[int_ent_current][0]++;

		if(iaux==0){ // list full - increase UTC_MAX_INT_ENT_NEIG
		  printf(" list full - increase UTC_MAX_INT_ENT_NEIG!!! Exiting");
		  exit(-1);
		}

	  }

	}

	  }

	}

  }

/*kbw*/
  printf("!!!!!!!!!!! Time for creating graph for coloring: %lf\n",
	 time_clock()-time_tmp);
/*kew*/


/*   for(i=0;i<Nr_dof_ent;i++){ */
/* /\*kbw */
/*   printf("dof_ent %d, number of int_ents %d, list: \n", */
/*   i, Nr_int_ent_loc[i]); */
/* /\*kew*\/ */
/*     for(j=0;j<Nr_int_ent_loc[i];j++){ */
/* /\*kbw */
/*   printf("%10d  ", L_int_ent_index[i][j]); */
/* /\*kew*\/ */
/*       int int_ent_current=L_int_ent_index[i][j]; */
/*       if(int_ent_current<Nr_elems){ //element */
/* 	for(k=j+1;k<Nr_int_ent_loc[i];k++){ */
/* 	  if(L_int_ent_index[i][k]<Nr_elems){ */
/* 	    flag=0; */
/* 	    for(l=0;l<int_ent_neig[int_ent_current][0];++l){ */
/* 	      if(int_ent_neig[int_ent_current][l+1]==L_int_ent_index[i][k]){flag=1;break;} */
/* 	    } */
/* 	    if(flag==0){ */
/*               ++int_ent_neig[int_ent_current][0]; */
/*               int_ent_neig[int_ent_current][int_ent_neig[int_ent_current][0]]=L_int_ent_index[i][k]; */
/*             } */
/*           } */
/*         } */

/*       }else{ //face */
/*         for(k=j+1;k<Nr_int_ent_loc[i];k++){ */
/*           if(L_int_ent_index[i][k]>=Nr_elems){ */
/*             flag=0; */
/*             for(l=0;l<int_ent_neig[int_ent_current][0];++l){ */
/*               if(int_ent_neig[int_ent_current][l+1]==(L_int_ent_index[i][k]-Nr_elems)){flag=1;break;} */
/*             } */
/*             if(flag==0){ */
/*               ++int_ent_neig[int_ent_current][0]; */
/*               int_ent_neig[int_ent_current][int_ent_neig[int_ent_current][0]]=L_int_ent_index[i][k]-Nr_elems; */
/*             } */
/*           } */
/*         } */
/*       } */

/*       if(int_ent_neig[int_ent_current][0]>=UTC_MAX_INT_ENT_NEIG){ // list full - increase UTC_MAX_INT_ENT_NEIG *\/ */
/* 	printf(" list full - increase UTC_MAX_INT_ENT_NEIG!!! Exiting"); */
/* 	exit(-1); */
/*       } */

/*     } */
/* /\*kbw */
/*       printf("\n"); */
/* /\*kew*\/ */

/*   } */

  /*kbw
  for(i=0;i<(Nr_faces+Nr_elems);i++){
	printf("%d=",i);
	for(j=0;j<=int_ent_neig[i][0];j++){
	  printf("%d ",int_ent_neig[i][j]);

	}printf("\n");
  }
  /*kew*/

  int* elem_color=(int*)malloc(Nr_elems * sizeof(int));
  int* face_color=(int*)malloc(Nr_faces * sizeof(int));

  time_tmp = time_clock();
  *Nr_colors_elems=utr_generate_int_ent_coloring(Nr_elems,int_ent_neig,elem_color);
  utv_time.color_elements += time_clock()-time_tmp;


/*kbw*/
  printf("\nNr_colors_for_elems= %d (creation time %lf)\n",
	 *Nr_colors_elems, time_clock()-time_tmp);
/*kew*/

  time_tmp = time_clock();
  *Nr_colors_faces=utr_generate_int_ent_coloring_old(Nr_faces,&int_ent_neig[Nr_elems],face_color);
  utv_time.color_faces += time_clock()-time_tmp;

/*kbw*/
  printf("Nr_colors_for_faces= %d (creation time %lf)\n",
	 *Nr_colors_faces, time_clock()-time_tmp);
/*kew*/

  for(i=0;i<nr_int_ent;i++) free(int_ent_neig[i]);
  free(int_ent_neig);


  *L_color_index_elems=(int*)malloc((*Nr_colors_elems+1) * sizeof(int));
  *L_color_index_faces=(int*)malloc((*Nr_colors_faces+1) * sizeof(int));

  int *int_ent_list_id=(int*)malloc(nr_int_ent * sizeof(int));
  int *int_ent_list_type=(int*)malloc(nr_int_ent * sizeof(int));
  for(i=0;i<nr_int_ent;i++){
	int_ent_list_id[i]=L_int_ent_id[i];
	int_ent_list_type[i]=L_int_ent_type[i];
  }

  time_tmp = time_clock();
  //convert colors vector for elements to crs
  (*L_color_index_elems)[0]=0;
  for(i=0,k=0;i< (*Nr_colors_elems);i++){
	for(j=0;j<Nr_elems;j++){
	  if(elem_color[j]==i){
		L_int_ent_id[k]=int_ent_list_id[j];
		L_int_ent_type[k]=int_ent_list_type[j];

		//int_ent_list_id[k]=L_int_ent_id[j];
		//int_ent_list_type[k]=L_int_ent_type[j];

		k++;
	  }
	}
	(*L_color_index_elems)[i+1]=k;
  }

  //convert colors vector for faces to crs
  (*L_color_index_faces)[0]=Nr_elems;
  for(i=0,k=Nr_elems;i< (*Nr_colors_faces);i++){
	for(j=0;j<Nr_faces;j++){
	  if(face_color[j]==i){
		L_int_ent_id[k]=int_ent_list_id[j+Nr_elems];
		L_int_ent_type[k]=int_ent_list_type[j+Nr_elems];

		//int_ent_list_id[k]=L_int_ent_id[j+Nr_elems];
		//int_ent_list_type[k]=L_int_ent_type[j+Nr_elems];

		k++;
	  }
	}
	(*L_color_index_faces)[i+1]=k;
  }

/*kbw*/
  printf("!!!!!!!!!!! Time for converting colors to CRS: %lf\n",
	 time_clock()-time_tmp);
/*kew*/

  free(elem_color);
  free(face_color);
  free(int_ent_list_id);
  free(int_ent_list_type);

  utv_time.total_coloring += time_clock()-time_begin;

	} // end if pdeg not linear

  } // end if quadratic module
  else  { // for DG approximation and others
/*kbw*/
  printf("!!!!!!!!!!! Coloring not supported for this approximation type !!!!!!!!!!!!!\n");
/*kew*/
  }


  return(0);
}

