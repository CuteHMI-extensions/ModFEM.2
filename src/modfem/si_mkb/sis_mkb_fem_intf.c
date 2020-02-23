/*************************************************************
File contains procedures:
 fem_proj_sol_lev - to L2 project solution dofs between mesh levels
 fem_vec_norm - to compute a norm of global vector in parallel
 fem_sc_prod - to compute a scalar product of two global vectors
 fem_exchange_dofs - to exchange dofs between processors

History:
	05.2001 - Krzysztof Banas, initial version

*************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <string.h>

/* problem dependent module with implementation of pdr_ routines */
#include <modfem/pdh_intf.h>

/* general solver interface */
#include <modfem/sih_intf.h>

/* linear solver (adapter) required interface */
#include "../../../include/modfem/si_mkb/sih_mkb_fem_intf.h"

/* local declarations and data structure for solver interface module */
#include "../../../include/modfem/si_mkb/sih_mkb.h"

/* LAPACK procedures */
#include <modfem/lin_alg_intf.h>

// for apr_module_introduce
#include <modfem/aph_intf.h>

/*---------------------------------------------------------
  fem_proj_sol_lev - to L2 project solution dofs between mesh levels
---------------------------------------------------------*/
int fem_proj_sol_lev( /* returns: >=0 - success; <0 - error code*/
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,        /* in: solver data structure to be used */
  int Ilev_from,    /* in: level number to project from */
  double* Vec_from, /* in: vector of values to project */
  int Ilev_to,      /* in: level number to project to */
  double* Vec_to    /* out: vector of projected values */
  )
{

#define COARSE_TO_FINE 0
#define FINE_TO_COARSE 1

#define MAX_NUM_ANC     8  /* maximal number of sons for element */
#define MAX_EL_DOFS     1000    /* maximal number of dofs for element */
#define MAX_PROJ_COMB   2 /* maximal number of pdeg combinations */

  sit_levels *level_coarse, *level_fine;
  int *l_bl_nrdof_coarse, *l_bl_posg_coarse;
  int *l_bl_nrdof_fine, *l_bl_posg_fine;
  int nrdofgl, nrblocks, max_nr_bl, problem_id;
  int mode; /* mode of operation: 0 (COARSE_TO_FINE) - from coarse to fine, */
			/*                    1 (FINE_TO_COARSE)- from fine to coarse */

  /* auxiliary variables */
  static int nr_precomp=0;
  static int precomp_ind[MAX_PROJ_COMB];
  static double
	shap_fun_precomp[MAX_PROJ_COMB][MAX_NUM_ANC][MAX_EL_DOFS][MAX_EL_DOFS];
  int i, iaux, iel_coarse, iel_fine;
  int posglob_coarse, posglob_fine, pdeg_coarse, pdeg_fine;
  int ndof_coarse, ndof_fine, ibl_coarse, ibl_fine;
  int idofs, jdofs, index, ison, fath, el_sons[MAX_NUM_ANC+1];
  double shap_fun_proj[MAX_EL_DOFS];
  int dof_struct_fine_id=0,dof_struct_temp_id=0;
  sit_dof_struct dof_struct_fine;
  sit_dof_struct dof_struct_temp;
/*++++++++++++++++ executable statements ++++++++++++++++*/

  problem_id = siv_solver[Solver_id].problem_id;

  // check the name of the field module
  char field_module_name[100];
  apr_module_introduce(field_module_name); // BROKEN INDEPENDENCE OF SID_
  // SHOULD BE CHANGED FOR INCLUDING TYPE OF APPROXIMATION IN SOLVER DATA

  // for standard continuous linear FEM approximation
  if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0 || strncmp(field_module_name,"STANDARD_QUADRATIC",18) == 0){

	printf("fem_proj_sol_lev NOT IMPLEMENTED for STD and QUAD approximations!");
	exit (-1);

  }
  else if( strncmp(field_module_name,"DG_SCALAR_PRISM ",15) == 0) { // for DG approximation

/* projecting from coarser to finer grid */
  if(Ilev_from<Ilev_to) {

	mode = COARSE_TO_FINE;
	level_coarse = &siv_solver[Solver_id].level[Ilev_from];
	level_fine = &siv_solver[Solver_id].level[Ilev_to];

  }
  else{

	mode = FINE_TO_COARSE;
	level_coarse = &siv_solver[Solver_id].level[Ilev_to];
	level_fine = &siv_solver[Solver_id].level[Ilev_from];

  }

  /* checking */
  if(level_fine->max_dofs_dof_ent > MAX_EL_DOFS){
	printf("Error 2374 in number of dofs in proj_sol!\n");
	exit(-1);
  }
  if(level_coarse->max_dofs_dof_ent > MAX_EL_DOFS){
	printf("Error 2374 in number of dofs in proj_sol!\n");
	exit(-1);
  }

  /*!!! DG: simple algorithm - DOF blocks are associated with elements only */
  /* for each OWNED element active in a coarse mesh*/
  // WARNING: number of DOF blocks is equal to the number of DOF entities

  // OLD VERSION
  //for(ibl_coarse=0;ibl_coarse<level_coarse->nr_dof_ent;ibl_coarse++){

  //int dof_struct_id = level_coarse->l_bl_to_struct[ibl_coarse];

  // NEW VERSION
  int dof_struct_id;
  for (dof_struct_id=0;dof_struct_id<level_coarse->nr_dof_ent;dof_struct_id++){

	sit_dof_struct dof_struct = level_coarse->l_dof_struct[dof_struct_id];
	iel_coarse = dof_struct.dof_ent_id;

	/* find list of descendents belonging to the next mesh */
	//mmr_el_fam(mesh_id,iel_coarse,el_sons,NULL);
	int is_owned = pdr_dof_ent_sons(problem_id, PDC_ELEMENT, iel_coarse, el_sons);

	if(is_owned==1){ // only for owned elements

	/* position in global vectors of unknowns */
	posglob_coarse = dof_struct.posglob;

	/* number of dofs associated with iel in the coarse mesh */
	ndof_coarse = dof_struct.nrdofs;

/*kbw
	printf("Coarse: iel %d, dof_struct %d, ibl %d, posg %d, ndof %d\n",
	   iel_coarse, dof_struct_id, ibl_coarse, posglob_coarse,
	   ndof_coarse);
/*kew*/


	if(level_coarse->pdeg_coarse == SIC_PDEG_COARSE_HIGH){
	  /* get sons pdegs and choose the highest */

	}
	else if(level_coarse->pdeg_coarse == SIC_PDEG_COARSE_LOW){
	  /* get sons pdegs and choose the lowest */

	}
	else {
	  /* force pdeg_coarse to be the same for all elements */
	  pdeg_coarse = level_coarse->pdeg_coarse;
	}


	/* first check whether the elements do not belong to both meshes */
	if(level_fine->l_dof_elem_to_struct[iel_coarse]>-1){

	  iel_fine = iel_coarse;
	  dof_struct_fine_id = level_fine->l_dof_elem_to_struct[iel_coarse];
	  dof_struct_fine = level_fine->l_dof_struct[dof_struct_id];
	  ibl_fine = dof_struct_fine.block_id;

	  if(dof_struct_fine.nrdofs!=ndof_coarse){

	printf("Dof_ent %d in two meshes but with different pdeg in proj_sol.\n",
		  iel_fine);
	printf("Not implemented yet. Exiting!\n");
	exit(-1);
	  }
	  else{


#ifdef DEBUG_SIM
	/* check */
	//if(level_fine->l_bl_to_struct[ibl_fine]!=dof_struct_id){
	//  printf("Error in renumbering of elements in proj_sol!\n");
	//  exit(-1);
	//}
	if(dof_struct_fine.nrdofs!=ndof_coarse){
	  printf("Error in number of dofs in proj_sol!\n");
	  exit(-1);
	}
#endif

	/* position in global vectors of unknowns */
	posglob_fine = dof_struct_fine.posglob;
	ndof_fine = ndof_coarse;


	  /* projecting from coarser to finer grid */
	if(mode == COARSE_TO_FINE){

/*kbw
	  printf("Projecting from iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		 iel_coarse, Ilev_from, ibl_coarse, posglob_coarse,
		 ndof_coarse);
	  printf("Projecting to iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		 iel_fine, Ilev_to, ibl_fine, posglob_fine,
		 ndof_fine);
/*kew*/

	  for(i=0;i<ndof_coarse;i++) Vec_to[posglob_fine+i]
					   = Vec_from[posglob_coarse+i];


/*kbw
	  for(i=0;i<ndof_coarse;i++) printf("%20.15lf",
						Vec_to[posglob_fine+i]);
	  printf("\n");
/*kew*/

	}
/* if projecting from  finer to coarser mesh */
	else if(mode == FINE_TO_COARSE){

/*kbw
	  printf("Projecting from iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		 iel_fine, Ilev_from, ibl_fine, posglob_fine,
		 ndof_fine);
	  printf("Projecting to iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		 iel_coarse, Ilev_to, ibl_coarse, posglob_coarse,
		 ndof_coarse);
/*kew*/
	  for(i=0;i<ndof_coarse;i++) Vec_to[posglob_coarse+i]
					   = Vec_from[posglob_fine+i];
/*kbw
	  for(i=0;i<ndof_coarse;i++) printf("%20.15lf",
						Vec_to[posglob_coarse+i]);
	  printf("\n");
/*kew*/
	}
	else{
	  printf("Unknown action specified in proj_sol!\n");
	  exit(-1);
	}

	  } // end if element in both meshes with the same pdeg
	}
	else{ /* element do not belong to both meshes */

	  /* when rewriting from fine mesh - initialize coarse dofs */
	  if(mode==FINE_TO_COARSE) {
	for(i=0;i<ndof_coarse;i++) Vec_to[posglob_coarse+i] = 0.0;
	  }


/*kbw
	  printf("sons:");
	  for(i=0;i<el_sons[0];i++){
	printf("%d  ", el_sons[i+1]);
	  }
	  printf("\n");
/*kew*/

	  /* loop over all sons */
	  for(ison=0;ison<el_sons[0];ison++){

	iel_fine=el_sons[ison+1];

/*kbw
	printf("Fine: iel %d (son %d), dof_struct_id %d\n",
		   iel_fine, ison, level_fine->l_dof_elem_to_struct[iel_fine]);
/*kew*/

	dof_struct_fine_id = level_fine->l_dof_elem_to_struct[iel_fine];
	dof_struct_fine =
	  level_fine->l_dof_struct[dof_struct_fine_id];
	ibl_fine = dof_struct_fine.block_id;


	/* position in global vectors of unknowns */
	posglob_fine = dof_struct_fine.posglob;

	/* number of dofs associated with iel in the coarse mesh */
	ndof_fine = dof_struct_fine.nrdofs;

	if(level_fine->pdeg_coarse == SIC_PDEG_FINEST){
	  pdeg_fine = pdr_get_ent_pdeg(problem_id, PDC_ELEMENT, iel_fine);
	  //pdeg_fine = apr_get_ent_pdeg(field_id, APC_ELEMENT, iel_fine);
	}
	else if(level_fine->pdeg_coarse == SIC_PDEG_COARSE_HIGH){
	  /* get sons pdegs and choose the highest */

	}
	else if(level_fine->pdeg_coarse == SIC_PDEG_COARSE_LOW){
	  /* get sons pdegs and choose the lowest */

	}
	else {
	  /* force pdeg_coarse to be the same for all elements */
	  pdeg_fine = level_fine->pdeg_coarse;
	}


/*kbw
	printf("el_c %d (bl %d, posg %d, nrdofs %d), p_c %d, \n",
		  iel_coarse, ibl_coarse, posglob_coarse, ndof_coarse, pdeg_coarse);
	printf("son %d, el_f %d (bl %d, posg %d, nrdofs %d), p_f %d\n",
		   ison, iel_fine, ibl_fine, posglob_fine, ndof_fine, pdeg_fine);
/*kew*/

	/* check whether coefficient arrays have been precomputed */
	index = -1;
	iaux = pdeg_coarse*1000+pdeg_fine;
	for(i=0;i<nr_precomp;i++){
	  if(precomp_ind[i]==iaux) {
		index = i;
		break;
	  }
	}

	/* coefficient matrix has not been computed yet */
	if(index<0){

	  int is, iel_temp, pdeg_temp;
	  double dofs_temp[MAX_EL_DOFS];

	  precomp_ind[nr_precomp]=iaux;
	  index = nr_precomp;
	  nr_precomp++;

	  /* loop over all sons */
	  for(is=0;is<el_sons[0];is++){

		iel_temp=el_sons[is+1];

		dof_struct_temp_id = level_fine->l_dof_elem_to_struct[iel_fine];
		dof_struct_temp = level_fine->l_dof_struct[dof_struct_id];
		ibl_fine = dof_struct_fine.block_id;


		if(level_fine->pdeg_coarse == SIC_PDEG_FINEST){
		  pdeg_temp = pdr_get_ent_pdeg(problem_id, PDC_ELEMENT, iel_temp);
		  //pdeg_temp = apr_get_ent_pdeg(field_id, APC_ELEMENT, iel_temp);
		}
		else if(level_fine->pdeg_coarse == SIC_PDEG_COARSE_HIGH){
		  /* get sons pdegs and choose the highest */

		}
		else if(level_fine->pdeg_coarse == SIC_PDEG_COARSE_LOW){
		  /* get sons pdegs and choose the lowest */

		}
		else {
		  /* force pdeg_coarse to be the same for all elements */
		  pdeg_temp = level_fine->pdeg_coarse;
		}


		/* for all degrees of freedom of a coarse element */
		for(idofs=0;idofs<ndof_coarse;idofs++) {

		  /* coefficients will be computed for every coarse shape function */
		  for(jdofs=0;jdofs<ndof_coarse;jdofs++) dofs_temp[jdofs]=0.0;
		  dofs_temp[idofs]=1.0;

/* L2 project: to fine element (iel_fine) with dofs in shap_fun_precomp */
/*             from coarse element (iel_coarse) with dofs dofs_temp */
/* 	      i = -1; */
/* 	      apr_L2_proj(field_id,i, */
/* 			  iel_temp,&pdeg_temp, */
/* 			  shap_fun_precomp[index][is][idofs], */
/* 			  &iel_coarse,&pdeg_coarse,dofs_temp, NULL); */
		  pdr_L2_proj_sol(problem_id,
				  iel_temp,&pdeg_temp,
				  shap_fun_precomp[index][is][idofs],
				  &iel_coarse,&pdeg_coarse,dofs_temp);

		} /* end for all coarse shape functions: idofs */
	  } /* end for all sons: is */
	} /* end if coefficinets have not been precomputed yet */

	/* if projecting from coarser to finer mesh */
	if(mode==COARSE_TO_FINE){
/*kbw
#ifdef DEBUG_SIM
	printf("Projecting from iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		   iel_coarse, Ilev_coarse, ibl_coarse, posglob_coarse,
		   ndof_coarse);
	printf("Projecting to iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		   iel_fine, Ilev_fine, ibl_fine, posglob_fine,
		   ndof_fine);
	printf("index %d, ison %d\n", index, ison);
#endif
/*kew*/
	  double* array = &shap_fun_precomp[index][ison][0][0];
	  for(jdofs=0;jdofs<ndof_fine;jdofs++){
		Vec_to[posglob_fine+jdofs] = 0.0;
	  }
	  for(idofs=0;idofs<ndof_coarse;idofs++) {
		for(jdofs=0;jdofs<ndof_fine;jdofs++){

		  Vec_to[posglob_fine+jdofs] += Vec_from[posglob_coarse+idofs]
									  * array[idofs*MAX_EL_DOFS+jdofs];
		  //* shap_fun_precomp[index][ison][idofs][jdofs];

/*kbw
#ifdef DEBUG_SIM
		  printf("i %d, j %d, v_f %15.12lf, v_c %15.12lf, coef %15.12lf\n",
			 idofs,jdofs,
			 Vec_fine[posglob_fine+jdofs],
			 Vec_coarse[posglob_coarse+idofs],
			 shap_fun_precomp[index][ison][idofs][jdofs]);
#endif
/*kew*/

		}
	  }
	}
	/* if projecting from  finer to coarser mesh */
	else if(mode == FINE_TO_COARSE){

/*kbw
#ifdef DEBUG_SIM
	printf("Projecting from iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		   iel_fine, Ilev_fine, ibl_fine, posglob_fine,
		   ndof_fine);
	printf("Projecting to iel %d, level %d, ibl %d, posg %d, ndof %d\n",
		   iel_coarse, Ilev_coarse, ibl_coarse, posglob_coarse,
		   ndof_coarse);
	printf("index %d, ison %d\n", index, ison);
#endif
/*kew*/
	  double* array = &shap_fun_precomp[index][ison][0][0];
	  for(idofs=0;idofs<ndof_coarse;idofs++) {
		for(jdofs=0;jdofs<ndof_fine;jdofs++){

		  Vec_to[posglob_coarse+idofs] += Vec_from[posglob_fine+jdofs]
									  * array[idofs*MAX_EL_DOFS+jdofs];
		  //* shap_fun_precomp[index][ison][idofs][jdofs];

/*kbw
#ifdef DEBUG_SIM
		  printf("i %d, j %d, v_f %15.12lf, v_c %15.12lf, coef %15.12lf\n",
			 idofs,jdofs,
			 Vec_coarse[posglob_coarse+idofs],
			 Vec_fine[posglob_fine+jdofs],
			 shap_fun_precomp[index][ison][idofs][jdofs]);
#endif
/*kew*/

		}
	  }
/*kbw
getchar();
/*kew*/

	} /* end if projecting from finer to coarse */

	  } /* end loop over sons: ison */

	} /* end if element does not exist in both meshes */

  } // end if owned element

  } /* end loop over all elements of a coarse mesh */

  // projection must be followed by exchange of dofs for Level_to
  pdr_exchange_dofs(problem_id, Solver_id, Ilev_to, Vec_to);


  } // end if DG approximation

  return(0);

}

/*---------------------------------------------------------
  fem_vec_norm - to compute a norm of global vector (in parallel)
---------------------------------------------------------*/
double fem_vec_norm( /* returns: L2 norm of global Vector */
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  int Nrdof,            /* in: number of vector components */
  double* Vector        /* in: local part of global Vector */
  )
{

  int IONE=1; int problem_id; double vec_norm;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(siv_solver[Solver_id].parallel){

	problem_id = siv_solver[Solver_id].problem_id;
	vec_norm = pdr_vec_norm(problem_id, Solver_id, Level_id, Nrdof, Vector);

  } else {

	vec_norm = dnrm2_(&Nrdof, Vector, &IONE);

  }

  return(vec_norm);

}


/*---------------------------------------------------------
  fem_sc_prod - to compute a scalar product of two global vectors
---------------------------------------------------------*/
double fem_sc_prod( /* retruns: scalar product of Vector1 and Vector2 */
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  int Nrdof,           /* in: number of vector components */
  double* Vector1,     /* in: local part of global Vector */
  double* Vector2      /* in: local part of global Vector */
  )
{

  int IONE=1; int problem_id; double sc_prod;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(siv_solver[Solver_id].parallel){

	problem_id = siv_solver[Solver_id].problem_id;
	sc_prod = pdr_sc_prod(problem_id,Solver_id,Level_id,Nrdof,Vector1,Vector2);

  } else {

	sc_prod = ddot_(&Nrdof, Vector1, &IONE, Vector2, &IONE);

  }

  return(sc_prod);
}


/*---------------------------------------------------------
  fem_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
int fem_exchange_dofs(
  int Solver_id,        /* in: solver data structure to be used */
  int Subsystem_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  double* Vec_dofs  /* in: vector of dofs to be exchanged */
)
{

  int problem_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

//!!! If renumbering used - should be applied here !!!

  if(siv_solver[Solver_id].parallel){

	problem_id = siv_solver[Solver_id].problem_id;
	pdr_exchange_dofs(problem_id, Solver_id, Level_id, Vec_dofs);

  }

  return(1);
}
