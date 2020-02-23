/*************************************************************
File contains procedures:
 lar_allocate_SM_and_LV - to allocate space for stiffness matrix and load vector
  lar_initialize_SM_and_LV - to initialize stiffness matrix and/or load vector
  lar_get_storage - to compute storage of SM, LV and preconditioner
  lar_fill_assembly_table_int_ent - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
  lar_assemble_SM_and_LV_with_table - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
  lar_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
  lar_allocate_preconditioner - to allocate space for preconditioner
  lar_fill_preconditioner - to fill preconditioner
  lar_free_preconditioner - to free space for a block structure
  lar_free_SM_and_LV - to free space for a block structure

lar_create_solver_structures_accel - utility to create data structures on GPU
lar_get_crs_data - utility to get CRS parameters and pointers

lar_compute_residual - to compute the residual of the not preconditioned
	system of equations, v = ( b - Ax )
lar_compute_preconditioned_residual - to compute the residual of the
	preconditioned system of equations, v = M^-1 * ( b - Ax )
		where M^-1 corresponds directly to the stored preconditioner matrix
lar_perform_BJ_or_GS_iterations - to perform one iteration of block Gauss-Seidel
	or block Jacobi algorithm:  v_out = v_in + M^-1 * ( b - A * v_in )
lar_perform_rhsub - to perform forward reduction and back-substitution for ILU
		   preconditioning

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

#include <modfem/uth_log.h>

// interface for linear algebra package
#include <modfem/ls_mkb/lah_intf.h>

/* internal information for implementation modules */
#include "./lad_block/lah_block.h"
#include "./lad_crs/lah_crs.h"
#include "./lad_bcrs/lah_bcrs.h"
#include "./lad_crs_generic/lah_crs_generic.h"
#ifdef AMG_LAR_EXTENSIONS
#include "./amg_ext/amg_ext.h"
#include "./amg_ext/lah_petsc_interface.h"
#endif

/* GLOBAL VARIABLES */
int   lav_nr_matrices;        /* the number of solvers in the problem */
lat_matrices lav_matrices[LAC_MAX_MATRICES];        /* array of solvers */


/*---------------------------------------------------------
  lar_allocate_SM_and_LV - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
int lar_allocate_SM_and_LV( // returns: matrix index in lav_matrices array
  int Storage_type, /* in: enforced storage type; if LAC_STORAGE_UNDEFINED */
					/*     storage type is decided based on Block_size */
  int Solver_type,  /* in: solver id  ==0 - direct, >0 iterative */
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
  int Block_size,  /* in: size of SM blocks (-1 - non-constant size) */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  int** L_offbl   /* in: list of lists of off diagonal blocks */
	)
{


/*++++++++++++++++ executable statements ++++++++++++++++*/

  lav_nr_matrices++;

  if(lav_nr_matrices>=LAC_MAX_MATRICES){
	printf("Too much (%d) matrices requested in las_intf! (correct lah_intf.h)\n",
	   lav_nr_matrices);
  }

  //printf("\n\n\n%d\n\n\n",Block_size);

  int storage_type_temp;
  if(Storage_type!=LAC_STORAGE_UNDEFINED) storage_type_temp = Storage_type;
  else {

	// no LAC_STORAGE_BLOCK - block storage must be enforced...
	if(Block_size == 1) storage_type_temp = LAC_STORAGE_CRS;
	else if(Block_size > 1) storage_type_temp = LAC_STORAGE_BCRS;
	else storage_type_temp = LAC_STORAGE_CRS_GENERIC;

	//if(Solver_type==0 && Block_size > 1) storage_type_temp = LAC_STORAGE_CRS_GENERIC;

  }

  //!!!!!!!!!!! FOR DEBUGGING !!!!!!!!!!
  //storage_type_temp = LAC_STORAGE_BLOCK;
  //storage_type_temp = LAC_STORAGE_CRS_GENERIC;


  lav_matrices[lav_nr_matrices-1].Storage_type = storage_type_temp;

/*ok_kbw*/
  int i,j;
  printf("\nEntering lar_allocate_SM_and_LV: storage_type: requested %d, final %d\n",
	 Storage_type, storage_type_temp);
  printf("nrblocks %d, max_sm_size %d, nrdof_glob %d, block_size %d\n",
	 Nrblocks, Max_SM_size, Nrdof_glob, Block_size);
//*kew*/
/*kbw
// !!! offset 1 numbering of blocks moved to lad_block !!!
  printf("block_id,\tnrdofbl,\tposglob,\tnroffbl\tneighbors\n");
  for(i=0; i<Nrblocks; i++){
	printf("%d\t\t%d\t\t%d\t\t%d\t",i,Nrdofbl[i], Posglob[i], Nroffbl[i]);
	for(j=0;j<Nroffbl[i]; j++){
	  printf("%10d",L_offbl[i][j]);
	}
	printf("\n");
  }
//*kew*/

  int internal_matrix_id;
  if(storage_type_temp==LAC_STORAGE_CRS){

	assert(Nrdof_glob==Nrblocks);
	internal_matrix_id = lar_allocate_SM_and_LV_crs(Nrdof_glob, Max_SM_size,
							Nroffbl, L_offbl);
  }
  else if(storage_type_temp==LAC_STORAGE_BCRS){

	assert(Nrdof_glob==Nrblocks*Block_size);
	internal_matrix_id = lar_allocate_SM_and_LV_bcrs(Nrdof_glob, Max_SM_size,
							 Block_size,  Nroffbl, L_offbl);
  }
  else if(storage_type_temp==LAC_STORAGE_BLOCK){

	internal_matrix_id = lar_allocate_SM_and_LV_block(Nrdof_glob, Max_SM_size,
							  Nrblocks, Nrdofbl, Posglob,
							  Nroffbl, L_offbl);

  }
  #ifdef AMG_LAR_EXTENSIONS
  else if(storage_type_temp==LAC_STORAGE_PETSC){

		internal_matrix_id = lar_allocate_SM_and_LV_petsc(Max_SM_size, Nrblocks, Nrdof_glob,
								  Nrdofbl, Posglob,
								  Nroffbl, L_offbl);
		//TODO: fix
		int i;
		for(i = 0; i<20 ;i++){
			lav_matrices[i].Internal_matrix_id = 19 - i;
			lav_matrices[i].Storage_type = storage_type_temp;
		}
  }
  #endif
  else{

	// CRS GENERIC as universal storage
	internal_matrix_id = lar_allocate_SM_and_LV_crs_generic(Nrdof_glob, Max_SM_size,
								Nrblocks, Nrdofbl, Posglob,
								Nroffbl, L_offbl);

  }

  lav_matrices[lav_nr_matrices-1].Internal_matrix_id = internal_matrix_id;

  return(lav_nr_matrices-1);

}

/*---------------------------------------------------------
  lar_initialize_SM_and_LV - to initialize stiffness matrix and/or load vector
---------------------------------------------------------*/
int lar_initialize_SM_and_LV(
  int Matrix_id,   /* in: matrix ID */
  int Scope    /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  )

{

  int info=-1;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*ok_kbw*/
  printf("\nEntering lar_initialize_SM_and_LV: storage_type: %d, scope %d\n",
	 storage_type, Scope);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_initialize_SM_and_LV_crs(internal_matrix_id, Scope);
	break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  info = lar_initialize_SM_and_LV_crs_generic(internal_matrix_id, Scope);
	break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_initialize_SM_and_LV_bcrs(internal_matrix_id, Scope);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  info = lar_initialize_SM_and_LV_block(internal_matrix_id, Scope);
	  break;
	}
#ifdef AMG_LAR_EXTENSIONS
	case LAC_STORAGE_PETSC:
	{
	  info = lar_initialize_SM_and_LV_petsc(internal_matrix_id, Scope);
	  break;
	}
#endif
	default:
	{
	  info = lar_initialize_SM_and_LV_block(internal_matrix_id, Scope);
	  break;
	}
  }//switch
  return(info);
}

/*---------------------------------------------------------
  lar_get_storage - to compute storage of SM, LV and preconditioner
---------------------------------------------------------*/
double lar_get_storage( /* returns: storage in MB */
  int Matrix_id   /* in: matrix ID */
			   )
{

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;
  double storage;

/*ok_kbw*/
  printf("\nEntering lar_get_storage: storage_type: %d\n",
	 storage_type);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  storage = lar_get_storage_crs( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  storage = lar_get_storage_crs_generic( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  storage = lar_get_storage_bcrs( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  storage = lar_get_storage_block( internal_matrix_id );
	  break;
	}
	default:
	{
	  storage = lar_get_storage_block( internal_matrix_id );
	  break;
	}
  }//switch
  return(storage);
}

/**-----------------------------------------------------------
  lar_fill_assembly_table_int_ent - to fill a part of the global assembly table
							  related to one integration entity, for which
							  lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
int lar_fill_assembly_table_int_ent(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int *Assembly_table_int_ent /* part of the global assembly table */
					 ){

  int info=-1;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*kbw
  printf("\nEntering lar_fill_assembly_table_int_ent: storage_type: %d, number of dof blocks %d\n",
	 storage_type, Nr_dof_bl);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_fill_assembly_table_int_ent_crs( internal_matrix_id,
			  Nr_dof_bl,  L_bl_id, L_bl_nrdof,
			  Assembly_table_int_ent);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	//printf("\n\n!!!!!!!! assembly with table !!!!!!!!\n");

	  info = lar_fill_assembly_table_int_ent_crs_generic( internal_matrix_id,
			  Nr_dof_bl,  L_bl_id, L_bl_nrdof,
			  Assembly_table_int_ent);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_fill_assembly_table_int_ent_bcrs( internal_matrix_id,
			   Nr_dof_bl,  L_bl_id, L_bl_nrdof,
			   Assembly_table_int_ent);
	  break;
	}
	case LAC_STORAGE_BLOCK:{
	  info = 0;
	  break;
	}
	default:
	{
	  info = 0;
	  break;
	}
  }//switch

  return(info);
}


/*------------------------------------------------------------
  lar_assemble_SM_and_LV_with_table - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector and assembly table
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_with_table(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Scope,         /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  int* Local_to_global_int_ent, /* part of the global local_to_global table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
					   ){

  int info=-1;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*kbw
  printf("\nEntering lar_assemble_SM_and_LV_with_table: storage_type: %d, scope %d, number of dof blocks %d\n",
	 storage_type, Scope, Nr_dof_bl);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_assemble_SM_and_LV_with_table_crs( internal_matrix_id, Scope,
		Nr_dof_bl, Assembly_table_int_ent, Local_to_global_int_ent,
				Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  info = lar_assemble_SM_and_LV_with_table_crs_generic( internal_matrix_id, Scope,
				Nr_dof_bl, Assembly_table_int_ent,  Local_to_global_int_ent,
				Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_assemble_SM_and_LV_with_table_bcrs( internal_matrix_id, Scope,
				 Nr_dof_bl, Assembly_table_int_ent,  Local_to_global_int_ent,
				 Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  info = 0;
	  break;
	}
	default:
	{
	  info = 0;
	  break;
	}
  }//switch

  return(info);
}

/*------------------------------------------------------------
  lar_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
						   and the global load vector using the provided local
						   stiffness matrix and load vector
!!!!!!!!!******************************************!!!!!!!!!!!!!!
lar_assemble_SM_and_LV allows for sending GHOST blocks on lists of blocks.
GHOST blocks  (indicated by block->Lngb==NULL) are blocks with associated block
structures, but without associated rows of GLOBAL stiffness matrix.
Such blocks are used for matrix-vector product with rows belonging to ordinary
(not ghost) blocks. In lar_assemble_SM_and_LV we nevertheless assume that provided
Stiff_mat is square, i.e. there are entries for all blocks. Procedures
calling lar_assemble_SM_and_LV must conform to this convention
!!!!!!!!!******************************************!!!!!!!!!!!!!!
------------------------------------------------------------*/
int lar_assemble_SM_and_LV(
						 /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Scope,         /* in: indicator for the scope of computations: */
				   /*   LAC_SCOPE_SM_AND_LV */
				   /*   LAC_SCOPE_LV - do not touch SM! */
  int Nr_dof_bl,         /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
  )
{

  int info=-1;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*kbw
  printf("\nEntering lar_assemble_SM_and_LV: storage_type: %d, scope %d, number of dof blocks %d\n",
	 storage_type, Scope, Nr_dof_bl);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_assemble_SM_and_LV_crs( internal_matrix_id, Scope,
				 Nr_dof_bl, L_bl_id, L_bl_nrdof,
				 Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  info = lar_assemble_SM_and_LV_crs_generic( internal_matrix_id, Scope,
				 Nr_dof_bl, L_bl_id, L_bl_nrdof,
				 Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_assemble_SM_and_LV_bcrs( internal_matrix_id, Scope,
			Nr_dof_bl, L_bl_id, L_bl_nrdof,
			Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  info = lar_assemble_SM_and_LV_block( internal_matrix_id, Scope,
			 Nr_dof_bl, L_bl_id, L_bl_nrdof,
			 Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
	#ifdef AMG_LAR_EXTENSIONS
	case LAC_STORAGE_PETSC:
	{
		info = lar_assemble_SM_and_LV_petsc( internal_matrix_id, Scope,
			   Nr_dof_bl, L_bl_id, L_bl_nrdof,
			   Stiff_mat, Rhs_vect, Rewr_dofs);
		break;
	}
	#endif
	default:
	{

	  info = lar_assemble_SM_and_LV_block( internal_matrix_id, Scope,
			 Nr_dof_bl, L_bl_id, L_bl_nrdof,
			 Stiff_mat, Rhs_vect, Rewr_dofs);
	  break;
	}
  }//switch
  return(info);
}


/*---------------------------------------------------------
lar_allocate_preconditioner - to create preconditioner blocks corresponding
					   to small subdomains of neighboring elements
---------------------------------------------------------*/
int lar_allocate_preconditioner( /* returns:   >0 number of diagonal blocks */
						  /*	       <=0 - error */
  int Matrix_id,   /* in: matrix ID */
  int Precon,      /* in: type of preconditioner (lah_block.h line circa 45) */
  int Second_arg // in: for ILU(k) - k;
				 //     for all other: Block_type for preconditioner blocks */
  )
{

  int info=-1;

  int Nrblocks_dia=0;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*ok_kbw*/
  printf("\nEntering lar_allocate_preconditioner: storage_type: %d, precon %d, block_type/ILU_k %d\n",
	 storage_type, Precon, Second_arg);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_allocate_preconditioner_crs(internal_matrix_id, Precon, Second_arg);
	   break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  info = lar_allocate_preconditioner_crs_generic(internal_matrix_id, Precon, Second_arg);
	   break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_allocate_preconditioner_bcrs(internal_matrix_id, Precon, Second_arg);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  Nrblocks_dia = lar_allocate_preconditioner_block(internal_matrix_id, Precon, Second_arg);
	  break;
	}
	#ifdef AMG_LAR_EXTENSIONS
	case LAC_STORAGE_PETSC:
	{

		break;
	}
	#endif
	default:
	{
	  Nrblocks_dia = lar_allocate_preconditioner_block(internal_matrix_id, Precon, Second_arg);
	  break;
	}
  }//switch
  return(Nrblocks_dia);
}

/*---------------------------------------------------------
  lar_fill_preconditioner - to fill preconditioner (here by factorizing
							diagonal blocks or by block ILU(0) factorization)
---------------------------------------------------------*/
int lar_fill_preconditioner(
  int Matrix_id   /* in: matrix ID */
	)
{

  int info=-1;;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*ok_kbw*/
  printf("\nEntering lar_fill_preconditioner: storage_type: %d\n",
	 storage_type);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_fill_preconditioner_crs( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  info = lar_fill_preconditioner_crs_generic( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_fill_preconditioner_bcrs( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  info = lar_fill_preconditioner_block( internal_matrix_id );
	  break;
	}
	#ifdef AMG_LAR_EXTENSIONS
	case LAC_STORAGE_PETSC:
	{
	  info = lar_fill_preconditioner_petsc( internal_matrix_id );
	  break;
	}
	#endif
	default:
	{
	  info = lar_fill_preconditioner_block( internal_matrix_id );
	  break;
	}
  }//switch
  return(0);

}


/*---------------------------------------------------------
  lar_free_preconditioner - to free space for a block structure
---------------------------------------------------------*/
int lar_free_preconditioner(
  int Matrix_id   /* in: matrix ID */
  )
{
  int info;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*ok_kbw*/
  printf("\nEntering lar_free_preconditioner: storage_type: %d\n",
	 storage_type);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_free_preconditioner_crs( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  info = lar_free_preconditioner_crs_generic( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_free_preconditioner_bcrs( internal_matrix_id );
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  info = lar_free_preconditioner_block( internal_matrix_id );
	  break;
	}
	#ifdef AMG_LAR_EXTENSIONS
	case LAC_STORAGE_PETSC:
	{
	  info = lar_free_preconditioner_petsc( internal_matrix_id );
	  break;
	}
	#endif
	default:
	{
	  info = lar_free_preconditioner_block( internal_matrix_id );
	  break;
	}
  }//switch
  return(info);
}

/**---------------------------------------------------------
  lar_create_solver_structures_accel - utility to create data structures on GPU
---------------------------------------------------------*/
extern int lar_create_solver_structures_accel(
  int Matrix_id   /* in: matrix ID */
)
{
  int info;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;
  switch(storage_type){
	case LAC_STORAGE_CRS:
	  {
	info = lar_create_solver_structures_accel_crs(internal_matrix_id);
	break;
	  }
	case LAC_STORAGE_CRS_GENERIC:
	  {
	info = lar_create_solver_structures_accel_crs_generic(internal_matrix_id);
	break;
	  }
	default:
	  {
	printf("\nRequested CRS data for not-CSR storage type matrix! Exiting!\n");
	exit(0);
	break;
	  }
  }//switch

  return(info);
}

/*---------------------------------------------------------
lar_get_crs_data - utility to get CRS parameters and pointers
---------------------------------------------------------*/
int lar_get_crs_data( /* returns: >0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int* Nrdof_glob_p, /* output pointers: nrdof, nnz etc. - self explanatory */
  int* Nnz_p,
  int** Crs_col_ind_p,
  int** Crs_row_ptr_p,
  double** Crs_val_p,
  double** Rhs_p
			 )
{

  int info;

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;
  switch(storage_type){
	case LAC_STORAGE_CRS:
	  {
	info = lar_get_crs_data_crs(internal_matrix_id, Nrdof_glob_p, Nnz_p,
				Crs_col_ind_p, Crs_row_ptr_p, Crs_val_p, Rhs_p );
	break;
	  }
	case LAC_STORAGE_CRS_GENERIC:
	  {
	info = lar_get_crs_data_crs_generic(internal_matrix_id, Nrdof_glob_p, Nnz_p,
				Crs_col_ind_p, Crs_row_ptr_p, Crs_val_p, Rhs_p );
	break;
	  }
	default:
	  {
	printf("\nRequested CRS data for not-CSR storage type matrix! Exiting!\n");
	exit(0);
	break;
	  }
  }//switch

  return(info);
}


/*---------------------------------------------------------
  lar_free_SM_and_LV - to free space for a block structure
---------------------------------------------------------*/
int lar_free_SM_and_LV(
  int Matrix_id   /* in: matrix ID */
  )
{

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*ok_kbw*/
  printf("\nEntering lar_free_SM_and_LV: storage_type: %d\n",
	 storage_type);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  lar_free_SM_and_LV_crs(internal_matrix_id);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  lar_free_SM_and_LV_crs_generic(internal_matrix_id);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  lar_free_SM_and_LV_bcrs(internal_matrix_id);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  lar_free_SM_and_LV_block(internal_matrix_id);
	  break;
	}
	#ifdef AMG_LAR_EXTENSIONS
	case LAC_STORAGE_PETSC:
	{
	  lar_free_SM_and_LV_petsc(internal_matrix_id);
	  break;
	}
	#endif
	default:
	{
	  lar_free_SM_and_LV_block(internal_matrix_id);
	  break;
	}
  }//switch

  if(Matrix_id != lav_nr_matrices-1){

	printf("lar_free_SM_and_LV: requested deallocation of matrix with id %d - which is not the last (of %d matrices). ",
	  Matrix_id, lav_nr_matrices);
	exit(-1);

  }
  lav_matrices[Matrix_id].Internal_matrix_id = -1;
  lav_nr_matrices--;

  return(0);
}


/*---------------------------------------------------------
lar_compute_residual - to compute the residual of the system of equations,
	v = ( b - Ax ) (used also to compute the product v = -Ax)
---------------------------------------------------------*/
void lar_compute_residual (
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: indicator whether to use RHS */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Ndof, 	/* in: number of unknowns (components of x) */
  double* X, 	/* in: input vector (may be NULL if Ini_zero==0) */
  double* B,	/* in:  the rhs vector, if NULL take rhs */
				/*      from block data structure ( B is not taken */
				/*      into account if Use_rhs!=1) */
  double* V 	/* out: v = b-Ax */
	)
{

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*kbw
  printf("\nEntering lar_compute_residual: storage_type: %d, Use_rhs %d, Ini_zero %d, Ndof %d\n",
	 storage_type, Use_rhs, Ini_zero, Ndof);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  lar_compute_residual_crs(internal_matrix_id, Use_rhs, Ini_zero, Ndof, X, B, V);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  lar_compute_residual_crs_generic(internal_matrix_id, Use_rhs, Ini_zero, Ndof, X, B, V);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  lar_compute_residual_bcrs(internal_matrix_id, Use_rhs, Ini_zero, Ndof, X, B, V);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  lar_compute_residual_block(internal_matrix_id, Use_rhs, Ini_zero, Ndof, X, B, V);
	  break;
	}
	default:
	{
	  lar_compute_residual_block(internal_matrix_id, Use_rhs, Ini_zero, Ndof, X, B, V);
	  break;
	}
  }//switch

  return;
}



/*---------------------------------------------------------
lar_perform_BJ_or_GS_iterations - to perform one iteration of block Gauss-Seidel algorithm
	(block Jacobi switched on by it_matrix->Precon==BLOCK_JACOBI)
		v_out = v_in + M^-1 * ( b - A * v_in )
---------------------------------------------------------*/
void lar_perform_BJ_or_GS_iterations(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: 0 - no rhs, 1 - with rhs */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Nr_prec,  /* in: number of preconditioner iterations */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* in,out: vector of unknowns updated */
		/* during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	)
{

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*kbw
  printf("\nEntering lar_perform_BJ_or_GS_iterations: storage_type: %d, Use_rhs %d, Ini_zero %d, Ndof %d\n",
	 storage_type, Use_rhs, Ini_zero, Ndof);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  lar_perform_BJ_or_GS_iterations_crs(internal_matrix_id,
			Use_rhs, Ini_zero, Nr_prec, Ndof, V, B);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  lar_perform_BJ_or_GS_iterations_crs_generic(internal_matrix_id,
			Use_rhs, Ini_zero, Nr_prec, Ndof, V, B);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  lar_perform_BJ_or_GS_iterations_bcrs(internal_matrix_id,
			 Use_rhs, Ini_zero, Nr_prec, Ndof, V, B);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  lar_perform_BJ_or_GS_iterations_block(internal_matrix_id,
			  Use_rhs, Ini_zero, Nr_prec, Ndof, V, B);
	  break;
	}
	#ifdef AMG_LAR_EXTENSIONS
	case LAC_STORAGE_PETSC:
	{
		//TODO: change to lav id
		lar_perform_BJ_or_GS_iterations_petsc(Matrix_id,
				Use_rhs, Ini_zero, Nr_prec, Ndof, V, B);
		break;
	}
	#endif
	default:
	{
	  lar_perform_BJ_or_GS_iterations_block(internal_matrix_id,
			  Use_rhs, Ini_zero, Nr_prec, Ndof, V, B);
	  break;
	}
  }//switch

 return;
}

/*---------------------------------------------------------
lar_perform_rhsub - to perform forward reduction and back-substitution for ILU
		   preconditioning: v_out = M^-1 * b
---------------------------------------------------------*/
void lar_perform_rhsub(
  int Matrix_id,   /* in: matrix ID */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* out: vector of unknowns updated */
		/*      during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	)
{

  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

/*kbw
  printf("\nEntering lar_perform_rhsub: storage_type: %d, Ndof %d\n",
	 storage_type, Ndof);
  printf("matrix id: in interface %d, internally in implementation %d\n",
	 Matrix_id, internal_matrix_id);
/*kew*/

  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  lar_perform_rhsub_crs(internal_matrix_id, Ndof, V, B);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  lar_perform_rhsub_crs_generic(internal_matrix_id, Ndof, V, B);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  lar_perform_rhsub_bcrs(internal_matrix_id, Ndof, V, B);
	  break;
	}
	case LAC_STORAGE_BLOCK:
	{
	  lar_perform_rhsub_block(internal_matrix_id, Ndof, V, B);
	  break;
	}
	default:
	{
	  lar_perform_rhsub_block(internal_matrix_id, Ndof, V, B);
	  break;
	}
  }//switch

  return;
}

/*---------------------------------------------------------
  lar_get_SM_and_LV_crs - convert all matrix format to crs
---------------------------------------------------------*/
int lar_get_SM_and_LV_crs( // returns: flag - 0 if not needed delete crs matrices
  int Matrix_id, /* in: matrix ID */
  int offset, /* in: offset in crs_row and crs_col matrix */
  int** crs_row,	   /* out: matrix of rows in crs */
  int** crs_col,	   /* out: matrix of column in crs */
  double** crs_val,	   /* out: matrix of value in crs */
  double** rhs	   /* out: rhs vector */
	)
{
  int info=-1;
  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;

  int i;


  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_get_SM_and_LV_crs_from_crs(internal_matrix_id, offset,
									 crs_row, crs_col, crs_val, rhs);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{
	  info = lar_get_SM_and_LV_crs_from_crs_generic(internal_matrix_id, offset,
										   crs_row, crs_col, crs_val,rhs);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_get_SM_and_LV_crs_from_bcrs(internal_matrix_id, offset, crs_row, crs_col, crs_val,rhs);
	  break;
	}
	default:{
	  printf("\n Error 82345: Converting from %d to crs not implemented jet!!",storage_type);
	  exit(0);
	}
  }

  return info;
}

/*---------------------------------------------------------
  lar_get_SM_and_LV_coo - convert all matrix format to coo
---------------------------------------------------------*/
int lar_get_SM_and_LV_coo( // returns: flag - 0 if not needed delete coo matrices
  int Matrix_id, /* in: matrix ID */
  int* nnz,		/* out: Number of nonzero elements in matrix */
  int offset, /* in: offset in coo_row and coo_col matrix */
  int** coo_row,	   /* out: matrix of rows in coo */
  int** coo_col,	   /* out: matrix of column in coo */
  double** coo_val,	   /* out: matrix of value in coo */
  double** rhs	   /* out: rhs vector */
	)
{
  int info=-1;
  int internal_matrix_id = lav_matrices[Matrix_id].Internal_matrix_id;
  int storage_type = lav_matrices[Matrix_id].Storage_type;


  switch(storage_type){
	case LAC_STORAGE_CRS:
	{
	  info = lar_get_SM_and_LV_coo_from_crs(internal_matrix_id, nnz, offset,
									 coo_row, coo_col, coo_val, rhs);
	  break;
	}
	case LAC_STORAGE_CRS_GENERIC:
	{

	  info = lar_get_SM_and_LV_coo_from_crs_generic(internal_matrix_id, nnz, offset,
											  coo_row, coo_col, coo_val,rhs);
	  break;
	}
	case LAC_STORAGE_BCRS:
	{
	  info = lar_get_SM_and_LV_coo_from_bcrs(internal_matrix_id, nnz, offset,
										   coo_row, coo_col, coo_val,rhs);

	  /*
	 lar_get_SM_and_LV_crs_from_bcrs(internal_matrix_id,nnz, offset, crs_row, crs_col, crs_val,rhs);
			info = lar_get_SM_and_LV_coo_from_crs(internal_matrix_id, nnz, offset,
									 coo_row, coo_col, coo_val, rhs);
	  */

	  break;
	}
	default:{
	  printf("\n Error 82345: Converting from %d to coo not implemented jet!!",storage_type);
	  exit(0);
	}
  }

  return info;
}
