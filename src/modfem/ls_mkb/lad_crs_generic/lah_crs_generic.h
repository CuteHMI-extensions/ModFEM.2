/************************************************************************
lah_intf.h - header for matrix storage and operations package 

Assumption: Each matrix is stored with the corresponding preconditioner 


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
  lar_free_preconditioner - to free space of preconditioner structure
  lar_free_SM_and_LV - to free space of matrix structure

  lar_create_solver_structures_accel - utility to create data structures on GPU
  lar_get_crs_data - utility to get CRS parameters and pointers
  lar_compute_residual - to compute the residual of the not preconditioned 
	system of equations, v = ( b - Ax )
  lar_compute_preconditioned_residual_crs_generic - to compute the residual of the  
	preconditioned system of equations, v = M^-1 * ( b - Ax )
        where M^-1 corresponds directly to the stored preconditioner matrix
  lar_perform_BJ_or_GS_iterations - to perform one iteration of block Gauss- 
	Seidel or block Jacobi algorithm, v_out = v_in + M^-1 * ( b - A * v_in )
  lar_perform_rhsub - to perform forward reduction and back-substitution for ILU
           preconditioning

*************************************************************************/

#ifndef _lah_crs_generic_
#define _lah_crs_generic_

/*** External Parameters ***/
#define ITC_MAX_MATRICES 20

/*** Data types ***/

typedef struct {


/* parameters */
  int Max_SM_size;      /* the maximal number of dofs in a stiffness matrix */
  int Nrblocks;		/* total number of small blocks */
  int Nrdofgl;		/* total number of degrees of freedom */
  int Nrdof_internal;	/* internal number of degrees of freedom */
  int Nnz;              /* total number of non-zero elements */
  int Half_bandwidth;   /* FEM matrices are banded... */
  
  int *posg;

  int Precon;
  int ILU_k;   // number of neighborhood rings for ILU(k)
  
  double *rhs;   /* the global right hand side vector */

/* crs matrix */
  double *crs_val;
  int *crs_col_ind;
  int *crs_row_ptr;
  
  /* preconditioner */
  int *diag_ptr;  
  double *diag_precon;  
  
} itt_crs_generic_matrices;


/* GLOBAL VARIABLES */
 extern int   itv_nr_crs_generic_matrices;   /* the number of matrices managed by module */
 extern int   itv_cur_crs_generic_matrix_id; /* ID of the current matrix */
 extern itt_crs_generic_matrices itv_crs_generic_matrices[LAC_MAX_MATRICES];
                                             /* array of CRS matrices */


 
/*---------------------------------------------------------
  lar_allocate_SM_and_LV_crs_generic - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
extern int lar_allocate_SM_and_LV_crs_generic( // returns: matrix index in itv_crs_generic_matrices array
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]==0 - iblock is a ghost block with no rows associated with
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
  );


/**--------------------------------------------------------
  lar_initialize_SM_and_LV_crs_generic - to initialize stiffness matrix and/or load vector
---------------------------------------------------------*/
extern int lar_initialize_SM_and_LV_crs_generic(
  int Matrix_id,   /* in: matrix ID */
  int Scope    /* in: indicator for the scope of computations: */
                   /*   LAC_SCOPE_SM_AND_LV */
                   /*   LAC_SCOPE_LV - do not touch SM! */
				    );

/**--------------------------------------------------------
  lar_get_storage_crs_generic - to compute storage of SM, LV and preconditioner
---------------------------------------------------------*/
extern double lar_get_storage_crs_generic( /* returns: storage in MB */
  int Matrix_id   /* in: matrix ID */
			       );

/**-----------------------------------------------------------
  lar_fill_assembly_table_int_ent_crs_generic - to fill a part of the global assembly table
                              related to one integration entity, for which
                              lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
extern int lar_fill_assembly_table_int_ent_crs_generic( 
                         /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int *Assembly_table_int_ent /* part of the global assembly table */
  );


/*------------------------------------------------------------
  lar_assemble_SM_and_LV_with_table_crs_generic - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local 
                           stiffness matrix and load vector and assembly table
------------------------------------------------------------*/
extern int lar_assemble_SM_and_LV_with_table_crs_generic( 
                         /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Scope,    /* in: indicator for the scope of computations: */
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
  );

/*------------------------------------------------------------
  lar_assemble_SM_and_LV_crs_generic - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local 
                           stiffness matrix and load vector
------------------------------------------------------------*/
extern int lar_assemble_SM_and_LV_crs_generic( 
                         /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Scope,    /* in: indicator for the scope of computations: */
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
  );


/**--------------------------------------------------------
lar_allocate_preconditioner_crs_generic - to allocate space for preconditioner 
---------------------------------------------------------*/
extern int lar_allocate_preconditioner_crs_generic( /* returns:   >0 number of diagonal blocks */
                          /*	       <=0 - error */
  int Matrix_id,   /* in: matrix ID */
  int Precon,      /* in: type of preconditioner (lah_block.h line circa 45) */
  int ILU_k // in: for ILU(k) - k;  
  );

/**--------------------------------------------------------
  lar_fill_preconditioner_crs_generic - to fill preconditioner
---------------------------------------------------------*/
extern int lar_fill_preconditioner_crs_generic( 
  int Matrix_id   /* in: matrix ID */
	);

/**--------------------------------------------------------
  lar_free_preconditioner_crs_generic - to free space for a block structure
---------------------------------------------------------*/
extern int lar_free_preconditioner_crs_generic(
  int Matrix_id   /* in: matrix ID */
  );

/**---------------------------------------------------------
  lar_create_solver_structures_accel_crs_generic - utility to create data structures on GPU
---------------------------------------------------------*/
extern int lar_create_solver_structures_accel_crs_generic(
  int Matrix_id   /* in: matrix ID */
							  );

/**---------------------------------------------------------
lar_get_crs_data_crs_generic - utility to get CRS parameters and pointers
---------------------------------------------------------*/
int lar_get_crs_data_crs_generic( /* returns: >0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int* Nrdof_glob_p, /* output pointers: nrdof, nnz etc. - self explanatory */
  int* Nnz_p,
  int** Crs_col_ind_p,
  int** Crs_row_ptr_p,
  double** Crs_val_p,
  double** Rhs_p
		      );

/*---------------------------------------------------------
  lar_free_SM_and_LV_crs_generic - to free space for a block structure
---------------------------------------------------------*/
extern int lar_free_SM_and_LV_crs_generic(
  int Matrix_id   /* in: matrix ID */
  );

/*---------------------------------------------------------
lar_compute_residual_crs_generic - to compute the residual of the system of equations,
	v = ( b - Ax ) (used also to compute the product v = -Ax)
---------------------------------------------------------*/
extern void lar_compute_residual_crs_generic( 
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: indicator whether to use RHS */
  int Ini_zero,	/* in: flag for zero initial guess */ 
  int Ndof, 	/* in: number of unknowns (components of x) */
  double* X, 	/* in: input vector (may be NULL if Ini_zero==0) */
  double* B,	/* in:  the rhs vector, if NULL take rhs */
                /*      from block data structure ( B is not taken */
                /*      into account if Use_rhs!=1) */
  double* V 	/* out: v = b-Ax */
	);
  
/**--------------------------------------------------------
lar_compute_preconditioned_residual_crs_generic - to compute the residual of the  
	preconditioned system of equations, v = M^-1 * ( b - Ax )
        where M^-1 corresponds directly to the stored preconditioner matrix
---------------------------------------------------------*/
extern void lar_compute_preconditioned_residual_crs_generic ( 
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: indicator whether to use RHS */
  int Ini_zero,	/* in: flag for zero initial guess */ 
  int Ndof, 	/* in: number of unknowns (components of x) */
  double* X, 	/* in: initial guess vector */
  double* B,	/* in:  the rhs vector, if NULL take rhs */
                /*      from block data structure */
  double* V 	/* out: preconditioned residual, v = M^-1*(b-Ax) */
				   );

/**--------------------------------------------------------
lar_perform_BJ_or_GS_iterations_crs_generic - to perform one iteration of block Gauss-Seidel
	or block Jacobi algorithm:  v_out = v_in + M^-1 * ( b - A * v_in )
        where M^-1 results from stored preconditioner matrix and the algorithm
---------------------------------------------------------*/
extern void lar_perform_BJ_or_GS_iterations_crs_generic(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: 0 - no rhs, 1 - with rhs */
  int Ini_zero,	/* in: flag for zero initial guess */ 
  int Nr_prec,  /* in: number of preconditioner iterations */
  int Ndof,	/* in: number of unknowns (components of v*) */ 
  double* V,	/* in,out: vector of unknowns updated */
                /* during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
                /*      from block data structure */
	);


/**--------------------------------------------------------
lar_perform_rhsub_crs_generic - to perform forward reduction and back-substitution for ILU
           preconditioning
---------------------------------------------------------*/
extern void lar_perform_rhsub_crs_generic(
  int Matrix_id,   /* in: matrix ID */
  int Ndof,	   /* in: number of unknowns (components of v*) */ 
  double* V,	   /* in,out: vector of unknowns updated */
                   /* during the loop over subdomains */
  double* B	   /* in:  the rhs vector, if NULL take rhs */
                   /*      from block data structure */
	);

/*---------------------------------------------------------
  lar_get_SM_and_LV_crs_from_crs_generic - convert crs_generic matrix format to crs
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_crs_from_crs_generic( // returns: flag - 0 if not needed delete crs matrices
  int Matrix_id, /* in: matrix ID */
  int offset, /* in: offset in crs_row and crs_col matrix */
  int** crs_row,	   /* out: matrix of rows in crs */
  int** crs_col,	   /* out: matrix of column in crs */
  double** crs_val,	   /* out: matrix of value in crs */
  double** rhs	   /* out: rhs vector */
			   );

#endif
