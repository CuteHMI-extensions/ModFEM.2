
#ifndef _lah_bcrs_
#define _lah_bcrs_

#ifdef __cplusplus
extern "C" {
#endif
/*** CONSTANTS ***/

// assumption: stiffness matrix is composed of blocks with CONSTANT size
// BLS = const > 1
// TODO: passed from compilation options
#define BLS 4

#define ITC_MAX_MATRICES 20

/*** Data types ***/

typedef struct {

/* parameters */
  int Max_SM_size;      /* the maximal number of dofs in a stiffness matrix */
  int Nrblocks;		/* total number of small blocks */
  int Nrdofgl;		/* total number of degrees of freedom */
  int Nrbl_internal;		/* number of internal row block */
  int Nr_sm_blocks;
  int Half_bandwidth;
  int Max_nrneig; // maximal number of neighbors for a single block
  int Precon;  // preconditioner
  int ILU_k;   // number of neighborhood rings for ILU(k)
  int ILU_nr_sm_blocks;
  int ILU_half_bandwidth;
  
  double *rhs;   /* the global right hand side vector */
   
/* bcrs matrix */
  double *block_val;
  int *block_col_ind;
  int *block_row_ptr;
  
/* block diagonal preconditioner */  
  int *block_diag_ptr;
  double *block_diag_precon;
  int *block_diag_ips;

/* ILU decomposition of stiffness matrix */
  int *ilu_diag_ptr;
  double *ilu_val;
  int *ilu_col_ind;
  int *ilu_row_ptr;

  double *pivots;
  
} itt_bcrs_matrices;


/* GLOBAL VARIABLES */
 extern int   itv_nr_bcrs_matrices;        /* the number of solvers in the problem */
 extern int   itv_cur_bcrs_matrix_id;                /* ID of the current problem */
 extern itt_bcrs_matrices itv_bcrs_matrices[ITC_MAX_MATRICES];        /* array of solvers */


 
/*---------------------------------------------------------
  lar_allocate_SM_and_LV_bcrs - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
extern int lar_allocate_SM_and_LV_bcrs( // returns: matrix index in itv_matrices array
  int Nrdof_glob,  /* in: total number of DOFs */
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Block_size,  /* in: size of SM blocks (-1 - non-constant size) */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]<=0 - iblock is a ghost block with no rows associated
  int** L_offbl   /* in: list of lists of off diagonal blocks */
  );


/**--------------------------------------------------------
  lar_initialize_SM_and_LV_bcrs - to initialize stiffness matrix and/or load vector
---------------------------------------------------------*/
extern int lar_initialize_SM_and_LV_bcrs(
  int Matrix_id,   /* in: matrix ID */
  int Scope    /* in: indicator for the scope of computations: */
                   /*   LAC_SCOPE_SM_AND_LV */
                   /*   LAC_SCOPE_LV - do not touch SM! */
				    );

/**--------------------------------------------------------
  lar_get_storage_bcrs - to compute storage of SM, LV and preconditioner
---------------------------------------------------------*/
extern double lar_get_storage_bcrs( /* returns: storage in MB */
  int Matrix_id   /* in: matrix ID */
			       );

/**-----------------------------------------------------------
  lar_fill_assembly_table_int_ent_bcrs - to fill a part of the global assembly table
                              related to one integration entity, for which
                              lists of DOF blocks (their global positions) are provided
------------------------------------------------------------*/
extern int lar_fill_assembly_table_int_ent_bcrs( 
                         /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Nr_dof_bl,         /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  int *Assembly_table_int_ent /* part of the global assembly table */
  );


/*------------------------------------------------------------
  lar_assemble_SM_and_LV_with_table_bcrs - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local 
                           stiffness matrix and load vector and assembly table
------------------------------------------------------------*/
extern int lar_assemble_SM_and_LV_with_table_bcrs( 
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
  lar_assemble_SM_and_LV_bcrs - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local 
                           stiffness matrix and load vector
------------------------------------------------------------*/
extern int lar_assemble_SM_and_LV_bcrs( 
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
lar_allocate_preconditioner_bcrs - to allocate space for preconditioner 
---------------------------------------------------------*/
extern int lar_allocate_preconditioner_bcrs( /* returns:   >0 number of diagonal blocks */
                          /*	       <=0 - error */
  int Matrix_id,   /* in: matrix ID */
  int Precon,      /* in: type of preconditioner (lah_block.h line circa 45) */
  int ILU_k // in: for ILU(k) - k;   
  );

/**--------------------------------------------------------
  lar_fill_preconditioner_bcrs - to fill preconditioner
---------------------------------------------------------*/
extern int lar_fill_preconditioner_bcrs( 
  int Matrix_id   /* in: matrix ID */
	);

/**--------------------------------------------------------
  lar_free_preconditioner_bcrs - to free space for a block structure
---------------------------------------------------------*/
extern int lar_free_preconditioner_bcrs(
  int Matrix_id   /* in: matrix ID */
  );

/*---------------------------------------------------------
  lar_free_SM_and_LV_bcrs - to free space for a block structure
---------------------------------------------------------*/
extern int lar_free_SM_and_LV_bcrs(
  int Matrix_id   /* in: matrix ID */
  );

  
/*---------------------------------------------------------
lar_compute_residual_bcrs - to compute the residual of the system of equations,
	v = ( b - Ax ) (used also to compute the product v = -Ax)
---------------------------------------------------------*/
extern void lar_compute_residual_bcrs( 
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
lar_compute_preconditioned_residual_bcrs - to compute the residual of the  
	preconditioned system of equations, v = M^-1 * ( b - Ax )
        where M^-1 corresponds directly to the stored preconditioner matrix
---------------------------------------------------------*/
extern void lar_compute_preconditioned_residual_bcrs ( 
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
lar_perform_BJ_or_GS_iterations_bcrs - to perform one iteration of block Gauss-Seidel
	or block Jacobi algorithm:  v_out = v_in + M^-1 * ( b - A * v_in )
        where M^-1 results from stored preconditioner matrix and the algorithm
---------------------------------------------------------*/
extern void lar_perform_BJ_or_GS_iterations_bcrs(
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
lar_perform_rhsub_bcrs - to perform forward reduction and back-substitution for ILU
           preconditioning
---------------------------------------------------------*/
extern void lar_perform_rhsub_bcrs(
  int Matrix_id,   /* in: matrix ID */
  int Ndof,	   /* in: number of unknowns (components of v*) */ 
  double* V,	   /* in,out: vector of unknowns updated */
                   /* during the loop over subdomains */
  double* B	   /* in:  the rhs vector, if NULL take rhs */
                   /*      from block data structure */
	);


/*---------------------------------------------------------
  lar_get_SM_and_LV_crs_from_bcrs - convert bcrs matrix format to crs
---------------------------------------------------------*/
extern int lar_get_SM_and_LV_crs_from_bcrs( // returns: flag - 0 if not needed delete crs matrices
  int Matrix_id, /* in: matrix ID */
  int offset, /* in: offset in crs_row and crs_col matrix */
  int** crs_row,	   /* out: matrix of rows in crs */
  int** crs_col,	   /* out: matrix of column in crs */
  double** crs_val,	   /* out: matrix of value in crs */
  double** rhs	   /* out: rhs vector */
			   );



///////////// utility
extern int lar_util_chk_list_bcrs(	/* returns: */
			/* >0 - position on the list */
            		/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	);

#ifdef __cplusplus
}
#endif

#endif
